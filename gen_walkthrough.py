"""Generate a forward/inverse kinematics walkthrough for a manipulator spec.

Reads a TOML spec (DH table, constants, sample configuration, optional
closed-form IK), does the symbolic and numeric work, and writes a JSON file
plus SVG figures that walkthrough.typ renders into an A4 document.
"""

import argparse
import json
import re
import tomllib
from pathlib import Path

import matplotlib
import numpy as np
import sympy as sp

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import kinematics_lib as kin  # noqa: E402
from render3d import render_3d  # noqa: E402
from schematic import JointInfo, draw_schematic  # noqa: E402
from singularities import analyse as analyse_singularities  # noqa: E402

COORDS = [sp.Symbol(n, real=True) for n in ("x", "y", "z")]

# figure colours: robotics x/y/z convention, then categorical slots for arms
AXIS_COLORS = {"x": "#e34948", "y": "#008300", "z": "#2a78d6"}
ARM_COLORS = ["#2a78d6", "#eb6834", "#1baf7a"]
INK = "#0b0b0b"
INK_MUTED = "#52514e"
GRID = "#d9d8d4"


# ---------------------------------------------------------------- parsing

def to_expr(text) -> sp.Expr:
    """Parse a spec value; every free name becomes a real symbol, `deg` = pi/180."""
    expr = sp.sympify(str(text), locals={"deg": sp.pi / 180})
    return expr.subs({s: sp.Symbol(s.name, real=True) for s in expr.free_symbols})


def to_float(text) -> float:
    return float(to_expr(text))


# ---------------------------------------------------------------- symbolic

def shorthand(expr, joints):
    """Replace cos/sin of (sums of) joint variables with c_{12}, s_{123}, ..."""
    index = {}
    for q in joints:
        m = re.fullmatch(r"[A-Za-z_]+(\d)", q.name)
        if m is None:
            return expr, {}
        index[q] = m.group(1)

    legend = {}
    replacements = {}
    for f in expr.atoms(sp.cos, sp.sin):
        terms = sp.Add.make_args(f.args[0])
        if not all(t in index for t in terms):
            continue
        idx = "".join(sorted(index[t] for t in terms))
        name = ("c" if f.func is sp.cos else "s") + f"_{{{idx}}}"
        replacements[f] = sp.Symbol(name, real=True)
        legend[name] = sp.latex(f)
    return expr.xreplace(replacements), legend


class Latexer:
    """Turns expressions into LaTeX, collecting the shorthand legend on the way."""

    def __init__(self, joints):
        self.joints = joints
        self.legend = {}

    def __call__(self, expr) -> str:
        short, legend = shorthand(sp.sympify(expr), self.joints)
        self.legend.update(legend)
        return sp.latex(short).replace(r"\operatorname{atan}_{2}", r"\operatorname{atan2}")


def simplify_entry(e: sp.Expr, joints) -> sp.Expr:
    """simplify() alone leaves e.g. a2*(c1*c2 - s1*s2); trigsimp each link-parameter coefficient too."""
    e = sp.simplify(e)
    params = sorted(e.free_symbols - set(joints), key=lambda s: s.name)
    if params:
        collected = sp.collect(sp.expand(e), params, evaluate=False)
        e = sp.Add(*(sp.trigsimp(c) * k for k, c in collected.items()))
    return sp.factor_terms(e)


def simplify_matrix(m: sp.Matrix, joints) -> sp.Matrix:
    return m.applyfunc(lambda e: simplify_entry(e, joints))


def joint_types(transforms, joints):
    types = {}
    for q, t in zip(joints, transforms):
        types[q] = "prismatic" if t.d.has(q) else "revolute"
    return types


# ---------------------------------------------------------------- numeric

def fmt(v: float, digits=4) -> str:
    return f"{v:.{digits}f}".rstrip("0").rstrip(".") if abs(v) >= 1e-4 or v == 0 else f"{v:.2e}"


def fmt_joint(v: float, kind: str) -> str:
    if kind == "revolute":
        deg = round(float(np.degrees(np.arctan2(np.sin(v), np.cos(v)))), 2) + 0.0   # no "-0.00"
        if deg <= -179.995:
            deg = 180.0                                              # pi and -pi: show +180
        return f"{deg:.2f}°"
    return fmt(v)


def damped_least_squares(p_num, j_num, target, q0, lam=0.05, tol=1e-10, max_iter=200):
    q = np.array(q0, dtype=float)
    history = []
    for k in range(max_iter + 1):
        e = target - np.asarray(p_num(*q), dtype=float).ravel()
        history.append((k, float(np.linalg.norm(e))))
        if history[-1][1] < tol:
            break
        J = np.asarray(j_num(*q), dtype=float)
        q = q + J.T @ np.linalg.solve(J @ J.T + lam**2 * np.eye(len(e)), e)
    return q, history


def sequential_evaluator(steps, inputs):
    """Compile [(symbol, expr), ...] where each expr may use the inputs and earlier symbols."""
    compiled = []
    known = list(inputs)
    for sym, expr in steps:
        compiled.append((sym, sp.lambdify(known, expr, "numpy")))
        known.append(sym)

    def evaluate(values):
        env = list(values)
        out = {}
        with np.errstate(all="ignore"):
            for sym, f in compiled:
                v = complex(f(*env))
                v = v.real if abs(v.imag) < 1e-12 else np.nan
                env.append(v)
                out[sym] = v
        return out

    return evaluate


# ---------------------------------------------------------------- figures

def style_axes(ax):
    for side in ("top", "right"):
        if side in ax.spines:
            ax.spines[side].set_visible(False)
    ax.tick_params(colors=INK_MUTED, labelsize=8)
    ax.grid(color=GRID, linewidth=0.6)


def is_planar(frames):
    origins_z = [f[2, 3] for f in frames]
    z_axes = [f[:3, 2] for f in frames]
    return np.ptp(origins_z) < 1e-9 and all(abs(abs(z[2]) - 1) < 1e-9 for z in z_axes)


def plot_frames(frames, path: Path):
    origins = np.array([f[:3, 3] for f in frames])
    reach = max(np.linalg.norm(np.diff(origins, axis=0), axis=1).sum(), 1e-3)
    L = 0.12 * reach
    planar = is_planar(frames)

    fig = plt.figure(figsize=(6, 4.5))
    ax = fig.add_subplot(111) if planar else fig.add_subplot(111, projection="3d")
    dims = 2 if planar else 3
    ax.plot(*origins[:, :dims].T, color=INK_MUTED, linewidth=2, marker="o",
            markersize=8, markerfacecolor="white", markeredgecolor=INK_MUTED, zorder=1)

    for i, f in enumerate(frames):
        o = f[:3, 3]
        for axis, col in zip("xyz", range(3)):
            d = f[:3, col] * L
            if planar and axis == "z":
                ax.scatter(*o[:2], s=30, marker="o" if d[2] > 0 else "x",
                           color=AXIS_COLORS["z"], zorder=3)
                continue
            if planar:
                ax.annotate("", xy=o[:2] + d[:2], xytext=o[:2],
                            arrowprops=dict(arrowstyle="-|>", color=AXIS_COLORS[axis], lw=1.6), zorder=2)
            else:
                ax.quiver(*o, *d, color=AXIS_COLORS[axis], linewidth=1.6, arrow_length_ratio=0.25)
        label_at = o[:dims] - 0.06 * reach * (np.array([1, 1, 0])[:dims])
        ax.text(*label_at, f"{{{i}}}", color=INK, fontsize=9, ha="center", va="center")

    handles = [plt.Line2D([], [], color=c, lw=2, label=f"${a}_i$") for a, c in AXIS_COLORS.items()]
    ax.legend(handles=handles, frameon=False, fontsize=8, loc="best")
    _equal_limits(ax, origins, planar, pad=L * 1.5)
    style_axes(ax)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_ik(arms, reference, target, path: Path):
    """arms: list of (label, origins array); the last one is the numeric solution (dashed).
    reference: origins at the sample configuration, drawn thin and grey underneath."""
    planar = all(np.ptp(o[:, 2]) < 1e-9 for _, o in arms + [("", reference)])
    dims = 2 if planar else 3
    fig = plt.figure(figsize=(6, 4.5))
    ax = fig.add_subplot(111) if planar else fig.add_subplot(111, projection="3d")
    ax.plot(*reference[:, :dims].T, color=GRID, linewidth=5, solid_capstyle="round",
            label="sample configuration", zorder=0)
    for k, ((label, origins), color) in enumerate(zip(arms, ARM_COLORS)):
        numeric = k == len(arms) - 1
        ax.plot(*origins[:, :dims].T, color=color, linewidth=2, linestyle="--" if numeric else "-",
                marker="o", markersize=7, markerfacecolor="white", markeredgecolor=color, label=label)
    ax.scatter(*target[:dims], s=140, marker="*", color=INK, zorder=5, label="target")
    ax.legend(frameon=False, fontsize=8, loc="best")
    all_pts = np.vstack([o for _, o in arms] + [reference, target[None, :]])
    _equal_limits(ax, all_pts, planar, pad=0.1 * np.ptp(all_pts, axis=0).max())
    style_axes(ax)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def _equal_limits(ax, pts, planar, pad):
    lo, hi = pts.min(axis=0) - pad, pts.max(axis=0) + pad
    if planar:
        ax.set_aspect("equal")
        ax.set_xlim(lo[0], hi[0])
        ax.set_ylim(lo[1], hi[1])
        ax.set_xlabel("$x_0$")
        ax.set_ylabel("$y_0$")
    else:
        mid, half = (lo + hi) / 2, (hi - lo).max() / 2
        for setter, m in zip((ax.set_xlim, ax.set_ylim, ax.set_zlim), mid):
            setter(m - half, m + half)
        ax.set_box_aspect((1, 1, 1))
        ax.set_xlabel("$x_0$")
        ax.set_ylabel("$y_0$")
        ax.set_zlabel("$z_0$")


# ---------------------------------------------------------------- singularities

def singular_doc(sg, joints, types, lx, coords) -> dict:
    """Singularity analysis (singularities.analyse) in the form walkthrough.typ renders."""
    def vec(v):
        return "(" + ", ".join(fmt(float(x) + 0.0, 3) for x in v) + ")"   # + 0.0 turns -0.0 into 0.0

    def example(ex):
        q, u, v = ex
        return {"q": [{"name": sp.latex(j), "value": fmt_joint(float(val), types[j])} for j, val in zip(joints, q)],
                "lost": vec(u) if u is not None else None, "still": vec(v) if v is not None else None}

    def family(it):
        solution = None
        if it["solutions"] is not None:
            solution = sp.latex(it["joints"][0]) + r" \in " + sp.latex(it["solutions"])
        return {"cond": lx(it["f"]), "joints": [sp.latex(q) for q in it["joints"]], "solution": solution,
                "examples": [example(e) for e in it["examples"]]}

    out = {"kind": sg["kind"], "m": sg["m"], "n": sg["n"], "na": len(sg["active"]), "coords": coords,
           "joint_names": [sp.latex(j) for j in joints],
           "idle": [sp.latex(q) for q in sg["idle"]],
           "active": [sp.latex(joints[j]) for j in sg["active"]],
           "factors": [family(it) for it in sg["factors"]]}
    if sg["kind"] == "square":
        out["det"] = lx(sg["det"])
    else:
        out.update(k=sg["k"], by=sg["by"], minors=[lx(mn) for mn in sg["minors"]],
                   example=example(sg["example"]) if sg["example"] is not None else None)
    return out


# ---------------------------------------------------------------- schematic

def render_schematic(spec, transforms, joints, types, consts, path: Path, totals=None):
    """Kinematic schematic at [schematic] pose (default: the sample configuration)."""
    opts = spec.get("schematic", {})
    pose = {**spec["sample"], **opts.get("pose", {})}
    q = [to_float(pose[j.name]) for j in joints]
    if totals is None:
        totals = [step.total_matrix for step in kin.kinematics_walkthrough(transforms)]
    frames = sp.lambdify(joints, [sp.eye(4)] + [T.subs(consts) for T in totals], "numpy")(*q)

    def label(e):
        return sp.latex(e) if e != 0 else None

    info = [JointInfo(kind=types[j], label=sp.latex(j), theta_label=sp.latex(t.theta),
                      d_label=label(t.d), a_label=label(t.l))
            for j, t in zip(joints, transforms)]
    draw_schematic(frames, info, path, view=opts.get("view", "auto"))
    render_3d(frames, [types[j] for j in joints], [sp.latex(j) for j in joints], path.with_name("model3d.svg"))


def load(spec_path: Path):
    spec = tomllib.loads(spec_path.read_text())
    joints = [sp.Symbol(n, real=True) for n in spec["joints"]]
    transforms = [
        kin.DenavitHartenbergFrameTransform(
            theta=to_expr(r["theta"]), d=to_expr(r["d"]), l=to_expr(r["a"]), alpha=to_expr(r["alpha"]))
        for r in spec["dh"]
    ]
    assert len(transforms) == len(joints), "one DH row per joint expected"
    types = joint_types(transforms, joints)
    consts = {sp.Symbol(k, real=True): to_float(v) for k, v in spec.get("constants", {}).items()}
    return spec, joints, transforms, types, consts


# ---------------------------------------------------------------- main

def build(spec_path: Path, out_dir: Path, rel_root: Path) -> dict:
    spec, joints, transforms, types, consts = load(spec_path)
    n = len(joints)
    lx = Latexer(joints)

    # --- forward kinematics, step by step (kinematics_lib does the products)
    steps = []
    totals = []
    for i, step in enumerate(kin.kinematics_walkthrough(transforms), start=1):
        A = simplify_matrix(step.current_matrix, joints)
        T = simplify_matrix(step.total_matrix, joints)
        totals.append(T)
        t = step.current_transform
        steps.append({
            "i": i,
            "row": [lx(t.theta), lx(t.d), lx(t.l), lx(t.alpha)],
            "A": lx(A),
            "T": lx(T),
        })

    T_n = totals[-1]
    p = T_n[:3, 3]
    task_rows = [k for k in range(3) if any(p[k].has(q) for q in joints)]
    p_task = sp.Matrix([p[k] for k in task_rows])
    task_coords = [COORDS[k] for k in task_rows]
    fixed_coords = [{"name": sp.latex(COORDS[k]), "value": lx(p[k])} for k in range(3) if k not in task_rows]
    m = len(task_rows)

    J = simplify_matrix(p_task.jacobian(joints), joints)

    # --- numeric model with constants substituted
    frames_num = sp.lambdify(joints, [sp.eye(4)] + [T.subs(consts) for T in totals], "numpy")
    p_full_num = sp.lambdify(joints, p.subs(consts), "numpy")
    p_num = sp.lambdify(joints, p_task.subs(consts), "numpy")
    j_num = sp.lambdify(joints, J.subs(consts), "numpy")

    def frames_at(q):
        return [np.asarray(f, dtype=float) for f in frames_num(*q)]

    def origins_at(q):
        return np.array([f[:3, 3] for f in frames_at(q)])

    q_sample = np.array([to_float(spec["sample"][q.name]) for q in joints])
    singular = singular_doc(analyse_singularities(J, joints, types, consts, q_sample), joints, types, lx,
                            [sp.latex(c) for c in task_coords])
    out_dir.mkdir(parents=True, exist_ok=True)
    plot_frames(frames_at(q_sample), out_dir / "frames.svg")
    render_schematic(spec, transforms, joints, types, consts, out_dir / "schematic.svg", totals)
    target_full = np.asarray(p_full_num(*q_sample), dtype=float).ravel()
    target = target_full[task_rows]

    ik_spec = spec.get("ik", {})
    arms = []

    # --- closed-form IK (optional): typeset, evaluate at target, verify on random targets
    closed = None
    if "solution" in ik_spec:
        intermediates = [(sp.Symbol(k, real=True), to_expr(v).subs(consts))
                         for k, v in ik_spec.get("intermediates", {}).items()]
        inter_tex = [{"name": sp.latex(to_expr(k)), "expr": lx(to_expr(v))}
                     for k, v in ik_spec.get("intermediates", {}).items()]
        rng = np.random.default_rng(0)
        solutions = []
        for sol in ik_spec["solution"]:
            assign = [(sp.Symbol(k, real=True), to_expr(v).subs(consts)) for k, v in sol.items() if k != "name"]
            missing = {q.name for q in joints} - {s.name for s, _ in assign}
            assert not missing, f"solution '{sol.get('name')}' does not assign {missing}"
            evaluate = sequential_evaluator(intermediates + assign, COORDS)

            def solve_for(xyz):
                vals = evaluate(xyz)
                return np.array([vals[sp.Symbol(q.name, real=True)] for q in joints])

            q_sol = solve_for(target_full)
            err = float(np.linalg.norm(np.asarray(p_num(*q_sol), dtype=float).ravel() - target))

            worst, failed, trials = 0.0, 0, 500
            for _ in range(trials):
                q_rand = np.array([rng.uniform(-np.pi, np.pi) if types[q] == "revolute"
                                   else rng.uniform(0, 2 * max(abs(v), 0.5))
                                   for q, v in zip(joints, q_sample)])
                xyz = np.asarray(p_full_num(*q_rand), dtype=float).ravel()
                q_back = solve_for(xyz)
                if np.isnan(q_back).any():
                    failed += 1
                    continue
                e = np.linalg.norm(np.asarray(p_num(*q_back), dtype=float).ravel() - xyz[task_rows])
                worst = max(worst, float(e))

            name = sol.get("name", f"Solution {len(solutions) + 1}")
            solutions.append({
                "name": name,
                "assign": [{"name": sp.latex(s), "expr": lx(to_expr(sol[s.name]))} for s, _ in assign],
                "at_target": [{"name": sp.latex(q), "value": fmt_joint(v, types[q])} for q, v in zip(joints, q_sol)],
                "error": fmt(err),
                "verify": {"trials": trials, "failed": failed, "max_error": f"{worst:.1e}",
                           "ok": worst < 1e-8 and failed == 0},
            })
            arms.append((name, origins_at(q_sol)))
        closed = {"intermediates": inter_tex, "solutions": solutions}

    # --- numeric IK (always)
    q0 = np.array([to_float(ik_spec.get("initial_guess", {}).get(q.name, 0.5)) for q in joints])
    lam = float(ik_spec.get("damping", 0.05))
    q_num, history = damped_least_squares(p_num, j_num, target, q0, lam=lam)
    if len(arms) >= len(ARM_COLORS):
        arms = arms[:len(ARM_COLORS) - 1]
    arms.append(("numeric (DLS)", origins_at(q_num)))
    plot_ik(arms, origins_at(q_sample), target_full, out_dir / "ik.svg")

    shown = history if len(history) <= 16 else history[:12] + [None] + history[-3:]
    numeric = {
        "damping": fmt(lam),
        "q0": [{"name": sp.latex(q), "value": fmt_joint(v, types[q])} for q, v in zip(joints, q0)],
        "iterations": [None if h is None else {"k": h[0], "err": f"{h[1]:.3e}"} for h in shown],
        "converged": history[-1][1] < 1e-10,
        "result": [{"name": sp.latex(q), "value": fmt_joint(v, types[q])} for q, v in zip(joints, q_num)],
    }

    # general DH matrix for the formula section
    th, d, a, al = sp.symbols("theta_i d_i a_i alpha_i", real=True)
    A_general = kin.DenavitHartenbergFrameTransform(theta=th, d=d, l=a, alpha=al).transformation_matrix()

    def rel(path: Path) -> str:
        return path.resolve().relative_to(rel_root.resolve()).as_posix()

    figure = spec.get("figure")
    return {
        "name": spec.get("name", spec_path.stem),
        "description": spec.get("description", ""),
        "figure": rel(spec_path.parent / figure) if figure else None,
        "joints": [{"name": sp.latex(q), "type": types[q]} for q in joints],
        "constants": [{"name": sp.latex(k), "value": fmt(v)} for k, v in consts.items()],
        "A_general": sp.latex(A_general),
        "steps": steps,
        "position": {"coords": sp.latex(sp.Matrix(COORDS)), "expr": lx(p)},
        "task_space": {"m": m, "n": n, "coords": [sp.latex(c) for c in task_coords], "fixed": fixed_coords},
        "jacobian": {"p": lx(p_task), "q": sp.latex(sp.Matrix(joints)), "J": lx(J)},
        "singular": singular,
        "sample": {
            "q": [{"name": sp.latex(q), "value": fmt_joint(v, types[q])} for q, v in zip(joints, q_sample)],
            "target": [{"name": sp.latex(c), "value": fmt(v)} for c, v in zip(COORDS, target_full)],
        },
        "figures": {"frames": rel(out_dir / "frames.svg"), "ik": rel(out_dir / "ik.svg"),
                    "schematic": rel(out_dir / "schematic.svg"), "model3d": rel(out_dir / "model3d.svg")},
        "ik": {"closed": closed, "numeric": numeric},
        "legend": [{"short": k, "full": v} for k, v in sorted(lx.legend.items(), key=lambda kv: (len(kv[0]), kv[0]))],
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("spec", type=Path, help="manipulator TOML spec")
    ap.add_argument("-o", "--out-dir", type=Path, required=True, help="directory for JSON and figures")
    ap.add_argument("--root", type=Path, default=Path("."), help="Typst project root (paths in JSON are relative to it)")
    ap.add_argument("--schematic-only", action="store_true", help="only draw the schematic (fast, no simplification)")
    args = ap.parse_args()

    if args.schematic_only:
        spec, joints, transforms, types, consts = load(args.spec)
        args.out_dir.mkdir(parents=True, exist_ok=True)
        render_schematic(spec, transforms, joints, types, consts, args.out_dir / "schematic.svg")
        print(f"wrote {args.out_dir / 'schematic.svg'}")
        return

    data = build(args.spec, args.out_dir, args.root)
    (args.out_dir / "walkthrough.json").write_text(json.dumps(data, indent=2, ensure_ascii=False))
    print(f"wrote {args.out_dir / 'walkthrough.json'}")


if __name__ == "__main__":
    main()
