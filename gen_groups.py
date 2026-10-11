"""Generate the data for the group discussion (definition of a group, SO(3), SE(3)).

Every group axiom is checked symbolically with sympy on generally parameterised
elements, the worked examples are computed, and the forward kinematics of one
manipulator spec is shown to be an element of SE(3). Writes a JSON file plus an
SVG figure that groups.typ renders into an A4 document.
"""

import argparse
import json
from pathlib import Path

import numpy as np
import sympy as sp

import kinematics_lib as kin
from gen_walkthrough import (AXIS_COLORS, INK, INK_MUTED, Latexer, fmt, fmt_joint, load, plt,
                             simplify_matrix, to_float)


# ---------------------------------------------------------------- elements

def euler(k: str) -> sp.Matrix:
    """General rotation: ZYX Euler angles alpha_k, beta_k, gamma_k."""
    a, b, g = sp.symbols(f"alpha{k} beta{k} gamma{k}", real=True)
    return kin.rotation_z(a) * kin.rotation_y(b) * kin.rotation_x(g)


def vector(k: str) -> sp.Matrix:
    return sp.Matrix(sp.symbols(f"x{k} y{k} z{k}", real=True))


def se3(R: sp.Matrix, p: sp.Matrix) -> sp.Matrix:
    T = sp.eye(4)
    T[:3, :3] = R
    T[:3, 3] = p
    return T


def se3_inverse(T: sp.Matrix) -> sp.Matrix:
    R, p = T[:3, :3], T[:3, 3]
    return se3(R.T, -R.T * p)


# ---------------------------------------------------------------- checks

def vanishes(M: sp.Matrix) -> bool:
    """True when every entry is zero for all angles.

    sin/cos of each angle become symbols s, c and the entry is reduced modulo
    s**2 + c**2 - 1. Those relations form a Groebner basis, so the remainder is
    zero exactly when the entry is identically zero.
    """
    angles = sorted({f.args[0] for f in M.atoms(sp.sin, sp.cos)}, key=str)
    subs, relations, gens = {}, [], []
    for k, a in enumerate(angles):
        s, c = sp.symbols(f"s_{k} c_{k}")
        subs[sp.sin(a)], subs[sp.cos(a)] = s, c
        relations.append(s**2 + c**2 - 1)
        gens += [s, c]
    for e in M:
        e = sp.expand(e.xreplace(subs))
        if e != 0 and (not relations or sp.reduced(e, relations, *gens)[1] != 0):
            return False
    return True


def is_identity(M: sp.Matrix) -> bool:
    return vanishes(M - sp.eye(M.rows))


def is_rotation(R: sp.Matrix) -> bool:
    return is_identity(R.T * R) and vanishes(sp.Matrix([R.det() - 1]))


def is_se3(T: sp.Matrix) -> bool:
    return is_rotation(T[:3, :3]) and T[3, :] == sp.Matrix([[0, 0, 0, 1]])


def short_latex(M: sp.Matrix) -> str:
    """LaTeX with c_alpha, s_alpha for cos(alpha), sin(alpha)."""
    repl = {f: sp.Symbol(("c" if f.func is sp.cos else "s") + "_" + str(f.args[0]), real=True)
            for f in M.atoms(sp.sin, sp.cos)}
    return sp.latex(M.xreplace(repl))


def num_latex(M: np.ndarray, digits=4) -> str:
    """LaTeX of a numeric matrix, entries rounded (and -0 shown as 0)."""
    rows = [" & ".join(fmt(round(float(v), digits) + 0.0, digits) for v in row) for row in np.atleast_2d(M)]
    return r"\left[\begin{matrix}" + r"\\".join(rows) + r"\end{matrix}\right]"


def numeric_example(raw: sp.Matrix, joints, types, consts, sample) -> dict:
    """The same statements with numbers: T at the sample configuration, its inverse, their product."""
    q = [to_float(sample[j.name]) for j in joints]
    T = np.array(sp.lambdify(joints, raw.subs(consts), "numpy")(*q), dtype=float)
    R, p = T[:3, :3], T[:3, 3]
    T_inv = np.eye(4)
    T_inv[:3, :3] = R.T
    T_inv[:3, 3] = -R.T @ p
    orth_err = float(np.abs(R.T @ R - np.eye(3)).max())
    det = float(np.linalg.det(R))
    prod_err = float(max(np.abs(T @ T_inv - np.eye(4)).max(), np.abs(T_inv @ T - np.eye(4)).max()))
    inv_err = float(np.abs(T_inv - np.linalg.inv(T)).max())
    tol = 1e-12
    return {
        "q": [{"name": sp.latex(j), "value": fmt_joint(v, types[j])} for j, v in zip(joints, q)],
        "constants": [{"name": sp.latex(k), "value": fmt(v)} for k, v in consts.items()],
        "T": num_latex(T), "R": num_latex(R), "p": num_latex(p[:, None]),
        "RtR": num_latex(R.T @ R), "orth_err": f"{orth_err:.1e}", "det": fmt(round(det, 10)),
        "member": orth_err < tol and abs(det - 1) < tol and np.array_equal(T[3], [0, 0, 0, 1]),
        "Rt": num_latex(R.T), "mRtp": num_latex((-R.T @ p)[:, None]),
        "T_inv": num_latex(T_inv), "product": num_latex(T @ T_inv),
        "prod_err": f"{prod_err:.1e}", "inverse_ok": prod_err < tol,
        "inv_err": f"{inv_err:.1e}", "matches_general": inv_err < tol,
    }


# ---------------------------------------------------------------- figure

def plot_order(Rx: sp.Matrix, Rz: sp.Matrix, path: Path):
    """The same two rotations composed in both orders, next to the starting frame."""
    panels = [("start: $I$", sp.eye(3)), ("$R_x R_z$", Rx * Rz), ("$R_z R_x$", Rz * Rx)]
    fig = plt.figure(figsize=(7.5, 2.8))
    for k, (title, R) in enumerate(panels, start=1):
        ax = fig.add_subplot(1, 3, k, projection="3d")
        R = np.array(R, dtype=float)
        for col, axis in enumerate("xyz"):
            base = np.eye(3)[col]
            ax.plot(*zip(-base, base), color=INK_MUTED, linewidth=0.6, linestyle=":")
            ax.text(*(-1.2 * base), f"$-{axis}_0$", color=INK_MUTED, fontsize=7, ha="center", va="center")
            ax.quiver(0, 0, 0, *R[:, col], color=AXIS_COLORS[axis], linewidth=2, arrow_length_ratio=0.18)
            ax.text(*(1.2 * R[:, col] + 0.12 * base), f"${axis}$", color=AXIS_COLORS[axis], fontsize=9,
                    ha="center", va="center")
        ax.set_title(title, fontsize=10, color=INK, pad=0)
        ax.set(xlim=(-1, 1), ylim=(-1, 1), zlim=(-1, 1))
        ax.set_box_aspect((1, 1, 1))
        ax.view_init(elev=22, azim=35)
        ax.set_axis_off()
    fig.subplots_adjust(left=0, right=1, bottom=0, top=0.92, wspace=0)
    fig.savefig(path)
    plt.close(fig)


# ---------------------------------------------------------------- main

def build(spec_path: Path, out_dir: Path, rel_root: Path) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- SO(3): general elements and the axioms
    R, R1, R2, R3 = euler(""), euler("_1"), euler("_2"), euler("_3")
    Rx, Rz = kin.rotation_x(sp.pi / 2), kin.rotation_z(sp.pi / 2)
    so3 = {
        "R": short_latex(R),
        "checks": {
            "member": is_rotation(R),
            "closure": is_rotation(R1 * R2),
            "associative": vanishes((R1 * R2) * R3 - R1 * (R2 * R3)),
            "identity": is_rotation(sp.eye(3)) and vanishes(R * sp.eye(3) - R) and vanishes(sp.eye(3) * R - R),
            "inverse": is_rotation(R.T) and is_identity(R * R.T) and is_identity(R.T * R),
        },
        "Rx": sp.latex(Rx), "Rz": sp.latex(Rz),
        "RxRz": sp.latex(Rx * Rz), "RzRx": sp.latex(Rz * Rx),
        "commute": Rx * Rz == Rz * Rx,
        "sum": sp.latex(Rx + Rz), "sum_det": sp.latex((Rx + Rz).det()),
        "sum_is_rotation": is_rotation(Rx + Rz),
    }
    plot_order(Rx, Rz, out_dir / "order.svg")

    # --- SE(3): general elements and the axioms
    T, T1, T2, T3 = (se3(euler(k), vector(k)) for k in ("", "_1", "_2", "_3"))
    product = T1 * T2
    Ta = se3(Rz, sp.Matrix([1, 0, 0]))
    Tb = se3(Rx, sp.Matrix([0, 2, 0]))
    se3_doc = {
        "checks": {
            "member": is_se3(T),
            "closure": is_se3(product),
            "product_rule": vanishes(product - se3(T1[:3, :3] * T2[:3, :3], T1[:3, :3] * T2[:3, 3] + T1[:3, 3])),
            "associative": vanishes((T1 * T2) * T3 - T1 * (T2 * T3)),
            "identity": is_se3(sp.eye(4)) and vanishes(T * sp.eye(4) - T) and vanishes(sp.eye(4) * T - T),
            "inverse": is_se3(se3_inverse(T)) and is_identity(T * se3_inverse(T))
                       and is_identity(se3_inverse(T) * T),
        },
        "Ta": sp.latex(Ta), "Tb": sp.latex(Tb),
        "TaTb": sp.latex(Ta * Tb), "TbTa": sp.latex(Tb * Ta),
        "commute": Ta * Tb == Tb * Ta,
        "Ta_inv": sp.latex(se3_inverse(Ta)),
        "Ta_inv_ok": Ta * se3_inverse(Ta) == sp.eye(4),
    }

    # --- forward kinematics of a manipulator as an element of SE(3)
    spec, joints, transforms, types, consts = load(spec_path)
    lx = Latexer(joints)
    raw = kin.kinematics(transforms)
    shown = simplify_matrix(raw, joints)
    shown_inv = simplify_matrix(se3_inverse(shown), joints)
    n = len(joints)
    manipulator = {
        "name": spec.get("name", spec_path.stem),
        "n": n,
        "T": lx(shown),
        "T_inv": lx(shown_inv),
        "factors_ok": all(is_se3(t.transformation_matrix()) for t in transforms),
        "member": is_se3(raw),
        # the displayed (simplified) matrices use cos(q2 + q3), so expand them before reducing
        "inverse_ok": is_identity(raw * sp.expand_trig(shown_inv)) and vanishes(sp.expand_trig(shown) - raw),
        "legend": [{"short": k, "full": v} for k, v in sorted(lx.legend.items(), key=lambda kv: (len(kv[0]), kv[0]))],
        "numeric": numeric_example(raw, joints, types, consts, spec["sample"]),
    }

    return {
        "so3": so3,
        "se3": se3_doc,
        "manipulator": manipulator,
        "figures": {"order": (out_dir / "order.svg").resolve().relative_to(rel_root.resolve()).as_posix()},
    }


def failed_checks(data) -> list[str]:
    bad = [f"{g}: {k}" for g in ("so3", "se3") for k, ok in data[g]["checks"].items() if not ok]
    bad += [f"manipulator: {k}" for k in ("factors_ok", "member", "inverse_ok") if not data["manipulator"][k]]
    numeric = data["manipulator"]["numeric"]
    bad += [f"numeric example: {k}" for k in ("member", "inverse_ok", "matches_general") if not numeric[k]]
    return bad


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("spec", type=Path, help="manipulator TOML spec used as the SE(3) example")
    ap.add_argument("-o", "--out-dir", type=Path, required=True, help="directory for JSON and figures")
    ap.add_argument("--root", type=Path, default=Path("."), help="Typst project root (paths in JSON are relative to it)")
    args = ap.parse_args()

    data = build(args.spec, args.out_dir, args.root)
    (args.out_dir / "groups.json").write_text(json.dumps(data, indent=2, ensure_ascii=False))
    print(f"wrote {args.out_dir / 'groups.json'}")
    for name in failed_checks(data):
        print(f"CHECK FAILED: {name}")


if __name__ == "__main__":
    main()
