"""Bundle every manipulator spec into the interactive 3-D app (a single HTML file).

The DH rows and closed-form IK solutions are turned into JavaScript with sympy's
jscode printer, so the browser evaluates exactly the kinematics in the specs.
Forward kinematics, the Jacobian (finite differences) and the damped
pseudo-inverse are computed in the page itself.
"""

import argparse
import json
from pathlib import Path

import numpy as np
import sympy as sp

import re

import kinematics_lib as kin
from gen_walkthrough import load, shorthand, simplify_matrix, to_expr, to_float
from singularities import analyse as analyse_singularities

TEMPLATE = Path(__file__).with_name("app_template.html")


def js(expr) -> str:
    return sp.jscode(sp.sympify(expr))


def short_text(expr, joints) -> str:
    """Plain-text expression with the c12/s123 shorthand (the page subscripts the digits)."""
    short, _ = shorthand(sp.sympify(expr), joints)
    return re.sub(r"_\{(\w+)\}", r"\1", str(short))


def singular_data(J, joints, types, consts, q_ref) -> dict:
    sg = analyse_singularities(J, joints, types, consts, q_ref)
    names = [q.name for q in joints]
    families = []
    for it in sg["factors"]:
        f_js = f"({', '.join(names)}) => {js(it['f'].subs(consts))}"
        values = [{"joint": joints.index(it["joints"][0]), "value": float(v)} for v in (it["values"] or [])]
        families.append(f'{{"text": {json.dumps(short_text(it["f"], joints) + " = 0")}, '
                        f'"joints": {json.dumps([joints.index(q) for q in it["joints"]])}, '
                        f'"values": {json.dumps(values)}, "f": {f_js}}}')
    head = {"kind": sg["kind"], "idle": [joints.index(q) for q in sg["idle"]],
            "detText": short_text(sg["det"], joints) if sg["kind"] == "square" else None,
            "k": sg.get("k"), "by": sg.get("by"), "nMinors": len(sg.get("minors", []))}
    return json.dumps(head)[:-1] + f', "families": [{", ".join(families)}]}}'


def manipulator(spec_path: Path) -> str:
    spec, joints, transforms, types, consts = load(spec_path)
    names = [q.name for q in joints]

    rows = [[js(e.subs(consts)) for e in (t.theta, t.d, t.l, t.alpha)] for t in transforms]
    dh_fn = f"({', '.join(names)}) => [" + ", ".join("[" + ", ".join(r) + "]" for r in rows) + "]"

    pose = {**spec["sample"], **spec.get("schematic", {}).get("pose", {})}
    q_pose = [to_float(pose[n]) for n in names]
    q_sample = [to_float(spec["sample"][n]) for n in names]

    limits = []
    for q, qs, qp in zip(joints, q_sample, q_pose):
        if types[q] == "revolute":
            limits.append((-np.pi, np.pi))
        else:
            hi = 2 * max(abs(qs), abs(qp), 0.5)
            limits.append((-hi if min(qs, qp) < 0 else 0.0, hi))

    totals = [step.total_matrix for step in kin.kinematics_walkthrough(transforms)]
    fk = sp.lambdify(joints, [T.subs(consts) for T in totals], "numpy")

    # coordinates no joint can change (z for planar arms) are left out of the IK task space
    rng = np.random.default_rng(1)
    ends = np.array([np.asarray(fk(*[rng.uniform(lo, hi) for lo, hi in limits])[-1], float)[:3, 3]
                     for _ in range(12)])
    task = [k for k in range(3) if np.ptp(ends[:, k]) > 1e-9]

    def chain_length(q):
        frames = [np.eye(4)] + [np.asarray(T, float) for T in fk(*q)]
        total = 0.0
        for i in range(1, len(frames)):
            o0, o1, x = frames[i - 1][:3, 3], frames[i][:3, 3], frames[i][:3, 0]
            m = o1 - np.dot(o1 - o0, x) * x
            total += np.linalg.norm(m - o0) + np.linalg.norm(o1 - m)
        return total

    # reach (prismatic joints fully extended) bounds the target sliders; size (as posed) scales the scene
    reach = chain_length([hi if types[q] == "prismatic" else v for q, v, (lo, hi) in zip(joints, q_pose, limits)])
    size = chain_length(q_pose)

    # symbolic position Jacobian over the task coordinates, for the singularity analysis
    p = kin.kinematics(transforms)[:3, 3]
    J = simplify_matrix(sp.Matrix([p[k] for k in task]).jacobian(joints), joints)
    singular = singular_data(J, joints, types, consts, q_pose)

    ik_spec = spec.get("ik", {})
    inter = [(k, js(to_expr(v).subs(consts))) for k, v in ik_spec.get("intermediates", {}).items()]
    branches = []
    for i, sol in enumerate(ik_spec.get("solution", [])):
        lines = [f"const {k} = {e};" for k, e in inter]
        lines += [f"const {k} = {js(to_expr(v).subs(consts))};" for k, v in sol.items() if k != "name"]
        body = " ".join(lines) + f" return [{', '.join(names)}];"
        name = json.dumps(sol.get("name", f"Solution {i + 1}"), ensure_ascii=False)
        branches.append(f'{{"name": {name}, "solve": (x, y, z) => {{ {body} }}}}')

    data = {
        "id": spec_path.stem,
        "name": spec.get("name", spec_path.stem),
        "description": spec.get("description", ""),
        "joints": [{"name": n, "type": types[q], "min": lo, "max": hi}
                   for n, q, (lo, hi) in zip(names, joints, limits)],
        "pose": q_pose,
        "task": task,
        "reach": float(max(reach, 1e-3)),
        "size": float(max(size, 1e-3)),
        "z0": 0.0,
        "dhText": [[str(r[k]) for k in ("theta", "d", "a", "alpha")] for r in spec["dh"]],
    }
    return (json.dumps(data, ensure_ascii=False)[:-1]
            + f', "dh": {dh_fn}, "ik": [{", ".join(branches)}], "singular": {singular}}}')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-o", "--out", type=Path, required=True, help="output HTML file")
    ap.add_argument("--specs", type=Path, default=Path("manipulators"), help="directory with *.toml specs")
    ap.add_argument("--default", default="trr", help="manipulator shown first")
    args = ap.parse_args()

    specs = sorted(p for p in args.specs.glob("*.toml") if not p.name.startswith("_"))
    entries = [manipulator(p) for p in specs]
    html = TEMPLATE.read_text()
    html = html.replace("/*@@MANIPULATORS@@*/[]", "[\n" + ",\n".join(entries) + "\n]")
    html = html.replace("/*@@DEFAULT@@*/null", json.dumps(args.default))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(html)
    print(f"wrote {args.out} ({len(entries)} manipulators)")


if __name__ == "__main__":
    main()
