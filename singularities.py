"""Singular configurations of the position Jacobian J = dp/dq.

Square J (as many joints as task coordinates): det J is factored, and every
factor that depends on the joints describes one family of singular
configurations. A factor in a single joint is solved for that joint; a factor
coupling several joints is a condition, solved numerically for an example.

Joints whose column of J is identically zero never move the end effector
(e.g. a joint rotating the last link about its own axis); they are reported and
left out. If the remaining J is not square (redundant: more joints than
coordinates, or the other way round) it loses rank only when every k x k minor
vanishes, k = min(m, n) (Cauchy-Binet). The minors are listed and an example
singular pose is found numerically.

For each example pose the smallest singular value of J gives the task direction
the end effector cannot move in (left singular vector) and, for square J, the
joint motion that leaves it still (right singular vector).
"""

from itertools import combinations

import numpy as np
import sympy as sp


def _wrap(v: float) -> float:
    return float(np.arctan2(np.sin(v), np.cos(v)))


def _direction(v: np.ndarray) -> np.ndarray:
    """Unit vector with tiny components zeroed and the first non-zero component positive."""
    v = np.where(np.abs(v) < 1e-9, 0.0, v)
    v = v / np.linalg.norm(v)
    nz = np.flatnonzero(np.abs(v) > 1e-9)
    return -v if nz.size and v[nz[0]] < 0 else v


def _lost_directions(J_num, q, k, want_u, want_v):
    """At a pose where rank J dropped below k: u is the task direction the effector cannot move in,
    v the motion of the analysed joints that leaves it still. Each only when it is unique."""
    U, S, Vt = np.linalg.svd(np.asarray(J_num(*q), dtype=float))
    u = _direction(U[:, k - 1]) if want_u else None
    v = _direction(Vt[k - 1]) if want_v else None
    return u, v


def _single_joint_roots(f, q, kind, consts):
    """Values of joint q (revolute: in (-pi, pi]) where f = 0, symbolic and numeric; None if not a finite set."""
    domain = sp.Interval(-sp.pi, sp.pi) if kind == "revolute" else sp.S.Reals
    try:
        sym = sp.solveset(f, q, domain)
        num = sp.solveset(f.subs(consts), q, domain)
    except (NotImplementedError, ValueError, TypeError):
        return None, None
    if not (isinstance(sym, sp.FiniteSet) and isinstance(num, sp.FiniteSet)):
        return None, None
    if kind == "revolute":
        sym = sp.FiniteSet(*[s for s in sym if s != -sp.pi])        # -pi and pi are the same angle
        unique = {}
        for w in (_wrap(float(s)) for s in num):
            w = np.pi if abs(abs(w) - np.pi) < 1e-9 else w      # -pi and pi: keep +pi
            unique.setdefault(round(w, 9), w)
        values = sorted(unique.values())
    else:
        values = sorted(float(s) for s in num)
    return sym, values


def _newton_one_joint(f_num, grads, q0, k):
    q = np.array(q0, dtype=float)
    for _ in range(60):
        v, d = float(f_num(*q)), float(grads[k](*q))
        if abs(v) < 1e-13 or abs(d) < 1e-14:
            break
        q[k] -= v / d
    return q if abs(float(f_num(*q))) < 1e-9 else None


def _volume(J, k):
    """Product of the k largest singular values: zero exactly when rank J < k."""
    return float(np.prod(np.linalg.svd(J, compute_uv=False)[:k]))


def _descend_to_singularity(J_num, q0, kinds, limits, k, steps=400):
    """Gradient descent on the product of the k largest singular values of J, starting from q0."""
    def w(q):
        return _volume(np.asarray(J_num(*q), dtype=float), k)

    q = np.array(q0, dtype=float)
    t = 0.3
    for _ in range(steps):
        f0 = w(q)
        if f0 < 1e-10:
            break
        g = np.array([(w(q + h) - w(q - h)) / 2e-6 for h in np.eye(len(q)) * 1e-6])
        if np.linalg.norm(g) < 1e-14:
            break
        d = -g / np.linalg.norm(g)
        while t > 1e-12:
            qn = np.clip(q + t * d, [lo for lo, _ in limits], [hi for _, hi in limits])
            if w(qn) < f0:
                q, t = qn, min(2 * t, 0.5)
                break
            t /= 2
        else:
            break
    return np.array([_wrap(v) if k == "revolute" else v for v, k in zip(q, kinds)])


def _undouble(expr):
    """sin(2x) -> 2 sin(x) cos(x), so the factors stay single-angle (s3*c3 rather than sin(2 q3)/2)."""
    return expr.replace(lambda e: isinstance(e, sp.sin) and (e.args[0] / 2).is_Symbol,
                        lambda e: 2 * sp.sin(e.args[0] / 2) * sp.cos(e.args[0] / 2))


def analyse(J, joints, types, consts, q_ref):
    """J: simplified symbolic position Jacobian (m x n); q_ref: numeric pose to build examples around."""
    m, n = J.shape
    kinds = [types[q] for q in joints]
    q_ref = np.asarray(q_ref, dtype=float)
    limits = [(-np.pi, np.pi) if k == "revolute" else (-1e9, 1e9) for k in kinds]

    active = [j for j in range(n) if any(sp.simplify(J[i, j]) != 0 for i in range(m))]
    idle = [joints[j] for j in range(n) if j not in active]
    Ja = J[:, active]
    ja = [joints[j] for j in active]
    na = len(active)
    Ja_num = sp.lambdify(ja, Ja.subs(consts), "numpy")

    def full(q_active, base):
        q = np.array(base, dtype=float)
        q[active] = q_active
        return q

    result = {"m": m, "n": n, "idle": idle, "active": active}
    square = m == na
    want_u = na >= m      # with fewer moving joints than coordinates, several directions are lost
    want_v = na <= m      # a redundant arm can always move its joints without moving the effector
    k_rank = min(m, na)
    q_ref_a = q_ref[active]

    def example(qe):
        u, v = _lost_directions(Ja_num, qe, k_rank, want_u, want_v)
        if v is not None:
            v = full(v, np.zeros(n))      # joint motion over all joints; idle joints get 0
        return full(qe, q_ref), u, v

    def factor_items(expr):
        """One item per joint-dependent factor of expr, with solutions and example poses."""
        items = []
        for f, _ in sp.factor_list(expr)[1]:
            involved = [q for q in ja if f.has(q)]
            if not involved:
                continue
            item = {"f": f, "joints": involved, "solutions": None, "values": None, "examples": []}
            if len(involved) == 1:
                q = involved[0]
                item["solutions"], item["values"] = _single_joint_roots(f, q, types[q], consts)
                for value in item["values"] or []:
                    qe = q_ref_a.copy()
                    qe[ja.index(q)] = value
                    item["examples"].append(qe)
            if not item["examples"]:
                f_num = sp.lambdify(ja, f.subs(consts), "numpy")
                grads = [sp.lambdify(ja, sp.diff(f.subs(consts), q), "numpy") for q in ja]
                k = max(range(na), key=lambda i: abs(float(grads[i](*q_ref_a))) if ja[i] in involved else -1)
                qe = _newton_one_joint(f_num, grads, q_ref_a, k)
                if qe is not None:
                    item["examples"].append(np.array([_wrap(v) if types[q] == "revolute" else v for v, q in zip(qe, ja)]))
            item["examples"] = [example(qe) for qe in item["examples"]]
            items.append(item)
        return items

    if square:
        det = sp.factor(_undouble(sp.simplify(Ja.det())))
        coeff, factor_list = sp.factor_list(det)
        constant = coeff * sp.Mul(*[f ** e for f, e in factor_list if not f.has(*joints)])
        result.update(kind="square", det=det, constant=constant, factors=factor_items(det))
        return result

    if na > m:   # redundant: square submatrices over column subsets
        minors = [sp.factor(_undouble(sp.simplify(Ja[:, list(c)].det()))) for c in combinations(range(na), m)]
    else:        # fewer moving joints than coordinates: row subsets
        minors = [sp.factor(_undouble(sp.simplify(Ja[list(r), :].det()))) for r in combinations(range(m), na)]
    minors = [mn for mn in minors if mn != 0]
    # a factor shared by every minor makes all of them vanish at once: a singular family of its own
    common = minors[0] if minors else sp.Integer(0)
    for mn in minors[1:]:
        common = sp.gcd(common, mn)
    kinds_a = [types[q] for q in ja]
    qe = _descend_to_singularity(Ja_num, q_ref_a, kinds_a, [limits[j] for j in active], k_rank)
    S = np.linalg.svd(np.asarray(Ja_num(*qe), dtype=float), compute_uv=False)
    found = example(qe) if S[k_rank - 1] < 1e-6 * max(S[0], 1e-12) else None
    result.update(kind="minors", k=k_rank, minors=minors, by="columns" if na > m else "rows",
                  factors=factor_items(common) if minors else [], example=found)
    return result
