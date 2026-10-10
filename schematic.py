"""Textbook-style kinematic schematic drawn from DH frames.

Geometry comes straight from the frames: joint i sits at origin o_{i-1} with its
axis along z_{i-1}, and link i runs o_{i-1} -> m_i (along z_{i-1}, the d_i offset)
-> o_i (along x_i, the a_i length).
"""

from dataclasses import dataclass
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Arc, Circle, Ellipse, Polygon  # noqa: E402

INK = "#0b0b0b"
INK_MUTED = "#52514e"

plt.rcParams["mathtext.fontset"] = "cm"


@dataclass
class JointInfo:
    kind: str          # "revolute" (bearing ]|[ or circle) | "prismatic" (crossed box)
    label: str         # mathtext for the joint variable, e.g. "q_{1}"
    theta_label: str   # mathtext for theta_i (used on angle arcs)
    d_label: str | None
    a_label: str | None


# ---------------------------------------------------------------- projections

def oblique(p):
    """Cabinet-style view as in the course sheets: y0 right, z0 up, x0 toward the viewer (down-left)."""
    k, phi = 0.5, np.radians(45)
    return np.array([p[1] - k * p[0] * np.cos(phi), p[2] - k * p[0] * np.sin(phi)])


def top(p):
    return np.array([p[0], p[1]])


VIEWS = {"oblique": oblique, "top": top}


# ---------------------------------------------------------------- primitives

def _unit(v):
    n = np.linalg.norm(v)
    return v / n if n > 1e-12 else v


def _perp(u):
    return np.array([-u[1], u[0]])


def _rect(ax, center, u, length, width, cross=False):
    w = _perp(u)
    c = [center + s1 * u * length / 2 + s2 * w * width / 2 for s1, s2 in ((-1, -1), (1, -1), (1, 1), (-1, 1))]
    ax.add_patch(Polygon(c, closed=True, facecolor="white", edgecolor=INK, lw=1.4, zorder=4))
    if cross:
        ax.plot(*np.array([c[0], c[2]]).T, color=INK, lw=1, zorder=5)
        ax.plot(*np.array([c[1], c[3]]).T, color=INK, lw=1, zorder=5)


def _bearing(ax, center, u, s, half):
    """Revolute joint with its axis in the drawing plane: the shaft runs along u between two plates, ]|[."""
    w = _perp(u)
    gap, tick = 0.35 * s, 0.35 * s
    # white gap hides the links running into the joint; the shaft is drawn through it
    corners = [center + a * half * u + b * gap * w for a, b in ((-1, -1), (1, -1), (1, 1), (-1, 1))]
    ax.add_patch(Polygon(corners, closed=True, facecolor="white", edgecolor="none", zorder=4))
    ax.plot(*np.array([center - half * u, center + half * u]).T, color=INK, lw=2, zorder=5)
    for sgn in (-1, 1):
        a, b = center + sgn * gap * w - half * u, center + sgn * gap * w + half * u
        ax.plot(*np.array([a, b]).T, color=INK, lw=1.6, zorder=5)
        for end in (a, b):
            ax.plot(*np.array([end, end + sgn * tick * w]).T, color=INK, lw=1.6, zorder=5)


def _leader(ax, label_at, target, target_r, s):
    """Thin leader from just outside the label text to the target (dot when it ends on a line)."""
    d = target - label_at
    length = np.linalg.norm(d)
    u = d / max(length, 1e-12)
    start, end = label_at + 0.6 * s * u, target - target_r * u
    if length - 0.6 * s - target_r < 0.25 * s:
        return
    ax.plot(*np.array([start, end]).T, color=INK_MUTED, lw=0.6, zorder=6)
    if target_r == 0:
        ax.plot(*end, marker="o", markersize=2.2, color=INK_MUTED, zorder=6)


def _label(ax, at, text, size=13, color=INK):
    ax.text(*at, f"${text}$", fontsize=size, color=color, ha="center", va="center", zorder=6)


def _ground(ax, at, s, direction=np.array([0.0, -1.0])):
    w = _perp(direction)
    a, b = at - 1.4 * s * w, at + 1.4 * s * w
    ax.plot(*np.array([a, b]).T, color=INK, lw=1.8, zorder=3)
    for t in np.linspace(0, 1, 7):
        p = a + t * (b - a)
        ax.plot(*np.array([p, p + 0.5 * s * (direction - w)]).T, color=INK, lw=0.9, zorder=3)


def _gripper(ax, at, u, s):
    w = _perp(u)
    base = at + 0.3 * s * u
    ax.plot(*np.array([at, base]).T, color=INK, lw=2, zorder=3)
    ax.plot(*np.array([base - 0.6 * s * w, base + 0.6 * s * w]).T, color=INK, lw=2, zorder=3)
    for sgn in (-1, 1):
        p = base + sgn * 0.6 * s * w
        ax.plot(*np.array([p, p + 0.6 * s * u]).T, color=INK, lw=2, zorder=3)


# ---------------------------------------------------------------- drawing

def draw_schematic(frames, joints: list[JointInfo], path: Path, view="auto", base_axes=True):
    frames = [np.asarray(f, dtype=float) for f in frames]
    planar = all(abs(f[2, 3] - frames[0][2, 3]) < 1e-9 and abs(abs(f[2, 2]) - 1) < 1e-9 for f in frames)
    if view == "auto":
        view = "top" if planar else "oblique"
    P = VIEWS[view]

    o = [f[:3, 3] for f in frames]
    n = len(joints)
    # m_i: foot of the common normal on z_{i-1}; o_i = m_i + a_i x_i
    m = [o[i] - np.dot(o[i] - o[i - 1], frames[i][:3, 0]) * frames[i][:3, 0] for i in range(1, n + 1)]

    reach = sum(np.linalg.norm(m[i - 1] - o[i - 1]) + np.linalg.norm(o[i] - m[i - 1]) for i in range(1, n + 1))
    s = 0.06 * max(reach, 1e-3)

    fig, ax = plt.subplots(figsize=(6.5, 5))
    ax.set_aspect("equal")
    ax.axis("off")

    # labels take the first candidate position that keeps clear of what is already drawn;
    # obstacles are (point, keep-out radius): joints, base, gripper, link lines, earlier labels
    u_end = P(o[n]) - P(m[n - 1])
    if np.linalg.norm(u_end) < 1e-9:
        u_end = P(o[n] + frames[n][:3, 2]) - P(o[n])
    u_end = _unit(u_end)
    obstacles = [(P(p), 1.2 * s) for p in o]
    obstacles += [(P(o[0]) + np.array([t * s, -0.4 * s]), 1.0 * s) for t in (-1.2, 0, 1.2)]
    obstacles += [(P(o[n]) + 0.6 * s * u_end, 1.3 * s)]
    for i in range(1, n + 1):
        for a, b in ((o[i - 1], m[i - 1]), (m[i - 1], o[i])):
            for t in np.linspace(0, 1, 7):
                obstacles.append((P(a + t * (b - a)), 0.7 * s))

    def place(candidates, text, anchor=None, anchor_r=0.0):
        """candidates: positions, or (position, anchor) pairs when each has its own anchor point.
        A thin leader line joins the label to its anchor (a dot on a link, or a joint symbol's edge)."""
        candidates = [c if isinstance(c, tuple) else (c, anchor) for c in candidates]

        def clearance(c):
            return min(np.linalg.norm(c[0] - b) - r for b, r in obstacles)
        pos, target = next((c for c in candidates if clearance(c) >= 0), max(candidates, key=clearance))
        obstacles.append((pos, 1.4 * s))
        _label(ax, pos, text)
        if target is not None:
            _leader(ax, pos, target, anchor_r, s)
            for t in (0.35, 0.65):
                obstacles.append((pos + t * (target - pos), 0.4 * s))

    # base coordinate axes (their arrows and labels become label obstacles)
    axis_marks = []
    if base_axes:
        L = 0.55 * reach
        names = ("x_0", "y_0", "z_0") if view == "oblique" else ("x_0", "y_0")
        for k, name in enumerate(names):
            e = np.zeros(3)
            e[k] = L
            a, b = P(o[0]), P(o[0] + e)
            ax.annotate("", xy=b, xytext=a, zorder=1, annotation_clip=False,
                        arrowprops=dict(arrowstyle="-|>", color=INK_MUTED, lw=0.9, mutation_scale=12))
            ax.update_datalim([a, b + 1.8 * s * _unit(b - a)])
            at = b + 0.9 * s * _unit(b - a) + 0.6 * s * _perp(_unit(b - a))
            _label(ax, at, name.upper(), size=12, color=INK_MUTED)
            axis_marks.append((at, 1.4 * s))
            axis_marks += [(a + t * (b - a), 0.5 * s) for t in np.linspace(0.2, 1, 6)]

    obstacles += axis_marks
    _ground(ax, P(o[0]) - (0.0 if view == "oblique" else 1.2 * s) * np.array([0, 1]), s)

    # links: o_{i-1} -> m_i -> o_i; d_i / a_i labels are placed last, beside the segment midpoints
    link_labels = []
    for i in range(1, n + 1):
        j = joints[i - 1]
        # x_{i-1}: reference of the angle arc at joint i; a_i labels go on the other side
        ref = P(o[i - 1] + frames[i - 1][:3, 0]) - P(o[i - 1])
        for (a, b), text, is_d in (((o[i - 1], m[i - 1]), j.d_label, True), ((m[i - 1], o[i]), j.a_label, False)):
            if np.linalg.norm(b - a) < 1e-9:
                continue
            pa, pb = P(a), P(b)
            ax.plot(*np.array([pa, pb]).T, color=INK, lw=2, solid_capstyle="round", zorder=2)
            if text and not (j.kind == "prismatic" and text == j.label):
                w = _perp(_unit(pb - pa))
                if (is_d and w[0] < 0) or (not is_d and np.dot(w, ref) > 0):
                    w = -w
                # preferred side first, sliding along the segment, then the other side, then further out
                along = [pa + t * (pb - pa) for t in (0.5, 0.7, 0.3)]
                candidates = ([(q + 1.6 * s * w, q) for q in along] + [(q - 1.6 * s * w, q) for q in along]
                              + [(along[0] + 2.4 * s * w, along[0]), (along[0] - 2.4 * s * w, along[0])])
                link_labels.append((candidates, text))

    # joints
    for i in range(1, n + 1):
        j = joints[i - 1]
        f = frames[i - 1]
        center3 = o[i - 1]
        axis2 = P(center3 + f[:3, 2]) - P(center3)
        along_d = m[i - 1] - o[i - 1]
        toward_viewer = np.linalg.norm(axis2) < 0.6

        if toward_viewer:
            c = P(center3)
            corners = [c + 1.3 * s * np.array(d) for d in ((-1, 1), (1, 1), (-1, -1), (1, -1))]
            ua = None
            if view == "oblique" and np.linalg.norm(axis2) > 1e-6:
                # foreshortened axis (e.g. along x0): draw it as a dashed line in its projected direction
                ua = _unit(axis2)
                ax.plot(*np.array([c - 2.4 * s * ua, c + 2.4 * s * ua]).T,
                        color=INK_MUTED, lw=1, ls=(0, (4, 3)), zorder=1)
                for t in np.linspace(-2.4, 2.4, 5):
                    obstacles.append((c + t * s * ua, 0.6 * s))
                corners = [c + 2.0 * s * ua + 1.2 * s * _perp(ua), c + 2.0 * s * ua - 1.2 * s * _perp(ua),
                           c - 2.0 * s * ua + 1.2 * s * _perp(ua), c - 2.0 * s * ua - 1.2 * s * _perp(ua)] + corners
            if j.kind == "revolute":
                if ua is not None:
                    # a disc seen at an angle: long side across the axis, as on the exercise sheets
                    tilt = np.degrees(np.arctan2(ua[1], ua[0])) + 90
                    ax.add_patch(Ellipse(c, 0.9 * s, 0.5 * s, angle=tilt,
                                         facecolor="white", edgecolor=INK, lw=1.4, zorder=4))
                else:
                    ax.add_patch(Circle(c, 0.55 * s, facecolor="white", edgecolor=INK, lw=1.4, zorder=4))
                arc = False
                if view == "top":
                    link = np.linalg.norm(P(o[i]) - P(m[i - 1]))
                    r = min(2.4 * s, max(0.4 * link, 1.1 * s))
                    arc = _angle_arc(ax, P, f, frames[i], c, s, j.theta_label, r, place)
                if not arc:
                    place(corners, j.label, anchor=c, anchor_r=0.5 * s)
            else:
                _rect(ax, c, np.array([1.0, 0.0]), 1.2 * s, 1.2 * s, cross=True)
                place(corners, j.label, anchor=c, anchor_r=0.7 * s)
        else:
            u = _unit(axis2)
            d_len = np.linalg.norm(along_d)
            d_sign = np.sign(np.dot(along_d, f[:3, 2])) or 1.0
            if j.kind == "revolute":
                # centred on the joint so the incoming link meets it in the middle (--]|[);
                # the first joint instead stands on the base
                half = min(s, 0.3 * d_len) if d_len > 1e-9 else s
                shift = d_sign * half if (i == 1 and d_len > 1e-9) else 0.0
                c = P(center3 + shift * f[:3, 2])
                _bearing(ax, c, u, s, half)
                u_d = u * d_sign
                side = _perp(u) * (-1 if _perp(u)[0] > 0 else 1)
                place([c - (half + 1.2 * s) * u_d, c + 1.8 * s * side, c - 1.8 * s * side], j.label,
                      anchor=c, anchor_r=0.6 * s)
            else:
                # slide the box along the d offset so it sits on the link, not on the previous joint
                shift = d_sign * min(1.6 * s, d_len / 2)
                c = P(center3 + shift * f[:3, 2])
                _rect(ax, c, u, 2.2 * s, 1.1 * s, cross=True)
                side = _perp(u) * (-1 if _perp(u)[0] > 0 else 1)
                place([c + 2.0 * s * side, c - 2.0 * s * side], j.label, anchor=c, anchor_r=0.7 * s)

    # end effector, pointing along the last link
    _gripper(ax, P(o[n]), u_end, s)

    for candidates, text in link_labels:
        place(candidates, text)

    ax.autoscale_view()
    ax.margins(0.08)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def _angle_arc(ax, P, f_prev, f_next, c, s, label, r, place):
    """Dashed x_{i-1} reference and an arc to x_i, for joints whose axis points at the viewer."""
    ref = _unit(P(f_prev[:3, 3] + f_prev[:3, 0]) - P(f_prev[:3, 3]))
    tgt = _unit(P(f_next[:3, 3] + f_next[:3, 0]) - P(f_next[:3, 3]))
    a0, a1 = np.degrees(np.arctan2(ref[1], ref[0])), np.degrees(np.arctan2(tgt[1], tgt[0]))
    sweep = (a1 - a0 + 180) % 360 - 180
    if abs(sweep) < 4:
        return False
    ax.plot(*np.array([c, c + 1.35 * r * ref]).T, color=INK_MUTED, lw=0.9, ls=(0, (4, 3)), zorder=1)
    start, end = (a0, a0 + sweep) if sweep > 0 else (a0 + sweep, a0)
    ax.add_patch(Arc(c, 2 * r, 2 * r, theta1=start, theta2=end, color=INK, lw=0.9, zorder=3))
    mid = np.radians(a0 + sweep / 2)
    # stay with the arc: along its bisector, or just past the dashed reference line
    outside = np.radians(a0 - np.sign(sweep) * 15)
    e_mid = np.array([np.cos(mid), np.sin(mid)])
    on_ref = c + r * ref
    candidates = [(c + (r + 1.2 * s) * e_mid, c + r * e_mid),
                  (c + (r + 1.9 * s) * e_mid, c + r * e_mid),
                  (c + (r + 1.2 * s) * np.array([np.cos(outside), np.sin(outside)]), on_ref)]
    place(candidates, label)
    return True
