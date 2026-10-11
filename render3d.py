"""Shaded 3-D view of a manipulator built from its DH frames.

Links are cylinders along o_{i-1} -> m_i -> o_i, revolute joints are wider
cylinders along z_{i-1}, prismatic joints are sleeves, and the gripper is
oriented by the last frame. Faces are shaded by hand (Lambert), back faces are
culled and everything goes into one Poly3DCollection so matplotlib depth-sorts
the faces together.
"""

from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import to_rgb  # noqa: E402
from mpl_toolkits.mplot3d.art3d import Poly3DCollection  # noqa: E402

INK = "#0b0b0b"
INK_MUTED = "#52514e"
FLOOR = "#d9d8d4"
SHADOW = "#c3c2bb"
COLORS = {"link": "#a7adb5", "revolute": "#2a78d6", "prismatic": "#eb6834", "base": "#52514e", "gripper": "#52514e"}

plt.rcParams["mathtext.fontset"] = "cm"

ELEV, AZIM = 22, 25   # viewer on the +x side, slightly toward +y: y0 to the right, z0 up
LIGHT = np.array([0.5, -0.3, 1.0]) / np.linalg.norm([0.5, -0.3, 1.0])


def _unit(v):
    n = np.linalg.norm(v)
    return v / n if n > 1e-12 else v


def _basis(d):
    """Two unit vectors perpendicular to d."""
    a = np.array([1.0, 0, 0]) if abs(d[0]) < 0.9 else np.array([0, 1.0, 0])
    u = _unit(np.cross(d, a))
    return u, np.cross(d, u)


def _cylinder(p0, p1, r, k=36):
    """Faces (vertices, outward normal) of a closed cylinder, split along its length for depth sorting."""
    d = p1 - p0
    length = np.linalg.norm(d)
    if length < 1e-9:
        return []
    d = d / length
    u, v = _basis(d)
    ang = np.linspace(0, 2 * np.pi, k, endpoint=False)
    ring = [np.cos(a) * u + np.sin(a) * v for a in ang]
    steps = max(1, int(np.ceil(length / (2.5 * r))))
    faces = []
    for t0, t1 in zip(np.linspace(0, length, steps + 1)[:-1], np.linspace(0, length, steps + 1)[1:]):
        for a in range(k):
            b = (a + 1) % k
            quad = [p0 + t0 * d + r * ring[a], p0 + t0 * d + r * ring[b],
                    p0 + t1 * d + r * ring[b], p0 + t1 * d + r * ring[a]]
            faces.append((quad, _unit(ring[a] + ring[b])))
    faces.append(([p0 + r * q for q in ring[::-1]], -d))
    faces.append(([p1 + r * q for q in ring], d))
    return faces


def _box(center, axes, half):
    """Faces of a box with orthonormal axes (3 vectors) and half-sizes (3 numbers)."""
    faces = []
    for k in range(3):
        i, j = [a for a in range(3) if a != k]
        for sgn in (-1, 1):
            c = center + sgn * half[k] * axes[k]
            quad = [c + si * half[i] * axes[i] + sj * half[j] * axes[j] for si, sj in ((-1, -1), (1, -1), (1, 1), (-1, 1))]
            faces.append((quad, sgn * axes[k]))
    return faces


def _shade(color, normal):
    base = np.array(to_rgb(color))
    return tuple(np.clip(base * (0.45 + 0.6 * max(0.0, float(np.dot(normal, LIGHT)))), 0, 1))


def render_3d(frames, kinds, labels, path: Path):
    frames = [np.asarray(f, dtype=float) for f in frames]
    n = len(kinds)
    o = [f[:3, 3] for f in frames]
    z = [f[:3, 2] for f in frames]
    m = [o[i] - np.dot(o[i] - o[i - 1], frames[i][:3, 0]) * frames[i][:3, 0] for i in range(1, n + 1)]
    chain = [o[0]]
    for i in range(1, n + 1):
        chain += [m[i - 1], o[i]]
    reach = sum(np.linalg.norm(b - a) for a, b in zip(chain, chain[1:]))
    r = 0.022 * max(reach, 1e-3)
    z0 = o[0][2]

    faces = []   # (vertices, normal, color)

    def add(fs, color):
        faces.extend((v, nrm, color) for v, nrm in fs)

    # base plate
    add(_box(o[0] - np.array([0, 0, 0.6 * r]), np.eye(3), (5 * r, 5 * r, 0.6 * r)), COLORS["base"])

    # links
    for a, b in zip(chain, chain[1:]):
        add(_cylinder(a, b, r), COLORS["link"])

    # joints
    label_at = []
    for i in range(1, n + 1):
        axis = _unit(z[i - 1])
        d_vec = m[i - 1] - o[i - 1]
        if kinds[i - 1] == "revolute":
            c = o[i - 1] + (3 * r * axis if i == 1 else 0)
            if i == 1 and np.dot(d_vec, axis) < 0:
                c = o[i - 1] - 3 * r * axis
            add(_cylinder(c - 2.6 * r * axis, c + 2.6 * r * axis, 2.0 * r), COLORS["revolute"])
        else:
            sgn = np.sign(np.dot(d_vec, axis)) or 1.0
            c = o[i - 1] + sgn * 4 * r * axis
            u, v = _basis(axis)
            add(_box(c, [u, v, axis], (2.0 * r, 2.0 * r, 4 * r)), COLORS["prismatic"])
        label_at.append(c)

    # gripper: stem along the approach direction, jaws across it (in the last frame, so it turns with q_n)
    approach = _unit(chain[-1] - next(p for p in reversed(chain[:-1]) if np.linalg.norm(chain[-1] - p) > 1e-9))
    R = frames[n][:3, :3]
    jaw = next(R[:, k] for k in (1, 0, 2) if abs(np.dot(R[:, k], approach)) < 0.5)
    jaw = _unit(jaw - np.dot(jaw, approach) * approach)
    side = np.cross(approach, jaw)
    end = o[n]
    bar = end + 1.5 * r * approach
    add(_cylinder(end, bar, 0.8 * r), COLORS["gripper"])
    add(_box(bar, [jaw, side, approach], (4 * r, 0.8 * r, 0.8 * r)), COLORS["gripper"])
    for s in (-1, 1):
        add(_box(bar + s * 3.4 * r * jaw + 3 * r * approach, [jaw, side, approach], (0.6 * r, 0.8 * r, 3 * r)),
            COLORS["gripper"])

    # camera: cull faces turned away, shade the rest
    view = np.array([np.cos(np.radians(ELEV)) * np.cos(np.radians(AZIM)),
                     np.cos(np.radians(ELEV)) * np.sin(np.radians(AZIM)),
                     np.sin(np.radians(ELEV))])
    visible = [(v, _shade(c, nrm)) for v, nrm, c in faces if np.dot(nrm, view) > -0.05]

    fig = plt.figure(figsize=(6.5, 5))
    ax = fig.add_subplot(111, projection="3d")
    ax.computed_zorder = False
    ax.set_proj_type("persp", focal_length=0.5)
    ax.view_init(elev=ELEV, azim=AZIM)
    ax.set_axis_off()

    pts = np.array(chain)
    lo, hi = pts.min(axis=0), pts.max(axis=0)
    span = max((hi - lo).max(), 1e-3)
    pad = 0.25 * span

    # floor grid and the arm's shadow trace: the depth cues
    gx = np.linspace(min(lo[0], 0) - pad, max(hi[0], 0) + pad, 7)
    gy = np.linspace(min(lo[1], 0) - pad, max(hi[1], 0) + pad, 7)
    for x in gx:
        ax.plot([x, x], [gy[0], gy[-1]], [z0, z0], color=FLOOR, lw=0.6)
    for y in gy:
        ax.plot([gx[0], gx[-1]], [y, y], [z0, z0], color=FLOOR, lw=0.6)
    ax.plot(pts[:, 0], pts[:, 1], np.full(len(pts), z0), color=SHADOW, lw=4, solid_capstyle="round")
    for p in (m[-1], o[-1]):
        ax.plot([p[0], p[0]], [p[1], p[1]], [z0, p[2]], color=INK_MUTED, lw=0.6, ls=(0, (4, 3)))

    # base axes
    L = 0.45 * reach
    for k, name in enumerate(("X_0", "Y_0", "Z_0")):
        e = np.zeros(3)
        e[k] = L
        ax.plot(*np.array([o[0], o[0] + e]).T, color=INK_MUTED, lw=0.9)
        ax.text(*(o[0] + 1.08 * e), f"${name}$", color=INK_MUTED, fontsize=11, ha="center", va="center")

    ax.add_collection3d(Poly3DCollection([v for v, _ in visible], facecolors=[c for _, c in visible],
                                         edgecolors=[c for _, c in visible], linewidths=0.7, zsort="average",
                                         rasterized=True))  # avoids hairline seams between faces in the SVG

    # labels to the screen-right of each joint, clear of the joint cylinder
    right = _unit(np.cross([0, 0, 1], view))
    for c, text in zip(label_at, labels):
        ax.text(*(c + 5.5 * r * right + 1.5 * r * np.array([0, 0, 1])), f"${text}$",
                color=INK, fontsize=13, ha="center", va="center", zorder=10)

    # equal scale on all axes, box fitted to the content
    lo3 = np.minimum(lo, [gx[0], gy[0], z0])
    hi3 = np.maximum(hi, [gx[-1], gy[-1], z0]) + np.array([0, 0, 2 * r])
    ax.set_xlim(lo3[0], hi3[0])
    ax.set_ylim(lo3[1], hi3[1])
    ax.set_zlim(lo3[2], hi3[2])
    ax.set_box_aspect(hi3 - lo3)
    fig.subplots_adjust(0, 0, 1, 1)
    fig.savefig(path, dpi=220)
    plt.close(fig)
