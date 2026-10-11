// Kinematics walkthrough document. Content comes from gen_walkthrough.py:
//   typst compile --input data=build/<name>/walkthrough.json walkthrough.typ out/<name>.pdf
#import "@preview/mitex:0.2.7": mitex, mi

#let data = json(sys.inputs.at("data", default: "build/planar2/walkthrough.json"))

#set document(title: data.name + " — kinematics walkthrough")
#set page(paper: "a4", margin: (x: 2cm, y: 2.2cm), numbering: "1 / 1")
#set text(size: 10.5pt, lang: "en")
#set par(justify: true)
#set heading(numbering: "1.1")
#show heading.where(level: 1): set block(above: 1.6em, below: 0.9em)
#show figure.caption: set text(size: 9pt)

#let note(body) = block(
  fill: luma(245), inset: 8pt, radius: 3pt, width: 100%, text(size: 9.5pt, body),
)
#let values-table(rows, header: ([Variable], [Value])) = table(
  columns: 2, align: (center, right), stroke: none,
  table.hline(), table.header(..header.map(strong)), table.hline(stroke: 0.5pt),
  ..rows.map(r => (mi(r.name), r.value)).flatten(),
  table.hline(),
)
#let check(ok) = if ok { text(fill: rgb("#008300"))[✔ passed] } else { text(fill: rgb("#e34948"))[✘ FAILED] }

// ---------------------------------------------------------------- title
#align(center)[
  #text(size: 18pt, weight: "bold", data.name) \
  #v(0.2em)
  #text(size: 11pt, fill: luma(90))[Forward and inverse kinematics walkthrough (Denavit–Hartenberg)]
]
#v(0.8em)
#data.description

// ---------------------------------------------------------------- manipulator
= Manipulator and DH table

#if data.figure != none [
  #figure(image(data.figure, width: 55%),
    caption: [The manipulator as given; compare with the generated drawings below.])
]

#figure(
  grid(columns: 2, gutter: 0.8em, align: horizon,
    image(data.figures.schematic, width: 100%), image(data.figures.model3d, width: 100%)),
  caption: [Left: kinematic schematic drawn from the DH table. Right: 3-D view of the same pose
    (revolute joints blue, prismatic orange; dashed lines drop to the $X_0 Y_0$ plane, the grey trace
    is the arm's projection onto it).],
)

The manipulator has #data.joints.len() joints. Joint variables:
#list(..data.joints.map(j => [#mi(j.name) — #j.type]))

#figure(
  table(
    columns: 5, align: center, stroke: none, inset: 6pt,
    table.hline(),
    table.header([$i$], [$theta_i$], [$d_i$], [$a_i$], [$alpha_i$]),
    table.hline(stroke: 0.5pt),
    ..data.steps.map(s => (str(s.i), ..s.row.map(mi))).flatten(),
    table.hline(),
  ),
  caption: [Denavit–Hartenberg parameters.],
)

#if data.constants.len() > 0 [
  Numeric link parameters (used only for the figures and the numeric IK):
  #data.constants.map(c => [#mi(c.name) $= #c.value$]).join([, ]).
]

// ---------------------------------------------------------------- DH formula
= Denavit–Hartenberg transformation

Each row of the table defines the homogeneous transformation from frame $i-1$ to frame $i$:
$ A_i = "Rot"_z (theta_i) dot "Trans"_z (d_i) dot "Trans"_x (a_i) dot "Rot"_x (alpha_i), $
#mitex("A_i = " + data.A_general)
and the pose of frame $i$ in the base frame is accumulated as
$ T_0^i = T_0^(i-1) A_i, quad T_0^0 = I. $

#if data.legend.len() > 0 [
  #note[
    *Notation.* To keep the matrices readable:
    #data.legend.map(l => mi(l.short + " = " + l.full)).join([, ]).
  ]
]

// ---------------------------------------------------------------- FK steps
= Forward kinematics, step by step

#for s in data.steps [
  == Joint #s.i
  DH row: #mi("\\theta_{" + str(s.i) + "} = " + s.row.at(0)),
  #mi("d_{" + str(s.i) + "} = " + s.row.at(1)),
  #mi("a_{" + str(s.i) + "} = " + s.row.at(2)),
  #mi("\\alpha_{" + str(s.i) + "} = " + s.row.at(3)).
  #mitex("A_{" + str(s.i) + "} = " + s.A)
  #if s.i == 1 [
    Since $T_0^0 = I$, the first cumulative transformation is $A_1$ itself.
  ] else [
    #mitex("T_0^{" + str(s.i) + "} = T_0^{" + str(s.i - 1) + "} A_{" + str(s.i) + "} = " + s.T)
  ]
]

// ---------------------------------------------------------------- position
= End-effector position

The position of the end effector is the translation part of $T_0^#data.steps.len()$:
#mitex("p(q) = " + data.position.coords + " = " + data.position.expr)

#if data.task_space.fixed.len() > 0 [
  The coordinate(s) #data.task_space.fixed.map(f => mi(f.name + " = " + f.value)).join([, ])
  do not depend on the joint variables, so the reachable positions are described by
  #data.task_space.coords.map(mi).join([, ]) only
  (task space of dimension $m = #data.task_space.m$).
]

#figure(
  grid(
    columns: (1fr, auto), gutter: 1.2em, align: horizon,
    image(data.figures.frames, width: 100%),
    stack(
      spacing: 1em,
      values-table(data.sample.q, header: ([Joint], [Sample])),
      values-table(data.sample.target, header: ([Coord.], [Value])),
    ),
  ),
  caption: [Frames $0 dots #data.steps.len()$ at the sample configuration
    ($x_i$ red, $y_i$ green, $z_i$ blue; in planar views $z_i$ points out of ($circle.filled.small$) or into ($times$) the page).
    Compare with the sketch to check the frame assignment.],
)

// ---------------------------------------------------------------- jacobian
= Position Jacobian

Differentiating the task-space position with respect to the joint variables gives the
Jacobian used by the numeric inverse kinematics:
#mitex("J(q) = \\frac{\\partial p}{\\partial q}, \\quad p = " + data.jacobian.p + ", \\quad q = " + data.jacobian.q)
#mitex("J(q) = " + data.jacobian.J)

// ---------------------------------------------------------------- singularities
= Singular configurations

#let sg = data.singular
A configuration is singular when $J(q)$ loses rank. The end effector then cannot move in some
direction, however fast the joints turn, and some joint motion leaves it in place.

#if sg.idle.len() > 0 [
  #note[
    *Joints that never move the end effector.* The column of $J$ for
    #sg.idle.map(mi).join([, ]) is identically zero, so these joints only change the orientation of
    the end effector, never its position. The analysis below uses #sg.active.map(mi).join([, ]).
  ]
]

#let tuple(names) = mi("(" + names.join(", ") + ")")
#let pose-table(rows) = table(
  columns: (auto, 1fr, 1fr), align: (left, center, center), stroke: none, inset: 5pt,
  table.hline(),
  table.header([*Configuration*], [*Cannot move along* #tuple(sg.coords)],
    [*Joint motion that leaves it in place* #tuple(sg.joint_names)]),
  table.hline(stroke: 0.5pt),
  ..rows.map(r => (
    r.q.map(v => [#mi(v.name) $=$ #v.value]).join([, ]),
    if r.lost == none { [not unique] } else { raw(r.lost) },
    if r.still == none { [not unique] } else { raw(r.still) },
  )).flatten(),
  table.hline(),
)
#let uniqueness-note = if sg.kind == "minors" and sg.by == "rows" [
  With fewer moving joints than coordinates, the effector already cannot reach every direction, and at
  a singularity it loses one more, so the lost direction is not unique. The joint motion is.
] else if sg.kind == "minors" [
  A redundant arm can always move its joints without moving the effector (self-motion), so that
  motion is not unique. The lost direction is.
]

#let families(fs) = for (i, f) in fs.enumerate() [
  == Family #(i + 1): #mi(f.cond + " = 0")
  #if f.solution != none [
    Solving for the joint: #mi(f.solution).
  ] else [
    This condition couples #f.joints.map(mi).join([, ]). An example is found numerically by solving
    it for one of them while the others keep their sample values.
  ]
  #if f.examples.len() > 0 { pose-table(f.examples) }
]

#uniqueness-note

#if sg.kind == "square" [
  $J$ is square, so it is singular exactly where its determinant vanishes:
  #mitex("\\det J(q) = " + sg.det)
  #if sg.factors.len() == 0 [
    The determinant has no factor that depends on the joints, so the manipulator has no singular
    configurations.
  ] else [
    Every factor that depends on the joints gives one family of singular configurations. In the
    examples, the other joints keep their sample values.
    #families(sg.factors)
  ]
] else [
  $J$ has #sg.m rows but #sg.na moving joints, so it is not square. It loses rank exactly when every
  #(sg.k)×#(sg.k) minor vanishes. By the Cauchy–Binet formula, the Gram determinant
  (#if sg.by == "columns" [$det(J J^T)$] else [$det(J^T J)$]) is the sum of their squares. The minors are:
  #for mn in sg.minors { mitex(mn) }
  #if sg.factors.len() > 0 [
    All the minors share a factor, so each of its zeros makes them vanish at once:
    #families(sg.factors)
  ]
  #if sg.example != none [
    *Numerical example.* Gradient descent on the product of the #sg.k largest singular values,
    started from the sample configuration, ends at
    #pose-table((sg.example,))
  ]
]

// ---------------------------------------------------------------- IK
= Inverse kinematics

Given a target position, find joint variables $q$ with $p(q) = p_"target"$.
The target used below is the position reached at the sample configuration, so it is
guaranteed to be reachable.

#if data.task_space.n > data.task_space.m [
  #note[
    *Redundancy.* The manipulator has $n = #data.task_space.n$ joints but the task space has only
    $m = #data.task_space.m$ coordinates, so a reachable position has infinitely many solutions
    ($n - m = #(data.task_space.n - data.task_space.m)$-dimensional self-motion). A unique closed form needs extra
    constraints; the numeric method below returns one solution close to the initial guess.
  ]
]

== Closed-form solution

#if data.ik.closed == none [
  No closed-form solution is given in the manipulator spec.
] else [
  #let closed = data.ik.closed
  #if closed.intermediates.len() > 0 [
    Auxiliary quantities:
    #for it in closed.intermediates { mitex(it.name + " = " + it.expr) }
  ]
  #for sol in closed.solutions [
    === #sol.name
    #for a in sol.assign { mitex(a.name + " = " + a.expr) }
    At the target this gives #sol.at_target.map(v => [#mi(v.name) $= $ #v.value]).join([, ]),
    with position error #sol.error.

    *Verification:* substituted back into forward kinematics for #sol.verify.trials random
    reachable targets, the largest position error is #sol.verify.max_error
    #if sol.verify.failed > 0 [ (#sol.verify.failed targets gave no real solution)]
    — #check(sol.verify.ok).
  ]
]

== Numeric solution (damped least squares)

Starting from $q_0$, iterate
$ e_k = p_"target" - p(q_k), quad
  q_(k+1) = q_k + J^T (J J^T + lambda^2 I)^(-1) e_k, quad
  lambda = #data.ik.numeric.damping, $
until $norm(e_k) < 10^(-10)$. The damping term keeps the step bounded near singular configurations.

#let num = data.ik.numeric
#grid(
  columns: (1fr, 1fr, 1fr), gutter: 1em,
  values-table(num.q0, header: ([Joint], [$q_0$])),
  table(
    columns: 2, align: (center, right), stroke: none,
    table.hline(), table.header([*$k$*], [*$norm(e_k)$*]), table.hline(stroke: 0.5pt),
    ..num.iterations.map(it => if it == none { ([⋮], []) } else { (str(it.k), it.err) }).flatten(),
    table.hline(),
  ),
  values-table(num.result, header: ([Joint], [Result])),
)
#if num.converged [Converged.] else [#text(fill: rgb("#e34948"))[*Did not converge* — try another initial guess.]]

#figure(
  image(data.figures.ik, width: 75%),
  caption: [Inverse-kinematics solutions reaching the target ($star$); the numeric solution is dashed, the sample configuration that defined the target is shown in grey.],
)
