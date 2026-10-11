// Group discussion document (definition of a group, SO(3), SE(3)). Content comes from gen_groups.py:
//   typst compile --input data=build/groups/groups.json groups.typ out/groups.pdf
#import "@preview/mitex:0.2.7": mitex, mi

#let data = json(sys.inputs.at("data", default: "build/groups/groups.json"))
#let so3 = data.so3
#let se3 = data.se3
#let man = data.manipulator

#set document(title: "Groups, SO(3) and SE(3)")
#set page(paper: "a4", margin: (x: 2cm, y: 2.2cm), numbering: "1 / 1")
#set text(size: 10.5pt, lang: "en")
#set par(justify: true)
#set heading(numbering: "1.1")
#show heading.where(level: 1): set block(above: 1.6em, below: 0.9em)
#show figure.caption: set text(size: 9pt)

#let note(body) = block(
  fill: luma(245), inset: 8pt, radius: 3pt, width: 100%, text(size: 9.5pt, body),
)
// every statement is checked by gen_groups.py; only a failed check is printed
#let check(ok) = if not ok { text(fill: rgb("#e34948"))[ — ✘ CHECK FAILED] }
#let SO = $"SO"(3)$
#let SE = $"SE"(3)$

// ---------------------------------------------------------------- title
#align(center)[
  #text(size: 18pt, weight: "bold")[Groups, SO(3) and SE(3)] \
  #v(0.2em)
  #text(size: 11pt, fill: luma(90))[Objects, operation, identity element, inverse element]
]
#v(0.8em)

// ---------------------------------------------------------------- definition
= How a group is defined

A *group* is a pair $(G, compose)$: a set $G$ of objects and one binary operation $compose$ on them,
such that the following four conditions (the group axioms) hold.

#table(
  columns: (auto, 1fr), stroke: none, inset: 6pt, align: (left, left),
  table.hline(),
  [*Closure*], [$a compose b in G$ for all $a, b in G$. The operation never leaves the set.],
  [*Associativity*], [$(a compose b) compose c = a compose (b compose c)$ for all $a, b, c in G$.],
  [*Identity element*], [There is $e in G$ with $e compose a = a compose e = a$ for every $a in G$.],
  [*Inverse element*], [Every $a in G$ has $a^(-1) in G$ with $a compose a^(-1) = a^(-1) compose a = e$.],
  table.hline(),
)

Commutativity ($a compose b = b compose a$) is *not* required. A group that has it is called
abelian. Two consequences of the axioms: the identity is unique, and so is the inverse of each
element.

#note[
  *Quick examples.* $(ZZ, +)$ is a group: identity $0$, inverse $-a$. $(RR without {0}, dot)$ is a
  group: identity $1$, inverse $1 slash a$. $(ZZ, dot)$ is not: $2$ has no integer inverse.
  Invertible $n times n$ matrices with matrix multiplication form the group $"GL"(n)$, and both
  groups below are subsets of such a group that are groups themselves (subgroups).
]

// ---------------------------------------------------------------- SO(3)
= The rotation group SO(3)

== Objects and operation

The objects are the $3 times 3$ real matrices that are orthogonal and have determinant $+1$, and
the operation is matrix multiplication:
$ SO = { R in RR^(3 times 3) : R^T R = I, quad det R = 1 }. $
"S" stands for special ($det R = 1$), "O" for orthogonal ($R^T R = I$).

Such a matrix describes the orientation of one coordinate frame with respect to another: its
columns are the $x$, $y$ and $z$ unit vectors of the rotated frame, written in the reference frame.
$R^T R = I$ says that these three vectors have unit length and are mutually perpendicular.
$det R = 1$ says that they form a right-handed frame. A matrix has 9 entries and $R^T R = I$ gives
6 independent equations, so a rotation has 3 degrees of freedom, for example the three angles in
#mitex("R(\\alpha, \\beta, \\gamma) = \\mathrm{Rot}_z(\\alpha)\\,\\mathrm{Rot}_y(\\beta)\\,\\mathrm{Rot}_x(\\gamma) = " + so3.R)
with $c_alpha = cos alpha$, $s_alpha = sin alpha$ and so on. Substituting this matrix into the
definition, sympy confirms $R^T R = I$ and $det R = 1$ for all angles#check(so3.checks.member).

== The four axioms

/ Closure: For $R_1, R_2 in SO$, the product is orthogonal,
  $(R_1 R_2)^T (R_1 R_2) = R_2^T R_1^T R_1 R_2 = R_2^T R_2 = I$, and
  $det(R_1 R_2) = det R_1 dot det R_2 = 1$. So $R_1 R_2 in SO$: two rotations one after the other
  are again a rotation. Checked symbolically with two general rotations (six angles)#check(so3.checks.closure).
/ Associativity: Matrix multiplication is associative, so $(R_1 R_2) R_3 = R_1 (R_2 R_3)$.
  Checked symbolically with three general rotations#check(so3.checks.associative).
/ Identity element: The identity matrix $I$, the rotation by zero angle. $I^T I = I$ and
  $det I = 1$, so $I in SO$, and $I R = R I = R$#check(so3.checks.identity).
/ Inverse element: $R^(-1) = R^T$, directly from $R^T R = I$. The transpose is again in #SO, because
  $(R^T)^T R^T = R R^T = I$ and $det R^T = det R = 1$. Inverting a rotation therefore costs
  nothing: transpose it#check(so3.checks.inverse).

== The group is not commutative

The order of rotations matters. Take the rotations by $90°$ about $x$ and about $z$:
#mitex("R_x = " + so3.Rx + ", \\qquad R_z = " + so3.Rz)
#mitex("R_x R_z = " + so3.RxRz + " \\;\\neq\\; R_z R_x = " + so3.RzRx)

#figure(
  image(data.figures.order, width: 88%),
  caption: [The same two $90°$ rotations composed in both orders ($x$ red, $y$ green, $z$ blue;
    the dotted lines are the reference frame). The resulting frames differ.],
)

The rotations about one fixed axis do commute:
$"Rot"_z (alpha) "Rot"_z (beta) = "Rot"_z (alpha + beta)$. They form the abelian subgroup $"SO"(2)$.

#note[
  *What breaks if the definition is changed.*
  - _Addition instead of multiplication._ The set is not closed:
    #mi("R_x + R_z = " + so3.sum) has determinant #mi(so3.sum_det) and is not orthogonal. The zero
    matrix, the identity of addition, is not a rotation either.
  - _Without $det R = 1$._ The orthogonal matrices alone still form a group, $"O"(3)$, which also
    contains reflections ($det R = -1$). The reflections alone do not form a group: the product of
    two reflections has determinant $+1$, and $I$ is not among them.
]

// ---------------------------------------------------------------- SE(3)
= The group of rigid motions SE(3)

== Objects and operation

The objects are the $4 times 4$ homogeneous transformation matrices built from a rotation and a
translation, and the operation is again matrix multiplication:
$ SE = { T = mat(delim: "[", R, p; 0_(1 times 3), 1) : R in SO, quad p in RR^3 }. $
"S" and "E" stand for special Euclidean. $T$ describes the pose of a rigid body or of a coordinate
frame: $R$ is its orientation and $p$ the position of its origin. It has $3 + 3 = 6$ degrees of
freedom. A point with coordinates $r$ in the moved frame has coordinates $R r + p$ in the
reference frame, which in homogeneous coordinates is the single product
$mat(delim: "[", R r + p; 1) = T mat(delim: "[", r; 1)$.
A general $T$ (three angles, three displacements) satisfies the definition; this is checked symbolically#check(se3.checks.member).

== The four axioms

/ Closure: Multiplying the blocks gives
  $ T_1 T_2 = mat(delim: "[", R_1, p_1; 0, 1) mat(delim: "[", R_2, p_2; 0, 1)
    = mat(delim: "[", R_1 R_2, R_1 p_2 + p_1; 0, 1). $
  $R_1 R_2 in SO$ by the closure of #SO, $R_1 p_2 + p_1 in RR^3$, and the last row is still
  $mat(delim: "[", 0, 0, 0, 1)$, so $T_1 T_2 in SE$. Both the block formula and the
  membership of the product are checked symbolically#check(se3.checks.product_rule and se3.checks.closure).
/ Associativity: Inherited from matrix multiplication. Checked symbolically with three general
  transformations#check(se3.checks.associative).
/ Identity element: The $4 times 4$ identity matrix $I_4$, that is $R = I$ and $p = 0$: no
  rotation and no translation#check(se3.checks.identity).
/ Inverse element: Solving $T T^(-1) = I_4$ with the block formula gives
  $ T^(-1) = mat(delim: "[", R^T, -R^T p; 0, 1), $
  which is in #SE because $R^T in SO$. Note that the translation part is $-R^T p$, not $-p$#check(se3.checks.inverse).

== Structure

#SE is not commutative either. For example
#mitex("T_a = " + se3.Ta + ", \\qquad T_b = " + se3.Tb)
#mitex("T_a T_b = " + se3.TaTb + " \\;\\neq\\; T_b T_a = " + se3.TbTa)
and the inverse of the first one, from the formula above, is
#mitex("T_a^{-1} = " + se3.Ta_inv)
with $T_a T_a^(-1) = I_4$#check(se3.Ta_inv_ok).

Two subgroups sit inside #SE: the pure rotations ($p = 0$), a copy of #SO, and the pure
translations ($R = I$), a copy of $(RR^3, +)$, since multiplying two translation matrices adds
their vectors. But #SE is not simply "rotation and translation side by side": in the product rule
the rotation $R_1$ acts on the translation $p_2$. This is the semidirect product
$SE = SO times.l RR^3$.

// ---------------------------------------------------------------- kinematics
= Why this matters for forward kinematics

Every Denavit–Hartenberg matrix is a product of four elementary rigid motions,
$ A_i = "Rot"_z (theta_i) dot "Trans"_z (d_i) dot "Trans"_x (a_i) dot "Rot"_x (alpha_i), $
each of them an element of #SE. The group axioms then say:

- *Closure:* $A_i in SE$, and so is $T_0^n = A_1 A_2 dots.c A_n$. The result of forward kinematics
  is always a valid pose: a rotation matrix and a position vector.
- *Associativity:* the product can be bracketed in any way, for example
  $T_0^n = T_0^k T_k^n$, which is what the step-by-step accumulation $T_0^i = T_0^(i-1) A_i$ uses.
- *Identity:* $T_0^0 = I_4$ is the starting value of that accumulation.
- *Inverse:* the base seen from the end effector is $T_n^0 = (T_0^n)^(-1)$, obtained with a
  transpose and one matrix–vector product instead of a general matrix inversion.
- *No commutativity:* the order of the DH rows, and of the four factors inside $A_i$, cannot be
  changed freely. (Factors along the same axis do commute:
  $"Rot"_z (theta) "Trans"_z (d) = "Trans"_z (d) "Rot"_z (theta)$.)

== Example: #man.name

For the #man.n DH rows of this manipulator, each $A_i$ is in #SE#check(man.factors_ok) and
their product is
#mitex("T_0^{" + str(man.n) + "} = " + man.T)
Its rotation block satisfies $R^T R = I$ and $det R = 1$ for all joint values#check(man.member). The inverse, from $T^(-1) = mat(delim: "[", R^T, -R^T p; 0, 1)$, is
#mitex("T_{" + str(man.n) + "}^0 = " + man.T_inv)
and multiplying the two gives $I_4$#check(man.inverse_ok).

#if man.legend.len() > 0 [
  #note[
    *Notation.* #man.legend.map(l => mi(l.short + " = " + l.full)).join([, ]).
  ]
]

== The same example with numbers

#let ex = man.numeric
Take the configuration #ex.q.map(v => [#mi(v.name) $=$ #v.value]).join([, ])#if ex.constants.len() > 0 [
  and the link parameters #ex.constants.map(c => [#mi(c.name) $= #c.value$]).join([, ])].
Substituting into $T_0^#man.n$ (values rounded to four decimals):
#mitex("T = " + ex.T)

*Is it in #SE?* The last row is $mat(delim: "[", 0, 0, 0, 1)$, and for the rotation block
#mitex("R^T R = " + ex.RtR + ", \\qquad \\det R = " + ex.det)
The largest deviation of $R^T R$ from $I$ is #ex.orth_err, which is floating-point round-off#check(ex.member).

*The inverse.* Transpose the rotation block and rotate the negated position:
#mitex("R^T = " + ex.Rt + ", \\quad p = " + ex.p + ", \\quad -R^T p = " + ex.mRtp)
#mitex("T^{-1} = " + ex.T_inv)

*Does it undo $T$?* Multiplying the two matrices:
#mitex("T \\, T^{-1} = " + ex.product)
The largest deviation of $T T^(-1)$ and $T^(-1) T$ from $I_4$ is #ex.prod_err#check(ex.inverse_ok).
A general-purpose numeric matrix inversion of $T$ gives the same matrix, with the largest
difference between entries #ex.inv_err#check(ex.matches_general).

// ---------------------------------------------------------------- summary
= Summary

#set par(justify: false)
#table(
  columns: (auto, 1fr, 1.2fr, 1.4fr), stroke: none, inset: 6pt, align: (left, center, center, center),
  table.hline(),
  table.header([], [*Any group* $(G, compose)$], [*#SO*], [*#SE*]),
  table.hline(stroke: 0.5pt),
  [*Objects*], [elements of $G$], [$R$ with $R^T R = I$, $det R = 1$],
    [$mat(delim: "[", R, p; 0, 1)$, $R in SO$, $p in RR^3$],
  [*Operation*], [$a compose b in G$], [matrix product $R_1 R_2$], [matrix product $(R_1 R_2, thick R_1 p_2 + p_1)$],
  [*Identity*], [$e$], [$I_3$], [$I_4$],
  [*Inverse*], [$a^(-1)$], [$R^T$], [$(R^T, thick -R^T p)$],
  [*Commutative*], [not required], [no], [no],
  [*Degrees of freedom*], [—], [3], [6],
  [*Describes*], [—], [orientation], [orientation and position],
  table.hline(),
)
