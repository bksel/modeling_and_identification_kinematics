# DH kinematics walkthrough

Generates an A4 PDF walkthrough of a manipulator's kinematics from its
Denavit–Hartenberg table. Each document contains:

1. A kinematic schematic drawn automatically from the DH table (and your own figure next to it, if given)
2. The DH table and the general DH transformation matrix
3. Forward kinematics step by step: every `A_i` and the cumulative `T_0^i`, simplified with
   the usual `c_12 = cos(q1 + q2)` shorthand
4. The end-effector position `p(q)` and a plot of the frames at a sample configuration
5. The position Jacobian `J(q) = ∂p/∂q`
6. Inverse kinematics:
   - your closed-form solution (optional), typeset and checked against forward kinematics on 500 random targets
   - a numeric solution by damped least squares, with an iteration table and a plot

## Requirements

- Python 3.11 or newer (for `tomllib`)
- [`just`](https://github.com/casey/just)
- [`typst`](https://github.com/typst/typst) 0.14
- Internet access on the first build, when Typst downloads the `mitex` package

## Setup

```sh
just setup    # creates .venv with sympy, numpy and matplotlib
```

## Usage

Run all commands from the repository root.

```sh
just                      # list recipes and the available manipulators
just walkthrough planar2  # manipulators/planar2.toml -> out/planar2.pdf
just all                  # build every spec in manipulators/
just schematic trr        # draw only the schematic (fast) -> build/trr/schematic.svg
just watch trr            # re-render the PDF whenever walkthrough.typ changes
just new mymanip          # start a new spec from the commented template
just clean                # remove build/ and out/
```

`just walkthrough <name>` runs in two stages. You can also run each one on its own:

1. `just generate <name>`: the Python stage does the symbolic and numeric work and writes
   `build/<name>/walkthrough.json` plus SVG figures.
2. `typst compile`: renders `walkthrough.typ` with that JSON into `out/<name>.pdf`.

A useful loop while writing a new DH table:

1. Edit the spec.
2. Run `just schematic <name>` and compare the drawing with your sketch.
3. When it matches, run `just walkthrough <name>`.

## Writing a manipulator spec

Specs are TOML files in `manipulators/`. Every value is parsed by sympy, so you can use
symbols (`q1`, `a1`), numbers, `pi`, `sqrt()`, `atan2()`, … The name `deg` means `pi/180`, so
`"30*deg"` works. `x`, `y` and `z` are reserved for the end-effector position. Avoid the
single-letter sympy names `S`, `E`, `I`, `N`, `O` and `Q`.

A minimal spec, a planar 2R arm:

```toml
name = "Planar 2R pendulum"
description = "Shown under the title."
# figure = "figures/planar2.png"   # optional: your own sketch, relative to this file

joints = ["q1", "q2"]              # joint variables, in DH-row order

[constants]                        # numeric link values: used only for figures and numeric IK
a1 = 1.0
a2 = 0.8

[sample]                           # configuration for the frame plot; its position is the IK target
q1 = "30*deg"
q2 = "45*deg"

[[dh]]                             # one block per joint: standard DH (theta, d, a, alpha)
theta = "q1"
d = "0"
a = "a1"
alpha = "0"

[[dh]]
theta = "q2"
d = "0"
a = "a2"
alpha = "0"
```

A joint is treated as prismatic when its variable appears in `d`, and as revolute when it
appears in `theta`.

### Schematic options

```toml
[schematic]
view = "auto"          # "auto" (top for planar arms, oblique otherwise), "top" or "oblique"

[schematic.pose]       # pose used only for the drawing (defaults to [sample])
q1 = "90*deg"
```

The oblique view matches the exercise sheets: X₀ points toward the viewer, Y₀ to the right and
Z₀ up. Pick a pose where the arm is not foreshortened. Joint symbols:

| Symbol | Meaning |
|---|---|
| crossed box | translational (prismatic) joint |
| `]\|[` bearing | revolute joint, axis in the drawing plane |
| ellipse with dashed axis | revolute joint, axis pointing toward X₀ (oblique view) |
| circle with angle arc | revolute joint, axis toward the viewer (top view) |

### Inverse kinematics

The numeric solution always runs. Optional settings:

```toml
[ik]
damping = 0.05

[ik.initial_guess]
q1 = "90*deg"
```

To include a closed-form solution, list any auxiliary quantities and then one block per solution
branch. Both are evaluated in order, so later entries can use earlier ones, including joint values:

```toml
[ik.intermediates]
D = "(x**2 + y**2 - a1**2 - a2**2) / (2*a1*a2)"

[[ik.solution]]
name = "Elbow branch q2 > 0"
q2 = "atan2(sqrt(1 - D**2), D)"
q1 = "atan2(y, x) - atan2(a2*sin(q2), a1 + a2*cos(q2))"
```

The document reports, for each branch, the largest position error over 500 random reachable
targets. A branch passes when that error is below 1e-8. A failure means the formula is wrong.

For redundant arms, where there are more joints than position coordinates (the planar 3R and 4R),
the document says that the solutions form a family, and only the numeric solution is shown.

## Included examples

| Spec | Manipulator |
|---|---|
| `planar2` | planar 2R pendulum, with closed-form IK |
| `planar3` | planar 3R pendulum (redundant) |
| `planar4` | planar 4R pendulum (redundant) |
| `elbow3r` | 3R elbow arm: base rotation plus two parallel horizontal joints, with closed-form IK |
| `trr` | TRR from exercise list 5: translational q1, then two revolute joints |

## Repository layout

```
Justfile             build recipes
kinematics_lib.py    SO(3)/SE(3) helpers and DH transforms (sympy)
gen_walkthrough.py   spec -> JSON + figures (symbolic FK, Jacobian, IK, plots)
schematic.py         kinematic schematic drawn from the DH frames
walkthrough.typ      Typst layout of the A4 document
manipulators/        specs; _template.toml documents every option
build/, out/         generated (ignored by git)
```

## Limitations

- **Orientation is out of scope.** Only the end-effector position is analysed: no RPY or Euler angles.
- **Closed-form IK is not derived automatically.** You write it; the tool typesets and verifies it.
- **Joints at the same point overlap in the schematic.** This happens when d = a = 0, for example in a spherical wrist.
