# Kinematics walkthroughs: manipulators/<name>.toml -> out/<name>.pdf

set shell := ["bash", "-euo", "pipefail", "-c"]

python := justfile_directory() / ".venv/bin/python"

# list recipes and available manipulators
default:
    @just --list
    @echo "manipulators: $(ls manipulators/*.toml | xargs -n1 basename | sed 's/\.toml$//' | grep -v '^_' | tr '\n' ' ')"

# one-time: create the virtualenv with sympy, numpy, matplotlib
setup:
    python3 -m venv .venv
    .venv/bin/pip install -q sympy numpy matplotlib

# symbolic + numeric work: spec -> build/<name>/ (JSON + SVG figures)
generate name:
    "{{python}}" gen_walkthrough.py "manipulators/{{name}}.toml" -o "build/{{name}}"

# quick schematic only (no symbolic simplification): build/<name>/schematic.svg
schematic name:
    "{{python}}" gen_walkthrough.py "manipulators/{{name}}.toml" -o "build/{{name}}" --schematic-only

# full pipeline: spec -> PDF
walkthrough name: (generate name)
    @mkdir -p out
    typst compile --input "data=build/{{name}}/walkthrough.json" walkthrough.typ "out/{{name}}.pdf"
    @echo '-> out/{{name}}.pdf'

# re-render the PDF on template edits (re-run `just generate` after spec edits)
watch name: (generate name)
    @mkdir -p out
    typst watch --input "data=build/{{name}}/walkthrough.json" walkthrough.typ "out/{{name}}.pdf"

# interactive 3-D app with every manipulator: out/app.html (open it in a browser)
app:
    "{{python}}" gen_app.py -o out/app.html

# build every manipulator spec
all:
    for f in manipulators/[!_]*.toml; do just walkthrough "$(basename "$f" .toml)"; done

# start a new spec from the commented template
new name:
    test ! -e "manipulators/{{name}}.toml"
    cp manipulators/_template.toml "manipulators/{{name}}.toml"
    @echo 'edit manipulators/{{name}}.toml, then: just walkthrough {{name}}'

clean:
    rm -rf build out
