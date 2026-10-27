# Phase Equilibrium Calculator

A Python tool that works out whether the data you have are enough to calculate the phase equilibrium of a mixture, then calculates it: vapor–liquid and solid–liquid equilibrium, phase diagrams, and pressure and composition profiles in gravity and centrifugal fields.

Course project for Chemical Engineering Thermodynamics. Reference text: Tester and Modell, *Thermodynamics and Its Applications*. Everything is derived from the fundamental equation in its Helmholtz form,

$$dA = -S\,dT - P\,dV + \sum_i \mu_i\,dN_i \quad [\,+\ \Phi\,dm \text{ in a field}\,]$$

and the criterion that each component's chemical potential is the same in coexisting phases.

**Documentation:** https://fang099.github.io/phase-equilibrium-calculator/ (theory, input format, worked examples and verification).

## Quick start

Requires Python 3.9+ with NumPy and SciPy.

```bash
git clone https://github.com/Fang099/phase-equilibrium-calculator.git
cd phase-equilibrium-calculator
pip install -r requirements.txt

python input_checker.py examples/sugar_water.txt        # is the input sufficient?
python equilibrium_solver.py examples/sugar_water.txt   # check and solve
```

Describe a system in plain text, one known quantity per line, with units:

```text
components = water, sucrose
P = 1 atm
sweep = x
water.Tb = 373.15 K
water.dHvap = 40.65 kJ/mol
water.Tm = 273.15 K
water.dHfus = 6.01 kJ/mol
sucrose.nonvolatile = yes
sucrose.Tm = 459 K
sucrose.dHfus = 46 kJ/mol
```

`python input_checker.py --vocab` lists every recognized quantity.

From Python or Colab:

```python
from thermo_models import build_system
from equilibrium_solver import solve, print_results

system = build_system(open("examples/sugar_water.txt").read())
print_results(system, solve(system))
```

## What it does now

| Module | Purpose | Status |
|---|---|---|
| `input_checker.py` | Parses free-form input, converts units, applies the phase rule (`F = n + 2 + m − π`, or Duhem's theorem with a known feed), and decides which equilibria can be calculated and what is missing | Built |
| `thermo_models.py` | Fugacity models for each phase: vdW, SRK and PR as a residual Helmholtz energy; ideal solution; ideal gas; pure solid from latent heats; gravity and centrifugal potentials | Built |
| `equilibrium_solver.py` | Bubble and dew points, flash, Pxy and Txy diagrams, solid–liquid curves and eutectic, field profiles with gas–liquid contacts | Built |
| Plotting | Diagrams, field profiles, 3D P–T–x surfaces | Next |
| Validation | Comparison with experimental data | Planned |
| Interface | Notebook or web front end | Planned |

## Tests

Each module runs its own verification when executed directly:

```bash
python thermo_models.py        # analytic vs numerical derivatives, identities, pure-propane Psat
python equilibrium_solver.py   # round trips, flash balances, hand calculations, hydrostatics
python input_checker.py        # verdicts on the demo inputs
```

## Repository layout

```text
input_checker.py         stage 1: is the input sufficient?
thermo_models.py         stage 2: fugacity models for each phase
equilibrium_solver.py    stage 3: solve the equilibrium conditions
examples/                input files, one per demo case
docs/                    documentation website (GitHub Pages)
scripts/make_docs_data.py  recomputes the numbers and charts on the website
```

## Documentation website

The site in `docs/` is a single static page. To publish it, go to the repository's **Settings → Pages**, choose **Deploy from a branch**, branch `main`, folder `/docs`.

All numbers and charts on the site are computed by the tool. After changing the code, run `python scripts/make_docs_data.py` to refresh them.

## Notes

Property values in `examples/` are typical literature values for demonstration. Take values for real work from the appendix of Tester and Modell or another cited source.
