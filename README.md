# Phase Equilibrium Calculator

A Python tool that takes whatever you know about a mixture, works out which phase equilibria that information can decide, then calculates and plots them: vapor–liquid and solid–liquid equilibrium, phase diagrams with azeotropes and eutectics, and pressure and composition profiles in gravity and centrifugal fields.

Course project for Chemical Engineering Thermodynamics. Everything is derived from the fundamental equation in its Helmholtz form,

$$dA = -S\,dT - P\,dV + \sum_i \mu_i\,dN_i \quad [\,+\ \Phi\,dm \text{ in a field}\,]$$

and the criterion that each component's chemical potential is the same in coexisting phases.

**Documentation:** https://fang099.github.io/phase-equilibrium-calculator/ (theory, input format, worked examples, plots and verification).

## The idea: the information you have decides the route

There is no fixed list of required inputs. Each kind of data opens a way to calculate a chemical potential, and the checker combines whatever is given:

| If you know | The tool can calculate |
|---|---|
| `Tc`, `Pc`, `omega` (or vdW `a`, `b`) | both phases from one cubic equation of state (vdW, SRK, PR) |
| Antoine constants, a typed `Psat(T)` relation, or `Tb` with `dHvap` | vapor–liquid equilibrium through the pure-liquid fugacity |
| activity-model parameters (Margules, van Laar, NRTL) | non-ideal liquids, including azeotropes |
| `Tm` with `dHfus` (and optionally `dCp_fus`) | solid–liquid equilibrium: crystallization curves and the eutectic |
| `nonvolatile = yes` | keeps a component out of the vapor |
| molar masses (plus liquid volumes without an EOS) | profiles in a gravity or centrifugal field |

It reports every equilibrium it can decide, names the calculation, and lists exactly what is missing for the rest.

## Quick start

Requires Python 3.9+ with NumPy, SciPy and Matplotlib.

```bash
git clone https://github.com/Fang099/phase-equilibrium-calculator.git
cd phase-equilibrium-calculator
pip install -r requirements.txt

python input_checker.py examples/benzene_naphthalene.txt      # is the input sufficient?
python equilibrium_solver.py examples/benzene_naphthalene.txt # check and solve
python plotter.py examples/benzene_naphthalene.txt            # solve and save plots to plots/
```

Describe a system in plain text, one known quantity per line, with units. This input has only boiling points, melting points and latent heats, and gives the complete phase diagram at 1 atm, eutectic included:

```text
components = benzene, naphthalene
P = 1 atm
sweep = x
benzene.Tb = 353.24 K
benzene.dHvap = 30.72 kJ/mol
benzene.Tm = 278.68 K
benzene.dHfus = 9.87 kJ/mol
naphthalene.Tb = 491.14 K
naphthalene.dHvap = 43.2 kJ/mol
naphthalene.Tm = 353.43 K
naphthalene.dHfus = 19.01 kJ/mol
```

Add an activity model with, for example:

```text
activity = vanlaar
A(benzene, naphthalene) = 0.10
A(naphthalene, benzene) = 0.08
```

`python input_checker.py --vocab` lists every recognized quantity.

From Python or Colab:

```python
from thermo_models import build_system
from equilibrium_solver import solve, print_results
from plotter import plot_results

system = build_system(open("examples/benzene_naphthalene.txt").read())
results = solve(system)
print_results(system, results)
figures = plot_results(system, results)   # {name: matplotlib Figure}
```

## Modules

| Module | Purpose | Status |
|---|---|---|
| `input_checker.py` | Parses free-form input, converts units, applies the phase rule (`F = n + 2 + m − π`, or Duhem's theorem with a known feed), and decides which equilibria can be calculated and what is missing | Built |
| `thermo_models.py` | Fugacity models for each phase: vdW, SRK and PR as a residual Helmholtz energy; liquid solution with Margules, van Laar or NRTL; ideal gas; pure solid from latent heats; gravity and centrifugal potentials | Built |
| `equilibrium_solver.py` | Bubble and dew points, flash, Pxy and Txy diagrams with azeotropes, solid–liquid curves and eutectic, field profiles with gas–liquid contacts | Built |
| `plotter.py` | Pxy/Txy and xy diagrams, complete T–x diagrams, composition charts, field profiles, 3D P–T–x and Gibbs-energy-of-mixing surfaces | Built |
| Validation | Comparison with experimental data | Next |
| Interface | Notebook or web front end | Planned |

## Tests

Each module runs its own verification when executed directly:

```bash
python thermo_models.py        # EOS derivatives, activity models, latent-heat route
python equilibrium_solver.py   # round trips, flash balances, azeotrope, eutectic, hydrostatics
python input_checker.py        # verdicts on the demo inputs
```

## Repository layout

```text
input_checker.py           stage 1: is the input sufficient?
thermo_models.py           stage 2: fugacity models for each phase
equilibrium_solver.py      stage 3: solve the equilibrium conditions
plotter.py                 stage 4: draw the results
examples/                  input files, one per demo case
docs/                      documentation website (GitHub Pages)
scripts/make_docs_data.py  recomputes the numbers, charts and images on the website
```

## Documentation website

The site in `docs/` is served by GitHub Pages from the `main` branch, folder `/docs`. All numbers, charts and images on it are computed by the tool. After changing the code, run `python scripts/make_docs_data.py` to refresh them.

## Notes

Property values in `examples/` are typical literature values for demonstration. Check them against a source you can cite before relying on the results.
