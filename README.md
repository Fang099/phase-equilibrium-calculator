# Phase Equilibrium Calculator

A Python program that calculates and plots the phase equilibrium of mixtures. It takes whatever you know about a mixture, works out which equilibria that information can decide, and solves them. A website runs the same Python code in the browser as an interactive calculator.

**Website:** https://fang099.github.io/phase-equilibrium-calculator/

- **Calculator:** pick any two of 24 compounds and a model, and plot their vapor–liquid equilibrium between any two variables.
- **Concepts:** chemical potential, fugacity, activity coefficients and the models, in plain language.
- **How it works:** one calculation followed through the code.
- **Examples** and **Reference:** worked cases, the input format, the compound data and its sources.

Course project for Chemical Engineering Thermodynamics. Everything is derived from the fundamental equation in its Helmholtz form,

$$dA = -S\,dT - P\,dV + \sum_i \mu_i\,dN_i \quad [\,+\ \Phi\,dm \text{ in a field}\,]$$

and the criterion that each compound's chemical potential (equivalently, its fugacity) is the same in coexisting phases.

## Quick start

Requires Python 3.9+ with NumPy and Matplotlib. SciPy is used if installed and is not required.

```bash
git clone https://github.com/Fang099/phase-equilibrium-calculator.git
cd phase-equilibrium-calculator
pip install -r requirements.txt

python input_checker.py examples/database_ethanol_water_nrtl.txt      # is the input sufficient?
python equilibrium_solver.py examples/database_ethanol_water_vanlaar.txt  # van Laar with A(T) from the database
python equilibrium_solver.py examples/database_ethanol_water_nrtl.txt # check and solve
python plotter.py examples/database_ethanol_water_nrtl.txt            # solve and save plots to plots/

python scripts/dev_server.py     # the website and calculator on http://localhost:8000/, no internet needed
```

An input file lists what is known, one quantity per line, with units. With the built-in database, two compound names are enough:

```text
components = ethanol, water
database = builtin
activity = nrtl
P = 1 atm
sweep = x
```

Without the database, give the properties yourself. This input has only boiling points, melting points and latent heats, and gives the complete phase diagram at 1 atm, eutectic included:

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

`python input_checker.py --vocab` lists every recognized quantity. Values typed in the input always take precedence over the database.

## The information you have decides the route

| If you know | The program can calculate |
|---|---|
| `Tc`, `Pc`, `omega` (or vdW `a`, `b`) | both phases from one cubic equation of state (vdW, SRK, PR) |
| Antoine constants, a typed `Psat(T)` relation, or `Tb` with `dHvap` | vapor–liquid equilibrium through the pure-liquid fugacity |
| activity-model parameters (Margules, van Laar, NRTL) | non-ideal liquids: azeotropes, and the three-phase line when the model predicts two liquid layers |
| `Tm` with `dHfus` (and optionally `dCp_fus`) | solid–liquid equilibrium: crystallization curves and the eutectic |
| `nonvolatile = yes` | keeps a compound out of the vapor |
| molar masses (plus liquid volumes without an EOS) | profiles in a gravity or centrifugal field |

The checker reports every equilibrium it can decide, names the calculation, and lists exactly what is missing for the rest.

## Modules

| Module | Purpose |
|---|---|
| `input_checker.py` | Parses free-form input, converts units, applies the phase rule (`F = n + 2 + m − π`, or Duhem's theorem with a known feed), and decides which equilibria can be calculated and what is missing |
| `compound_db.py` | The built-in database (`data/compounds.json`) and `database = builtin` |
| `thermo_models.py` | Fugacity models for each phase: vdW, SRK and PR as a residual Helmholtz energy; liquid solution with Margules, van Laar or NRTL; ideal gas; pure solid from latent heats; gravity and centrifugal potentials |
| `equilibrium_solver.py` | Bubble and dew points, flash, Pxy/Txy and P–T diagrams, azeotropes, liquid stability and the three-phase line, activity coefficients and fugacities along a diagram, solid–liquid curves and eutectic, field profiles |
| `plotter.py` | Matplotlib figures: diagrams, fugacity and activity plots, complete T–x diagrams, field profiles, 3D surfaces |
| `web_api.py` | JSON interface used by the website's calculator; also usable from Python |
| `numerics.py` | Root finders used when SciPy is not installed (as in the browser) |

Next stages: systematic comparison with experimental data, liquid–liquid equilibrium as a calculation of its own, ternary diagrams.

## Compound data

`data/compounds.json` holds 24 compounds, NRTL parameters for 66 pairs, temperature-dependent van Laar parameters for 65 pairs and Peng–Robinson k₁₂ for 23 pairs. It is generated by `scripts/build_compound_db.py` and `scripts/fit_vanlaar.py`, never typed by hand, and every value carries a tag naming its source:

- **Critical constants and acentric factors:** Horstmann, Jabłoniec, Krafczyk, Fischer, Gmehling, *Fluid Phase Equilibria* 227 (2005) 157–164 (PSRK revision IV).
- **Antoine constants:** Poling, Prausnitz, O'Connell, *The Properties of Gases and Liquids*, 5th ed., 2001, Appendix A.
- **NRTL parameters and Peng–Robinson k₁₂:** ChemSep interaction parameter data (H. Kooijman and R. Taylor, Artistic License 2.0).
- **van Laar parameters** `A_ij(T) = a_ij + b_ij/T`: regressed by `scripts/fit_vanlaar.py` from the DECHEMA NRTL parameters of the same pair (ln γ over all compositions, between the two boiling points). No open van Laar table covers these pairs; the fits agree with commonly quoted constants (see the website's Reference page).
- All read from the data files of the open-source `chemicals` and `thermo` libraries (C. Bell and contributors, MIT License).

To add compounds, add rows to `data/compound_list.csv` (name, formula, CAS number) and run the script again, then `python scripts/fit_vanlaar.py`; it accepts any compound those tables cover.

## Tests

Each module verifies itself when run directly:

```bash
python thermo_models.py        # EOS derivatives, activity models, latent-heat route
python equilibrium_solver.py   # round trips, flash balances, azeotrope, eutectic, hydrostatics
python web_api.py              # calculator requests: azeotropes, equal fugacities, three-phase line
python input_checker.py        # verdicts on the demo inputs
```

## Repository layout

```text
input_checker.py, compound_db.py, thermo_models.py,
equilibrium_solver.py, plotter.py, web_api.py, numerics.py    the program
data/                      compound database and the list it is built from
examples/                  input files, one per demo case
site/                      page sources and the shared layout of the website
docs/                      the built website (GitHub Pages serves this folder)
scripts/build_site.py      rebuilds docs/: pages, numbers, charts, images, and the copy of the code the browser runs
scripts/precompute_surfaces.py   solves the 3D surfaces in advance into docs/surfaces/
scripts/build_compound_db.py   regenerates data/compounds.json from the source tables
scripts/dev_server.py      local preview of the website with a local calculator engine
```

After changing any Python module or page, run `python scripts/build_site.py`, then `python scripts/precompute_surfaces.py` to bring the precomputed 3D surfaces (`docs/surfaces/`) up to date. No number on the site is typed by hand. On GitHub, the workflow `.github/workflows/surfaces.yml` runs both automatically on every push to `main` that changes the code, the data or the pages, and commits the result.

The first version of the project (before the calculator and database) is tagged `v0.1`.

## Notes

Database values are reproduced from the cited sources for teaching use. Check them against the original publication before relying on them. Property values typed in `examples/` are typical literature values for demonstration.
