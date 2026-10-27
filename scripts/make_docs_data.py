"""
make_docs_data.py -- compute the numbers shown on the documentation site.

Every chart and results table in docs/index.html is drawn from data that this
script computes with the project's own modules. Run it after changing the
code so the site stays truthful:

    python scripts/make_docs_data.py

It rewrites the block between the DOCS-DATA markers in docs/index.html.
"""

import json
import math
import os
import re
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from input_checker import DEMOS, check_input                      # noqa: E402
from thermo_models import CubicEOS, build_system, R               # noqa: E402
import equilibrium_solver as es                                   # noqa: E402


def r(v, d=5):
    return None if v is None or not np.isfinite(v) else round(float(v), d)


def example(name):
    with open(os.path.join(ROOT, "examples", name)) as fh:
        return fh.read()


data = {}

# ---------------------------------------------------------------- Pxy (PR)
sysPR = build_system(example("propane_butane_pr_pxy.txt"))
L, V = sysPR.equilibria[0].phases
d = es.vle_diagram(L, V, T=300.0, npts=51)
data["pxy"] = {"T": 300.0, "x": [r(v) for v in d.x1], "y": [r(v) for v in d.y1],
               "P": [r(v / 1e5) for v in d.var]}

# ----------------------------------------------- pure propane Psat, 3 EOS
psat = {}
for model in ("vdW", "SRK", "PR"):
    eos = CubicEOS(model, [369.8], [42.48e5], [0.152])
    P = 9e5
    for _ in range(200):
        lL = eos.ln_phi(300.0, P, [1.0], "liquid")[0]
        lV = eos.ln_phi(300.0, P, [1.0], "vapor")[0]
        P *= math.exp(lL - lV)
        if abs(lL - lV) < 1e-12:
            break
    psat[model] = r(P / 1e5, 3)
data["psat_propane_300K"] = psat

# -------------------------------------------------------- sugar-water
sysS = build_system(example("sugar_water.txt"))
res = es.solve(sysS)
boil = res[0].data
ice = res[1].data
suc = res[2].data
eut = [x for x in res if x.equilibrium == "Eutectic"][0].data
# extend the sucrose curve from x_sucrose = 0.02 down to the eutectic
Ls_, suc_phase = sysS.equilibria[2].phases
x_eu_suc = 1 - eut["x"]
extra_x, extra_T = [], []
for xs in np.geomspace(suc.x_solid[-1], x_eu_suc, 16)[1:]:
    extra_x.append(xs)
    extra_T.append(es.sle_temperature(Ls_, suc_phase, 101325.0, [1 - xs, xs]))
suc_x = np.concatenate([suc.x_solid, extra_x])
suc_T = np.concatenate([suc.T, extra_T])
data["sugar"] = {
    "boil": {"x": [r(v) for v in boil.x1], "T": [r(v, 3) for v in boil.var]},
    "ice": {"x": [r(v) for v in ice.x_solid], "T": [r(v, 3) for v in ice.T]},
    "sucrose": {"x": [r(1 - v, 6) for v in suc_x], "T": [r(v, 3) for v in suc_T]},
    "eutectic": {"x_water": r(eut["x"], 6), "T": r(eut["T"], 3)},
}
Ls, Vs = sysS.equilibria[0].phases
ice_phase = sysS.equilibria[1].phases[1]
data["sugar_point"] = {
    "Tb": r(es.bubble_T(Ls, Vs, 101325.0, [0.95, 0.05]).T, 3),
    "Tf": r(es.sle_temperature(Ls, ice_phase, 101325.0, [0.95, 0.05]), 3),
}

# ------------------------------------------------------- gravity column
sysG = build_system(example("gravity_column_gas_cap.txt"))
LG, VG = sysG.equilibria[0].phases
pt = es.flash_TP(LG, VG, 320.0, 25e5, [0.4, 0.6])
prof = es.column_profile(sysG, LG, VG, pt, npts=101)
contact = next(s for s in prof.notes if "contact" in s)
zc = float(re.search(r"contact at ([\d.]+) m", contact).group(1))
Pc = float(re.search(r"P = ([\d.]+) bar", contact).group(1))
data["gravity"] = {"z": [r(v, 2) for v in prof.position],
                   "P": [r(v / 1e5) for v in prof.P],
                   "x1": [r(v[0], 6) for v in prof.x],
                   "rho": [r(v, 2) for v in prof.density],
                   "phase": prof.phase, "contact_z": zc, "contact_P": Pc}

# -------------------------------------- centrifuge (bubble point x = 0.5)
sysC = build_system(example("ethanol_water_centrifuge.txt"))
LC, VC = sysC.equilibria[0].phases
pc = es.saturation(LC, VC, [0.5, 0.5], 0.0, T=350.0)
pc.phase = "two-phase"
profC = es.column_profile(sysC, LC, VC, pc, npts=11)
data["centrifuge"] = {"r": [r(v, 4) for v in profC.position],
                      "P": [r(v / 1e5, 3) for v in profC.P],
                      "x1": [r(v[0], 6) for v in profC.x]}

# ------------------------------------------------------ verification
b = es.bubble_P(L, V, 300.0, [0.4, 0.6])
dw = es.dew_P(L, V, 300.0, b.y, guess=b.P)
Pd = es.dew_P(L, V, 300.0, [0.4, 0.6]).P
fl = es.flash_TP(L, V, 300.0, 0.5 * (b.P + Pd), [0.4, 0.6])
mb = np.max(np.abs((1 - fl.beta) * fl.x + fl.beta * fl.y - np.array([0.4, 0.6])))
lf = np.max(np.abs(L.ln_f(300.0, 0.5 * (b.P + Pd), fl.x) - V.ln_f(300.0, 0.5 * (b.P + Pd), fl.y)))

sys100 = build_system(example("gravity_column_gas_cap.txt").replace("0, 500 m", "0, 100 m"))
L1, V1 = sys100.equilibria[0].phases
p100 = es.column_profile(sys100, L1, V1, es.flash_TP(L1, V1, 320.0, 25e5, [0.4, 0.6]), npts=21)
dz = np.diff(p100.position)
P_int = p100.P[0] - np.concatenate(
    [[0], np.cumsum(0.5 * (p100.density[1:] + p100.density[:-1]) * 9.80665 * dz)])
hydro = np.max(np.abs(P_int - p100.P)) / 1e5

data["verification"] = {
    "bubble_P_x04": r(b.P / 1e5, 5), "y_x04": r(b.y[0], 5),
    "roundtrip_err": float(f"{abs(dw.P / b.P - 1):.1e}"),
    "flash_beta": r(fl.beta, 4), "flash_mb": float(f"{mb:.1e}"),
    "flash_lnf": float(f"{lf:.1e}"),
    "hydro_err_bar": float(f"{hydro:.1e}"),
    "hydro_drop_bar": r((p100.P[0] - p100.P[-1]) / 1e5, 3),
}

# ------------------------------------------------- checker verdicts
verdicts = []
for fn in sorted(os.listdir(os.path.join(ROOT, "examples"))):
    rep = check_input(example(fn))
    verdicts.append({"file": fn, "status": rep.status})
data["verdicts"] = verdicts

# ------------------------------------------------------ write into the page
page = os.path.join(ROOT, "docs", "index.html")
html = open(page, encoding="utf-8").read()
blob = json.dumps(data, separators=(",", ":"))
new = re.sub(r"(<!-- DOCS-DATA:START -->).*?(<!-- DOCS-DATA:END -->)",
             lambda m: (m.group(1) + '\n<script id="docs-data" type="application/json">'
                        + blob + "</script>\n" + m.group(2)),
             html, flags=re.S)
if new == html and "DOCS-DATA:START" not in html:
    sys.exit("docs/index.html has no DOCS-DATA markers")
open(page, "w", encoding="utf-8").write(new)
print(f"wrote {len(blob) / 1024:.1f} kB of data into docs/index.html")
print(json.dumps({k: v for k, v in data.items()
                  if k in ("psat_propane_300K", "sugar_point", "verification", "verdicts")},
                 indent=1))
print("eutectic", data["sugar"]["eutectic"], "contact", zc, Pc)
print("centrifuge", data["centrifuge"]["P"][0], data["centrifuge"]["P"][-1],
      data["centrifuge"]["x1"][-1])
