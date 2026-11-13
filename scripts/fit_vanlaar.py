"""
fit_vanlaar.py -- add temperature-dependent van Laar parameters to the
compound database.

    python scripts/fit_vanlaar.py            # updates data/compounds.json

No open table of van Laar constants covers these pairs, so the parameters are
regressed here from the DECHEMA NRTL parameters already in the database
(ChemSep file 'DECHEMA NRTL at P = 1 atm'), which were themselves fitted to
measured isobaric vapor-liquid data. For each pair with NRTL parameters:

  * activity coefficients are generated with NRTL on a grid of liquid
    compositions (x1 = 0.01 ... 0.99) and temperatures spanning the two
    normal boiling points (at least 30 K wide), the range of the data;
  * the four constants of

        A_ij(T) = a_ij + b_ij / T          (T in K)

    are fitted to ln(gamma_1) and ln(gamma_2) at every grid point at once
    (least squares). A_ij(T) is ln(gamma_i) at infinite dilution in j, and
    the 1/T form is the temperature dependence of the original van Laar
    theory (constant excess energy);
  * the fit is kept with its root-mean-square and largest error in
    ln(gamma), and its temperature range. A pair whose two infinite-dilution
    ln(gamma) have opposite signs cannot be described by van Laar and gets no
    entry.

Each entry carries the source tag 'vanlaar_fit'. Re-run after
scripts/build_compound_db.py.
"""

import json
import os
import sys
from datetime import date

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from thermo_models import NRTL                                        # noqa: E402

try:
    from scipy.optimize import least_squares
except ImportError:                                                     # pragma: no cover
    least_squares = None

DB = os.path.join(ROOT, "data", "compounds.json")
SOURCE = ("van Laar parameters A_ij(T) = a_ij + b_ij/T regressed in this project "
          "(scripts/fit_vanlaar.py) to ln(gamma) from the DECHEMA NRTL parameters of the "
          "same pair (ChemSep file 'DECHEMA NRTL at P = 1 atm', see chemsep_nrtl), over "
          "x1 = 0.01 to 0.99 and the temperature range given with each pair.")
X = np.linspace(0.01, 0.99, 49)


def vl_ln_gamma(A12, A21, x1):
    x2 = 1 - x1
    d = A12 * x1 + A21 * x2
    return A12 * (A21 * x2 / d) ** 2, A21 * (A12 * x1 / d) ** 2


def gauss_newton(res, p0, iters=200):
    """Small fallback when SciPy is not installed."""
    p = np.array(p0, float)
    for _ in range(iters):
        r = res(p)
        J = np.empty((len(r), len(p)))
        for k in range(len(p)):
            h = 1e-6 * max(1.0, abs(p[k]))
            q = p.copy(); q[k] += h
            J[:, k] = (res(q) - r) / h
        step = np.linalg.lstsq(J, -r, rcond=None)[0]
        lam = 1.0
        while lam > 1e-4 and np.sum(res(p + lam * step) ** 2) > np.sum(r ** 2):
            lam /= 2
        p = p + lam * step
        if np.max(np.abs(lam * step)) < 1e-12:
            break
    return p


def fit_pair(nrtl, Tb_a, Tb_b):
    lo, hi = sorted((Tb_a, Tb_b))
    if hi - lo < 30:
        mid = 0.5 * (lo + hi)
        lo, hi = mid - 15, mid + 15
    Ts = np.linspace(lo, hi, 7)
    R = 8.314462618
    m = NRTL(np.array([[0, nrtl["b_ab"] * R], [nrtl["b_ba"] * R, 0]]),
             np.array([[0, nrtl["alpha"]], [nrtl["alpha"], 0]]))
    target = np.array([[m.ln_gamma(T, [x, 1 - x]) for x in X] for T in Ts])   # (T, x, 2)
    inf12 = np.array([m.ln_gamma(T, [0.0, 1.0])[0] for T in Ts])
    inf21 = np.array([m.ln_gamma(T, [1.0, 0.0])[1] for T in Ts])
    if np.any(inf12 * inf21 <= 0):
        return None, "infinite-dilution ln(gamma) of opposite sign: outside the van Laar form"
    s = np.sign(inf12[0])

    def res(p):
        a12, b12, a21, b21 = p
        out = []
        for k, T in enumerate(Ts):
            A12, A21 = a12 + b12 / T, a21 + b21 / T
            if A12 * s <= 0 or A21 * s <= 0:
                return np.full(2 * len(Ts) * len(X), 1e3)
            g1, g2 = vl_ln_gamma(A12, A21, X)
            out += [g1 - target[k, :, 0], g2 - target[k, :, 1]]
        return np.concatenate(out)

    # start from the infinite-dilution values, with their own 1/T slope
    c12 = np.polyfit(1 / Ts, inf12, 1)
    c21 = np.polyfit(1 / Ts, inf21, 1)
    p0 = [c12[1], c12[0], c21[1], c21[0]]
    if least_squares is not None:
        p = least_squares(res, p0, method="lm", xtol=1e-12, ftol=1e-12).x
    else:
        p = gauss_newton(res, p0)
    r = res(p)
    return {"a12": round(float(p[0]), 5), "b12": round(float(p[1]), 3),
            "a21": round(float(p[2]), 5), "b21": round(float(p[3]), 3),
            "Tmin": round(float(Ts[0]), 1), "Tmax": round(float(Ts[-1]), 1),
            "rms_lngamma": round(float(np.sqrt(np.mean(r ** 2))), 4),
            "max_lngamma": round(float(np.max(np.abs(r))), 4),
            "source": "vanlaar_fit"}, None


def main():
    with open(DB, encoding="utf-8") as fh:
        db = json.load(fh)
    Tb = {c["id"]: c["Tb"]["value"] for c in db["compounds"]}
    n_ok, skipped = 0, []
    for p in db["pairs"]:
        p.pop("vanlaar", None)
        if "nrtl" not in p:
            continue
        entry, why = fit_pair(p["nrtl"], Tb[p["a"]], Tb[p["b"]])
        if entry is None:
            skipped.append(f"{p['a']} / {p['b']}: {why}")
            continue
        p["vanlaar"] = entry
        n_ok += 1
        print(f"{p['a']:>18s} / {p['b']:<18s} A12 = {entry['a12']:+.4f} {entry['b12']:+9.2f}/T   "
              f"A21 = {entry['a21']:+.4f} {entry['b21']:+9.2f}/T   "
              f"{entry['Tmin']:.0f}-{entry['Tmax']:.0f} K   rms {entry['rms_lngamma']:.4f}")
    db["meta"]["sources"]["vanlaar_fit"] = SOURCE
    db["meta"]["vanlaar_fitted"] = str(date.today())
    with open(DB, "w", encoding="utf-8") as fh:
        json.dump(db, fh, indent=1)
    print(f"\n{n_ok} pairs fitted; {len(skipped)} without a van Laar entry:")
    for s_ in skipped:
        print("  " + s_)


if __name__ == "__main__":
    main()
