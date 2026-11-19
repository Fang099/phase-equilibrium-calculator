"""Compare the models with measured and reference data.

    python scripts/validate_real_data.py            writes validation/results/*.csv
                                                     and prints a summary

Data used (all in validation/):
  reference_saturation.csv   vapor pressure and coexisting densities of ten
                             fluids from their reference equations of state
                             (IAPWS-95 for water; Span-Wagner for CO2; Lemmon
                             et al. for propane, ...; sources in
                             reference_eos_sources.json). These equations
                             reproduce the measured data to within their
                             experimental uncertainty and are the basis of the
                             NIST Chemistry WebBook fluid tables.
  reference_field.csv        the same reference equations integrated through
                             the field (hydrostatic columns, the coexistence
                             band, a CO2 centrifuge).
  vle_*.csv                  measured binary VLE from the literature, as
                             archived in the NIST TRC ThermoML archive.
"""
import csv
import math
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import pure_fluid as pf                                         # noqa: E402
import web_api                                                  # noqa: E402
from thermo_models import build_system                          # noqa: E402
from equilibrium_solver import bubble_T, bubble_P, vle_diagram  # noqa: E402

VAL = os.path.join(ROOT, "validation")
OUT = os.path.join(VAL, "results")
R = 8.314462618
G0 = 9.81


def read(name):
    with open(os.path.join(VAL, name)) as fh:
        rows = [r for r in fh if not r.startswith("#")]
    return list(csv.DictReader(rows))


def write(name, rows):
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, name), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in r.items()})
    return os.path.join(OUT, name)


# --------------------------------------------------------------------------
# A. pure-fluid saturation against the reference equations of state
# --------------------------------------------------------------------------
def pure_saturation():
    ref = read("reference_saturation.csv")
    points, summary = [], []
    for comp in dict.fromkeys(r["compound"] for r in ref):
        rows = [r for r in ref if r["compound"] == comp]
        for model in ("vdW", "SRK", "PR"):
            eos, M, _ = pf.make_eos(comp, model)
            dP, dL, dV = [], [], []
            for r in rows:
                T = float(r["T_K"])
                try:
                    P, VL, VV = pf.saturation_pressure(eos, T)
                except ValueError:
                    continue
                p_ref, l_ref, v_ref = float(r["Psat_Pa"]), float(r["rhoL_kg_m3"]), float(r["rhoV_kg_m3"])
                eP, eL, eV = (100 * (P / p_ref - 1), 100 * (M / VL / l_ref - 1), 100 * (M / VV / v_ref - 1))
                dP.append(eP); dL.append(eL); dV.append(eV)
                points.append({"compound": comp, "model": model, "T_K": T, "Psat_ref_bar": p_ref / 1e5,
                               "Psat_model_bar": P / 1e5, "rhoL_ref": l_ref, "rhoL_model": M / VL,
                               "rhoV_ref": v_ref, "rhoV_model": M / VV,
                               "dev_Psat_pct": eP, "dev_rhoL_pct": eL, "dev_rhoV_pct": eV})
            aad = lambda a: float(np.mean(np.abs(a))) if a else float("nan")
            summary.append({"compound": comp, "model": model, "n": len(dP), "T_min_K": float(rows[0]["T_K"]),
                            "T_max_K": float(rows[-1]["T_K"]), "AAD_Psat_pct": aad(dP),
                            "AAD_rhoL_pct": aad(dL), "AAD_rhoV_pct": aad(dV)})
    write("pure_saturation_points.csv", points)
    write("pure_saturation_summary.csv", summary)
    return summary


# --------------------------------------------------------------------------
# B. binary VLE at 1 atm with NRTL (and the ideal solution) against measurements
# --------------------------------------------------------------------------
NRTL_SETS = [("vle_ethanol_water_101kPa_Kamihama2012.csv", "ethanol", "water"),
             ("vle_acetone_ethanol_101kPa_KuTu2005.csv", "acetone", "ethanol"),
             ("vle_chloroform_methanol_101kPa_Li2014.csv", "chloroform", "methanol")]


def _phases(req, first):
    s = build_system(web_api.build_input(req))
    eq = [e for e in s.equilibria if e.phase_names[1] == "V"][0]
    L, V = eq.phases
    flip = s.names[0] != first           # the input may list the pair in the other order
    return L, V, flip


def _exp_azeotrope(x, y, T):
    d = np.array(y) - np.array(x)
    for i in range(len(d) - 1):
        if 0 < x[i] and x[i + 1] < 1 and d[i] * d[i + 1] < 0:
            f = d[i] / (d[i] - d[i + 1])
            return x[i] + f * (x[i + 1] - x[i]), T[i] + f * (T[i + 1] - T[i])
    return None


def binary_activity():
    points, summary = [], []
    for fname, a, b in NRTL_SETS:
        data = read(fname)
        P = float(data[0]["P_kPa"]) * 1e3
        xs = [float(r["x1"]) for r in data]
        ys = [float(r["y1"]) for r in data]
        Ts = [float(r["T_K"]) for r in data]
        az_exp = _exp_azeotrope(xs, ys, Ts)
        for model in ("nrtl", "ideal"):
            L, V, flip = _phases({"components": [a, b], "model": model, "hold": "P", "P": P / 1e5}, a)
            dT, dy = [], []
            for x, y, T in zip(xs, ys, Ts):
                if x <= 0 or x >= 1:
                    continue
                xv = [1 - x, x] if flip else [x, 1 - x]
                pt = bubble_T(L, V, P, xv)
                ym = pt.y[1] if flip else pt.y[0]
                dT.append(pt.T - T); dy.append(ym - y)
                points.append({"system": f"{a}+{b}", "model": model, "x1": x, "T_exp_K": T, "T_model_K": pt.T,
                               "y1_exp": y, "y1_model": ym})
            d = vle_diagram(L, V, P=P, npts=81)
            az = d.azeotrope
            az_m = None
            if az is not None:
                az_m = (1 - az[0], az[1]) if flip else (az[0], az[1])
            summary.append({"system": f"{a}+{b}", "model": model, "P_kPa": P / 1e3, "n": len(dT),
                            "AAD_T_K": float(np.mean(np.abs(dT))), "max_dT_K": float(np.max(np.abs(dT))),
                            "AAD_y1": float(np.mean(np.abs(dy))), "max_dy1": float(np.max(np.abs(dy))),
                            "x_az_exp": az_exp[0] if az_exp else "none", "T_az_exp_K": az_exp[1] if az_exp else "none",
                            "x_az_model": az_m[0] if az_m else "none", "T_az_model_K": az_m[1] if az_m else "none",
                            "source": fname})
    write("binary_activity_points.csv", points)
    write("binary_activity_summary.csv", summary)
    return summary


# --------------------------------------------------------------------------
# C. methane + ethane with Peng-Robinson against measurements
# --------------------------------------------------------------------------
def binary_eos():
    data = read("vle_methane_ethane_Han2012.csv")
    variants = [("PR", "E-PPR78 k12(T)", {}), ("PR", "DECHEMA k12 = -0.0059", {"kij": -0.0059}),
                ("PR", "k12 = 0", {"kij": 0.0}), ("SRK", "k12 = 0", {"kij": 0.0}),
                ("vdW", "k12 = 0", {"kij": 0.0})]
    points, summary = [], []
    for model, label, par in variants:
        dP, dy = [], []
        cache = {}
        for r in data:
            T, Pe, x, y = float(r["T_K"]), float(r["P_kPa"]) * 1e3, float(r["x1"]), float(r["y1"])
            if T not in cache:
                cache[T] = _phases({"components": ["methane", "ethane"], "model": model, "hold": "T",
                                    "T": T, "params": par}, "methane")
            L, V, flip = cache[T]
            xv = [1 - x, x] if flip else [x, 1 - x]
            try:
                pt = bubble_P(L, V, T, xv)
            except Exception:
                continue
            ym = pt.y[1] if flip else pt.y[0]
            dP.append(100 * (pt.P / Pe - 1)); dy.append(ym - y)
            points.append({"model": model, "kij": label, "T_K": T, "x_CH4": x, "P_exp_kPa": Pe / 1e3,
                           "P_model_kPa": pt.P / 1e3, "y_exp": y, "y_model": ym})
        summary.append({"model": model, "kij": label, "n": len(dP), "AAD_P_pct": float(np.mean(np.abs(dP))),
                        "bias_P_pct": float(np.mean(dP)), "AAD_y": float(np.mean(np.abs(dy)))})
    write("methane_ethane_points.csv", points)
    write("methane_ethane_summary.csv", summary)
    return summary


# --------------------------------------------------------------------------
# D. the field calculations against the reference equations of state
# --------------------------------------------------------------------------
def field():
    ref = read("reference_field.csv")
    out = []
    by = lambda case: [r for r in ref if r["case"].startswith(case)]
    for model in ("SRK", "PR"):
        e, M, _ = pf.make_eos("propane", model)
        for r in by("band_propane_100m"):
            T, gf = float(r["T_K"]), float(r["factor"])
            fld = pf.PureField("gravity", G0 * gf, 0.0, 100.0, 0.0)
            try:
                _, Ps, P1, P2 = pf.band(e, M, fld, [T])
            except Exception:
                continue
            out.append({"case": "band width, propane, 100 m", "model": model, "T_K": T, "g_factor": gf,
                        "reference": float(r["value"]) / 1e5, "model_value": float(P2[0] - Ps[0]) / 1e5,
                        "unit": "bar", "dev_pct": 100 * (float(P2[0] - Ps[0]) / float(r["value"]) - 1)})
        # interface height with 10 bar at the bottom, 300 K
        r = by("interface_propane")[0]
        col = pf.column(e, M, 300.0, 10e5, pf.PureField("gravity", G0, 0.0, 20.0, 0.0))
        zi = col.interface if col.interface is not None else float("nan")   # nan: no liquid at all
        out.append({"case": "interface height, propane, 10 bar at bottom", "model": model, "T_K": 300.0,
                    "g_factor": 1, "reference": float(r["value"]), "model_value": zi, "unit": "m",
                    "dev_pct": 100 * (zi / float(r["value"]) - 1)})
        # the same, at the same excess pressure above each model's own Psat
        Ps_ref = float(r["Psat_Pa"])
        Ps_m = pf.saturation_pressure(e, 300.0)[0]
        col2 = pf.column(e, M, 300.0, Ps_m + (10e5 - Ps_ref), pf.PureField("gravity", G0, 0.0, 20.0, 0.0))
        out.append({"case": "interface height, propane, same P - Psat at bottom", "model": model, "T_K": 300.0,
                    "g_factor": 1, "reference": float(r["value"]), "model_value": col2.interface, "unit": "m",
                    "dev_pct": 100 * (col2.interface / float(r["value"]) - 1)})
        # CO2 centrifuge: I = int dP / rho_V from the inner pressure to Psat
        e2, M2, _ = pf.make_eos("carbon-dioxide", model)
        Ps2 = pf.saturation_pressure(e2, 280.0)[0]
        ln_fs = pf._ln_f(e2, 280.0, Ps2, "vapor")
        for r in by("centrifuge_CO2"):
            frac = float(r["factor"])
            I_m = R * 280.0 * (ln_fs - pf._ln_f(e2, 280.0, frac * Ps2, "vapor")) / M2
            w = 20000 * 2 * math.pi / 60
            r_ref = math.sqrt(0.02 ** 2 + 2 * float(r["value"]) / w ** 2)
            r_m = math.sqrt(0.02 ** 2 + 2 * I_m / w ** 2)
            out.append({"case": f"interface radius, CO2 centrifuge 20000 rpm, P_in = {frac} Psat", "model": model,
                        "T_K": 280.0, "g_factor": "", "reference": 100 * r_ref, "model_value": 100 * r_m,
                        "unit": "cm", "dev_pct": 100 * (r_m / r_ref - 1)})
    # water column: PR and SRK only (the activity route has no liquid EOS)
    r = by("water_top")[0]
    for model in ("SRK", "PR"):
        e, M, _ = pf.make_eos("water", model)
        col = pf.column(e, M, 300.0, 10e5, pf.PureField("gravity", G0, 0.0, 50.0, 0.0))
        out.append({"case": "water column 50 m, 10 bar at bottom: top pressure", "model": model, "T_K": 300.0,
                    "g_factor": 1, "reference": float(r["value"]) / 1e5, "model_value": float(col.P[-1]) / 1e5,
                    "unit": "bar", "dev_pct": 100 * ((10e5 - col.P[-1]) / (10e5 - float(r["value"])) - 1)})
    write("field_vs_reference.csv", out)
    return out


def main():
    A = pure_saturation()
    print("A. pure fluids, AAD % against reference EOS (Psat / rhoL / rhoV)")
    for s in A:
        print(f"  {s['compound']:15s} {s['model']:4s} {s['AAD_Psat_pct']:7.2f} {s['AAD_rhoL_pct']:7.2f} "
              f"{s['AAD_rhoV_pct']:7.2f}")
    B = binary_activity()
    print("B. binary VLE at 1 atm")
    for s in B:
        print(f"  {s['system']:20s} {s['model']:6s} AAD T {s['AAD_T_K']:.2f} K, AAD y {s['AAD_y1']:.4f}, "
              f"az exp {s['x_az_exp']}, model {s['x_az_model']}")
    C = binary_eos()
    print("C. methane + ethane")
    for s in C:
        print(f"  {s['model']:4s} {s['kij']:24s} AAD P {s['AAD_P_pct']:.2f} %, bias {s['bias_P_pct']:.2f} %, "
              f"AAD y {s['AAD_y']:.4f}")
    D = field()
    print("D. field")
    for s in D:
        print(f"  {s['case'][:55]:55s} {s['model']:4s} {s.get('T_K')} g{s.get('g_factor')}: "
              f"ref {s['reference']:.4g}, model {s['model_value']:.4g} {s['unit']} ({s['dev_pct']:+.1f} %)")


if __name__ == "__main__":
    main()
