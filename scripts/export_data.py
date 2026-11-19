"""Write the input and output data of every scenario the project handles, as
plain text and CSV files, for the submission's Data folder.

    python scripts/export_data.py              # writes submission/Data
    python scripts/export_data.py --out DIR

Every scenario gets its own numbered folder holding

    input.txt            the input file, exactly as the program reads it
    checker_report.txt   what the input checker concluded (sufficient, missing data, ...)
    results.txt          the solver's printed results (for the example input files)
    *.csv                the computed numbers, one row per point, units in the column names

Part A runs every file in examples/ through the command-line program (input
checker, then solver). Parts B to E reproduce what the website calculates:
binary diagrams, P-T envelopes, external fields, mixtures of three or four
components, and the validation tables. A README.txt index is written at the
top of the output folder.
"""
import argparse
import contextlib
import csv
import io
import math
import os
import shutil
import sys
import warnings
from dataclasses import is_dataclass, fields

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
warnings.filterwarnings("ignore")
np.seterr(all="ignore")

import compound_db                                    # noqa: E402
import equilibrium_solver as es                       # noqa: E402
import web_api                                        # noqa: E402
import pure_fluid                                     # noqa: E402
from input_checker import check_input                 # noqa: E402
from thermo_models import build_system                # noqa: E402

INDEX = []                                            # (folder, description, files)


# --------------------------------------------------------------------------
# writing helpers
# --------------------------------------------------------------------------
def _num(v):
    if v is None:
        return ""
    if isinstance(v, (bool, np.bool_)):
        return str(bool(v))
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    if isinstance(v, (float, np.floating)):
        return "" if not math.isfinite(float(v)) else f"{float(v):.8g}"
    return str(v)


def write_csv(path, columns):
    """columns: dict name -> list/array (all the same length) or scalar."""
    n = max((len(v) for v in columns.values() if hasattr(v, "__len__") and not isinstance(v, str)),
            default=1)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(columns.keys())
        for i in range(n):
            row = []
            for v in columns.values():
                if hasattr(v, "__len__") and not isinstance(v, str):
                    row.append(_num(v[i]) if i < len(v) else "")
                else:
                    row.append(_num(v) if i == 0 else "")
            w.writerow(row)
    return os.path.basename(path)


def write_kv(path, rows):
    """Two-column CSV of quantity, value."""
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["quantity", "value"])
        for k, v in rows:
            w.writerow([k, _num(v)])
    return os.path.basename(path)


def write_text(path, text):
    with open(path, "w") as fh:
        fh.write(text if text.endswith("\n") else text + "\n")
    return os.path.basename(path)


def slug(s):
    return "".join(c if c.isalnum() else "_" for c in s).strip("_").lower()


class Scenario:
    counter = 0

    def __init__(self, out, part, name, description):
        Scenario.counter += 1
        self.folder = f"{Scenario.counter:02d}_{part}_{slug(name)}"
        self.path = os.path.join(out, self.folder)
        os.makedirs(self.path, exist_ok=True)
        self.description = description
        self.files = []

    def file(self, fname):
        return os.path.join(self.path, fname)

    def add(self, fname):
        self.files.append(fname)

    def close(self):
        INDEX.append((self.folder, self.description, self.files))
        print(f"  {self.folder}: {len(self.files)} files")


# --------------------------------------------------------------------------
# A. the example input files through the command-line program
# --------------------------------------------------------------------------
UNITS = {"P": "P_bar", "T": "T_K", "f1": "f1_bar", "f2": "f2_bar", "P_dew": "P_dew_bar",
         "P_bubble": "P_bubble_bar", "position": "position_m", "rho": "density_kg_m3",
         "P_sat": "P_sat_bar", "density": "density_kg_m3"}


def _result_tables(sc, res, k):
    """CSV files for one Result of the command-line solver."""
    d = res.data
    tag = f"{k}_{slug(res.equilibrium)}"
    if isinstance(d, es.Diagram):
        L, V = res.phases
        tab = es.diagram_properties(L, V, d)
        cols = {}
        for key, v in tab.items():
            if key in ("P", "f1", "f2"):
                v = np.asarray(v) / 1e5
            cols[UNITS.get(key, key)] = v
        sc.add(write_csv(sc.file(f"{tag}_{d.kind}.csv"), cols))
        summ = [("diagram", d.kind), ("held constant", f"{d.fixed / 1e5:.6g} bar" if d.kind == "Txy"
                                      else f"{d.fixed:.6g} K"),
                ("component 1", d.names[0]), ("component 2", d.names[1])]
        if d.azeotrope:
            summ += [("azeotrope x1", d.azeotrope[0]),
                     ("azeotrope " + ("P (bar)" if d.kind == "Pxy" else "T (K)"),
                      d.azeotrope[1] / 1e5 if d.kind == "Pxy" else d.azeotrope[1])]
        if d.three_phase:
            for kk, vv in d.three_phase.items():
                summ.append((f"three-phase {kk}", vv if not isinstance(vv, (list, tuple, np.ndarray))
                             else ", ".join(_num(x) for x in np.ravel(vv))))
        if d.unstable:
            summ.append(("unstable liquid between x1", f"{d.unstable[0]:.4g} and {d.unstable[1]:.4g}"))
        for note in d.notes:
            summ.append(("note", note))
        sc.add(write_kv(sc.file(f"{tag}_summary.csv"), summ))
    elif isinstance(d, es.Envelope):
        cols = {"T_K": d.T, "P_bubble_bar": d.P_bubble / 1e5, "P_dew_bar": d.P_dew / 1e5}
        for i, nm in enumerate(d.names):
            cols[f"y_{nm}_at_bubble"] = d.y_bubble[:, i]
            cols[f"x_{nm}_at_dew"] = d.x_dew[:, i]
        sc.add(write_csv(sc.file(f"{tag}_PT_envelope.csv"), cols))
    elif isinstance(d, es.SLECurve):
        sc.add(write_csv(sc.file(f"{tag}_solubility_curve.csv"),
                         {f"x_{d.solid}_liquid": d.x_solid, "T_K": d.T, "P_bar": d.P / 1e5}))
    elif isinstance(d, dict):
        rows = []
        for key, v in d.items():
            if isinstance(v, (list, tuple, np.ndarray)):
                v = ", ".join(_num(x) for x in np.ravel(v))
            elif key == "P" and isinstance(v, (int, float)):
                key, v = "P_bar", v / 1e5
            rows.append((key, v))
        sc.add(write_kv(sc.file(f"{tag}_values.csv"), rows))
    elif is_dataclass(d):
        rows = [(f.name, getattr(d, f.name)) for f in fields(d)]
        sc.add(write_kv(sc.file(f"{tag}_values.csv"), [(a, b if not isinstance(b, np.ndarray) else
                                                         ", ".join(_num(x) for x in np.ravel(b)))
                                                        for a, b in rows]))
    for j, prof in enumerate(res.profiles):
        cols = {"position_m": prof.position, "P_bar": prof.P / 1e5}
        for i, nm in enumerate(prof.names):
            cols[f"x_{nm}"] = prof.x[:, i]
        cols["phase"] = prof.phase
        cols["density_kg_m3"] = prof.density
        sc.add(write_csv(sc.file(f"{tag}_{prof.kind}_profile_{j + 1}.csv"), cols))


def part_a(out):
    print("A. example input files")
    ex_dir = os.path.join(ROOT, "examples")
    for fname in sorted(os.listdir(ex_dir)):
        if not fname.endswith(".txt"):
            continue
        text = open(os.path.join(ex_dir, fname)).read()
        first = next((ln.lstrip("# ").strip() for ln in text.splitlines() if ln.startswith("#")), "")
        rep = check_input(text)
        sc = Scenario(out, "input", fname[:-4],
                      f"examples/{fname} through the command-line program. Checker verdict: "
                      f"{rep.status}. {first}")
        sc.add(write_text(sc.file("input.txt"), text))
        sc.add(write_text(sc.file("checker_report.txt"), str(rep)))
        if rep.status == "SUFFICIENT" and len(rep.inputs.components) == 1:      # a pure fluid
            txt, charts = pure_fluid.solve_input(rep)
            sc.add(write_text(sc.file("results.txt"), txt))
            for k, ch in enumerate(charts, 1):
                cols = {}
                for se in ch["series"]:
                    cols[f"{se['name']}: {ch['xTitle']}"] = [q[0] for q in se["pts"]]
                    cols[f"{se['name']}: {ch['yTitle']}"] = [q[1] for q in se["pts"]]
                sc.add(write_csv(sc.file(f"{k}_{slug(ch['title'])}.csv"), cols))
        elif rep.status == "SUFFICIENT":
            system = build_system(text)
            results = es.solve(system)
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                es.print_results(system, results)
            sc.add(write_text(sc.file("results.txt"), buf.getvalue()))
            for k, res in enumerate(results, 1):
                _result_tables(sc, res, k)
        sc.close()


# --------------------------------------------------------------------------
# B-E. what the website calculates
# --------------------------------------------------------------------------
def _web(sc, r, stem):
    """input.txt, checker_report.txt and a CSV of columns for one compute()."""
    if not r.get("ok"):
        sc.add(write_text(sc.file(f"{stem}_error.txt"), r.get("error", "failed")))
        return
    if "input_text" in r:
        sc.add(write_text(sc.file(f"{stem}_input.txt"), r["input_text"]))
    if "report" in r:
        sc.add(write_text(sc.file(f"{stem}_checker_report.txt"), r["report"]))
    names = r.get("names", ["1", "2"])
    cols = {}
    for k, v in r["columns"].items():
        col = UNITS.get(k, k)
        col = col.replace("x1", f"x_{names[0]}").replace("y1", f"y_{names[0]}")
        cols[col] = v
    sc.add(write_csv(sc.file(f"{stem}.csv"), cols))
    summ = [("components", " + ".join(names)), ("route", r.get("route")),
            ("calculation", r.get("calculation"))]
    if r.get("azeotrope"):
        for k, v in r["azeotrope"].items():
            summ.append((f"azeotrope {UNITS.get(k, k).replace('x1', 'x_' + names[0])}", v))
    if r.get("three_phase"):
        for k, v in r["three_phase"].items():
            summ.append((f"three-phase {k}", v))
    if r.get("unstable"):
        summ.append(("unstable liquid between x1", ", ".join(_num(v) for v in r["unstable"])))
    for c in r.get("contacts") or []:
        summ.append(("phase boundary in the field", ", ".join(f"{k} = {_num(v)}" for k, v in c.items())))
    for n_ in r.get("notes", []):
        summ.append(("note", n_))
    sc.add(write_kv(sc.file(f"{stem}_summary.csv"), summ))


DIAGRAMS = [
    # (name, request, description)
    ("ethanol_water_nrtl_Txy_1atm",
     {"components": ["ethanol", "water"], "model": "nrtl", "hold": "P", "P": 1.01325},
     "Ethanol + water, NRTL (DECHEMA parameters), Txy at 1 atm: minimum-boiling azeotrope."),
    ("ethanol_water_ideal_Txy_1atm",
     {"components": ["ethanol", "water"], "model": "ideal", "hold": "P", "P": 1.01325},
     "Ethanol + water as an ideal solution (Raoult's law), Txy at 1 atm, for comparison with NRTL."),
    ("ethanol_water_nrtl_Pxy_350K",
     {"components": ["ethanol", "water"], "model": "nrtl", "hold": "T", "T": 350},
     "Ethanol + water, NRTL, Pxy at 350 K, with activity coefficients and fugacities."),
    ("acetone_chloroform_nrtl_Txy_1atm",
     {"components": ["acetone", "chloroform"], "model": "nrtl", "hold": "P", "P": 1.01325},
     "Acetone + chloroform, NRTL, Txy at 1 atm: maximum-boiling azeotrope (negative deviation)."),
    ("benzene_toluene_ideal_Pxy_360K",
     {"components": ["benzene", "toluene"], "model": "ideal", "hold": "T", "T": 360},
     "Benzene + toluene, ideal solution, Pxy at 360 K: a nearly ideal pair."),
    ("benzene_toluene_nrtl_Pxy_360K",
     {"components": ["benzene", "toluene"], "model": "nrtl", "hold": "T", "T": 360},
     "Benzene + toluene, NRTL, Pxy at 360 K (compare with the ideal solution)."),
    ("methanol_water_nrtl_Txy_1atm",
     {"components": ["methanol", "water"], "model": "nrtl", "hold": "P", "P": 1.01325},
     "Methanol + water, NRTL, Txy at 1 atm: no azeotrope."),
    ("propane_nbutane_PR_Pxy_300K",
     {"components": ["propane", "n-butane"], "model": "PR", "hold": "T", "T": 300},
     "Propane + n-butane, Peng-Robinson with E-PPR78 k12(T), Pxy at 300 K."),
    ("propane_nbutane_vdW_Pxy_300K",
     {"components": ["propane", "n-butane"], "model": "vdW", "hold": "T", "T": 300},
     "Propane + n-butane, van der Waals (k12 = 0), Pxy at 300 K, for comparison with Peng-Robinson."),
    ("co2_ethane_PR_PPR78_Pxy_263K",
     {"components": ["carbon-dioxide", "ethane"], "model": "PR", "hold": "T", "T": 263.15},
     "CO2 + ethane, Peng-Robinson with E-PPR78 k12(T), Pxy at 263.15 K: the azeotrope from the course slides."),
    ("co2_ethane_PR_k0_Pxy_263K",
     {"components": ["carbon-dioxide", "ethane"], "model": "PR", "hold": "T", "T": 263.15,
      "params": {"kij": 0.0}},
     "CO2 + ethane, Peng-Robinson with k12 = 0, Pxy at 263.15 K: no azeotrope."),
    ("methane_nbutane_PR_Pxy_250K",
     {"components": ["methane", "n-butane"], "model": "PR", "hold": "T", "T": 250},
     "Methane + n-butane, Peng-Robinson with E-PPR78, Pxy at 250 K: methane is supercritical, "
     "the curves close at a mixture critical point."),
    ("nhexane_noctane_PR_Txy_1atm",
     {"components": ["n-hexane", "n-octane"], "model": "PR", "hold": "P", "P": 1.01325},
     "n-Hexane + n-octane, Peng-Robinson with E-PPR78, Txy at 1 atm."),
]

ENVELOPES = [
    ("propane_nbutane_PR_PT_envelope_x05",
     {"components": ["propane", "n-butane"], "model": "PR", "hold": "x", "x1": 0.5},
     "Propane + n-butane (50/50), Peng-Robinson: bubble and dew pressure against temperature."),
    ("ethanol_water_nrtl_PT_envelope_x05",
     {"components": ["ethanol", "water"], "model": "nrtl", "hold": "x", "x1": 0.5},
     "Ethanol + water (50/50), NRTL: bubble and dew pressure against temperature."),
]

FIELDS = [
    ("ethanol_water_gravity_500m",
     {"components": ["ethanol", "water"], "model": "nrtl", "hold": "field", "T": 350, "P": 1.01325,
      "z1": 0.4, "field": {"kind": "gravity", "g": 9.81, "height": 500}},
     "Ethanol + water (NRTL, liquid volumes from the Rackett equation) in a 500 m column under "
     "gravity at 350 K: pressure, composition and phase along the height."),
    ("propane_nbutane_gravity_500m",
     {"components": ["propane", "n-butane"], "model": "PR", "hold": "field", "T": 320, "P": 25.0,
      "z1": 0.4, "field": {"kind": "gravity", "g": 9.80665, "height": 500}},
     "Propane + n-butane (Peng-Robinson, E-PPR78) in a 500 m column under gravity at 320 K, 25 bar at "
     "the bottom: liquid below, a gas cap above."),
    ("ethanol_water_centrifuge",
     {"components": ["ethanol", "water"], "model": "nrtl", "hold": "field", "T": 350, "P": 0.9,
      "z1": 0.4, "field": {"kind": "centrifugal", "rpm": 20000, "r_in": 0.05, "r_out": 0.5}},
     "Ethanol + water (NRTL) in a centrifuge at 20 000 rpm between radii 0.05 and 0.5 m, 350 K."),
]


def part_b(out):
    print("B. binary phase diagrams (as on the website)")
    for name, req, desc in DIAGRAMS:
        sc = Scenario(out, "diagram", name, desc)
        _web(sc, web_api.compute(dict(req, npts=req.get("npts", 61))), "diagram")
        sc.close()
    print("C. P-T envelopes")
    for name, req, desc in ENVELOPES:
        sc = Scenario(out, "envelope", name, desc)
        _web(sc, web_api.compute(dict(req, npts=41)), "envelope")
        sc.close()
    print("D. external fields")
    for name, req, desc in FIELDS:
        sc = Scenario(out, "field", name, desc)
        _web(sc, web_api.compute(req), "profile")
        sc.close()


PURE = ["water", "carbon-dioxide", "nitrogen", "methane", "propane", "n-butane", "n-hexane",
        "benzene", "methanol", "ethanol"]


def part_pure(out):
    print("P. pure fluids")
    # vapor-pressure curves and coexisting densities, every compound and EOS
    for cid in PURE:
        rec = compound_db.find(cid)
        sc = Scenario(out, "pure", f"{cid}_saturation",
                      f"Pure {rec['name']}: vapor pressure, coexisting densities and Clapeyron enthalpy "
                      f"of vaporization from van der Waals, SRK and Peng-Robinson, with the Antoine "
                      f"vapor pressure (Poling et al.) where it applies.")
        rows = []
        for m in ("vdW", "SRK", "PR"):
            r = web_api.pure({"compound": cid, "model": m, "T": 0.7 * rec["Tc"]["value"], "P": 1.01325})
            c = r["curve"]
            sc.add(write_csv(sc.file(f"saturation_{m}.csv"),
                             {"T_K": c["T"], "Psat_bar": c["P"], "rho_liquid_kg_m3": c["rhoL"],
                              "rho_vapor_kg_m3": c["rhoV"], "dHvap_kJ_mol": c["dHvap"],
                              "Psat_Antoine_bar": c["P_antoine"]}))
            rows += [(f"{m}: critical point Tc (K)", r["Tc"]), (f"{m}: Pc (bar)", r["Pc"]),
                     (f"{m}: critical density (kg/m3)", r["rhoc"]),
                     (f"{m}: normal boiling point (K)", r["Tb_eos"])]
        rows.append(("normal boiling point from the Antoine data (K)", r["Tb_data"]))
        sc.add(write_kv(sc.file("summary.csv"), rows))
        sc.close()
    # isotherms with the equal-area line
    for cid, m in (("carbon-dioxide", "vdW"), ("carbon-dioxide", "PR")):
        sc = Scenario(out, "pure", f"{cid}_{m}_isotherms",
                      f"Pure CO2, {m}: isotherms at T/Tc = 0.8 ... 1.1 in reduced units (P/Pc against "
                      f"V/Vc), the saturated volumes and the equal-area (Maxwell) pressure of each.")
        r = web_api.pure({"compound": cid, "model": m, "T": 280, "P": 40})
        for it in r["isotherms"]:
            sc.add(write_csv(sc.file(f"isotherm_Tr_{it['Tr']:.2f}.csv"), {"V_over_Vc": it["V"], "P_over_Pc": it["P"]}))
        sc.add(write_csv(sc.file("maxwell_lines.csv"),
                         {"T_over_Tc": [it["Tr"] for it in r["isotherms"]],
                          "Psat_over_Pc": [it["Psat"] for it in r["isotherms"]],
                          "vL_over_Vc": [it["vL"] for it in r["isotherms"]],
                          "vV_over_Vc": [it["vV"] for it in r["isotherms"]],
                          "area_mismatch_fraction": [it["area"] for it in r["isotherms"]]}))
        sc.add(write_csv(sc.file("two_phase_dome.csv"), {"vL_over_Vc": r["dome"]["vL"],
                                                           "vV_over_Vc": r["dome"]["vV"],
                                                           "Psat_over_Pc": r["dome"]["P"]}))
        sc.close()
    # fields
    cases = [("propane_gravity_20m", "propane", 300.0, 10.0, {"kind": "gravity", "g": 9.81, "height": 20.0},
              "Pure propane at 300 K in a 20 m column (g = 9.81 m/s2), 10.0 bar at the bottom: liquid at the "
              "bottom, the interface, vapor above."),
             ("propane_gravity_100m_10g", "propane", 300.0, 10.0, {"kind": "gravity", "g": 98.1, "height": 100.0},
              "Pure propane at 300 K in a 100 m column at g = 98.1 m/s2 (10 g)."),
             ("water_gravity_100m", "water", 373.15, 1.5, {"kind": "gravity", "g": 9.81, "height": 100.0},
              "Pure water at 373.15 K in a 100 m column, 1.5 bar at the bottom."),
             ("co2_centrifuge", "carbon-dioxide", 280.0, 35.0,
              {"kind": "centrifugal", "rpm": 20000, "r_in": 0.02, "r_out": 0.15},
              "Pure CO2 at 280 K in a centrifuge at 20 000 rpm between radii 0.02 and 0.15 m, 35 bar at the "
              "inner radius: vapor inside, liquid at the rim."),
             ("nitrogen_gravity_1000m", "nitrogen", 300.0, 1.0, {"kind": "gravity", "g": 9.81, "height": 1000.0},
              "Pure nitrogen at 300 K (above Tc) in a 1000 m column: one phase, close to the barometric formula.")]
    for name, cid, T, P, f, desc in cases:
        sc = Scenario(out, "pure_field", name, desc)
        r = web_api.pure({"compound": cid, "model": "PR", "T": T, "P": P, "field": f})
        F = r["field"]
        cl = F["column"]
        sc.add(write_csv(sc.file("profile.csv"), {"position_m": cl["position"], "P_bar": cl["P"],
                                                   "density_kg_m3": cl["rho"], "phase": cl["phase"]}))
        rows = [("T_K", T), ("P at the reference position (bar)", P), ("reference position (m)", F["reference"]),
                ("vapor pressure at T (bar)", cl["Psat"]), ("interface position (m)", cl["interface"])]
        rows += [("note", n_) for n_ in cl["notes"]]
        sc.add(write_kv(sc.file("summary.csv"), rows))
        for b in F["bands"]:
            sc.add(write_csv(sc.file(f"phase_boundary_band_{b['fraction']:.2f}_strength.csv"),
                             {"T_K": b["T"], "Psat_bar": b["Psat"], "P_first_liquid_bar": b["P_first"],
                              "P_last_vapor_bar": b["P_last"]}))
        if "interface_curve" in F:
            ic = F["interface_curve"]
            sc.add(write_csv(sc.file("interface_vs_field_strength.csv"),
                             {"g_m_s2" if f["kind"] == "gravity" else "speed_rpm": ic["strength"],
                              "interface_position_m": ic["position"]}))
        sc.close()


def part_mix(out):
    print("E. mixtures of three and four components")
    for ex in web_api.MIXTURE_EXAMPLES:
        n = len(ex["components"])
        sc = Scenario(out, f"{n}comp", ex["id"] + "_" + ex["model"],
                      f"{ex['label']}: flash at {ex['T']:g} K and {ex['P']:g} bar"
                      + (", and the bubble points over the triangle of liquid compositions"
                         if n == 3 else "") + ".")
        r = web_api.mixture(dict(ex, map=n == 3))
        if not r["ok"]:
            sc.add(write_text(sc.file("error.txt"), r["error"]))
            sc.close()
            continue
        sc.add(write_text(sc.file("input.txt"), r["input_text"]))
        sc.add(write_text(sc.file("checker_report.txt"), r["report"]))
        cols = {"component": r["names"], "z_feed": r["z"], "x_liquid": r["x"], "y_vapor": r["y"]}
        if r.get("K"):
            cols["K"] = r["K"]
        if r.get("gamma"):
            cols["gamma"] = r["gamma"]
        if r.get("phiL"):
            cols["phi_liquid"], cols["phi_vapor"] = r["phiL"], r["phiV"]
        sc.add(write_csv(sc.file("flash_compositions.csv"), cols))
        rows = [("T_K", r["T"]), ("P_bar", r["P"]), ("phase", r["phase"]),
                ("vapor fraction beta", r["beta"]), ("route", r["route"])]
        rows += [(k + ("_bar" if k.startswith("P") else "_K"), v) for k, v in r["bounds"].items()]
        for k, v in (r.get("kij") or {}).items():
            rows.append((f"kij {k.replace('|', ' / ')} at T", v))
        sc.add(write_kv(sc.file("flash_summary.csv"), rows))
        if r.get("map"):
            m = r["map"]
            key = m["kind"]
            nm = r["names"]
            pts = sorted(m["points"], key=lambda p: (-p["x"][0], -p["x"][1]))
            cols = {f"x_{nm[i]}": [p["x"][i] for p in pts] for i in range(3)}
            cols["bubble_T_K" if key == "T" else "bubble_P_bar"] = [p[key] for p in pts]
            for i in range(3):
                cols[f"y_{nm[i]}"] = [p["y"][i] if p["y"] else None for p in pts]
            sc.add(write_csv(sc.file("triangle_bubble_points.csv"), cols))
        sc.close()


def part_validation(out):
    print("F. validation")
    # E-PPR78 against the course slides and DECHEMA
    sc = Scenario(out, "validation", "eppr78_kij",
                  "E-PPR78 k12(T) against the values in the course slides (CO2 + ethane, CO2 + "
                  "methane) and against the fixed DECHEMA k12 of each pair that has one.")
    rows = []
    for T in (230.0, 263.15, 273.15, 300.0, 350.0):
        rows.append({"T_K": T, "kij_CO2_ethane": compound_db.PPR78Kij("carbon-dioxide", "ethane")(T),
                     "kij_CO2_methane": compound_db.PPR78Kij("carbon-dioxide", "methane")(T)})
    sc.add(write_csv(sc.file("kij_vs_T_course_pairs.csv"), {k: [r_[k] for r_ in rows] for k in rows[0]}))
    db = compound_db.load()
    ids = [c["id"] for c in db["compounds"] if "ppr78" in c]
    a_, b_, k1, k2, kd = [], [], [], [], []
    for i, a in enumerate(ids):
        for b in ids[i + 1:]:
            if not compound_db.ppr78_ready(a, b):
                continue
            p = compound_db.pair(a, b)
            a_.append(a); b_.append(b)
            k1.append(compound_db.PPR78Kij(a, b)(298.15)); k2.append(compound_db.PPR78Kij(a, b)(250.0))
            kd.append(p.get("kij_pr"))
    sc.add(write_csv(sc.file("kij_all_pairs.csv"), {"compound_1": a_, "compound_2": b_,
                                                     "kij_EPPR78_298K": k1, "kij_EPPR78_250K": k2,
                                                     "kij_DECHEMA_fixed": kd}))
    sc.close()

    # azeotropes of every NRTL pair at 1 atm
    sc = Scenario(out, "validation", "nrtl_azeotropes_1atm",
                  "Every pair of the activity-model compounds with NRTL parameters: Txy at 1 atm, "
                  "whether there is an azeotrope, where, and whether the liquid splits.")
    from itertools import combinations
    import json
    ed = json.load(open(os.path.join(ROOT, "site", "edition.json")))
    act = ed["editions"][ed["active"]]["compound_sets"]["activity"]
    cols = {k: [] for k in ("compound_1", "compound_2", "azeotrope", "x1_azeotrope",
                             "T_azeotrope_K", "liquid_split", "Tb1_K", "Tb2_K")}
    for a, b in combinations(act, 2):
        if "dg_ab" not in compound_db.pair(a, b):
            continue
        r = web_api.compute({"components": [a, b], "model": "nrtl", "hold": "P", "P": 1.01325,
                             "npts": 61})
        az = r.get("azeotrope") if r.get("ok") else None
        cols["compound_1"].append(a); cols["compound_2"].append(b)
        cols["azeotrope"].append("yes" if az else ("failed" if not r.get("ok") else "no"))
        cols["x1_azeotrope"].append(az["x1"] if az else None)
        cols["T_azeotrope_K"].append(az["T"] if az else None)
        cols["liquid_split"].append("yes" if (r.get("three_phase") or r.get("unstable")) else "no")
        cols["Tb1_K"].append(compound_db.find(a)["Tb"]["value"])
        cols["Tb2_K"].append(compound_db.find(b)["Tb"]["value"])
    sc.add(write_csv(sc.file("azeotropes.csv"), cols))
    sc.close()

    # pure fluids: each equation of state against the Antoine vapor pressure
    sc = Scenario(out, "validation", "pure_vapor_pressure",
                  "Normal boiling points and vapor pressures at T = 0.7 Tc from van der Waals, SRK and "
                  "Peng-Robinson against the Antoine correlation (a fit to measured vapor pressures, "
                  "Poling et al.), for ten compounds. The comparison with reference equations of state "
                  "is in the validation_real_data scenarios.")
    cols = {k: [] for k in ("compound", "Tb_Antoine_K", "Tb_vdW_K", "Tb_SRK_K", "Tb_PR_K", "T_07Tc_K",
                             "Psat_Antoine_bar", "Psat_vdW_bar", "Psat_SRK_bar", "Psat_PR_bar")}
    for cid in PURE:
        rec = compound_db.find(cid)
        a_ = rec["antoine"]
        T7 = 0.7 * rec["Tc"]["value"]
        cols["compound"].append(cid)
        cols["Tb_Antoine_K"].append(rec["Tb"]["value"] if rec["Tb"] else None)
        cols["T_07Tc_K"].append(T7)
        ok7 = a_ and a_["Tmin"] - 5 <= T7 <= a_["Tmax"] + 5
        cols["Psat_Antoine_bar"].append(10 ** (a_["A"] - a_["B"] / (T7 + a_["C"])) / 1e5 if ok7 else None)
        for m in ("vdW", "SRK", "PR"):
            e_, M_, _ = pure_fluid.make_eos(cid, m)
            try:
                cols[f"Tb_{m}_K"].append(pure_fluid.saturation_temperature(e_, 101325.0))
            except ValueError:
                cols[f"Tb_{m}_K"].append(None)
            cols[f"Psat_{m}_bar"].append(pure_fluid.saturation_pressure(e_, T7)[0] / 1e5)
    sc.add(write_csv(sc.file("vapor_pressure_check.csv"), cols))
    sc.close()

    # flash checks: material balance and K = y/x
    sc = Scenario(out, "validation", "flash_checks",
                  "For each 3- and 4-component flash: the material balance (1 - beta) x + beta y - z "
                  "and the bubble/dew bracket, which should be ~0 and T_bubble < T < T_dew.")
    cols = {k: [] for k in ("mixture", "component", "z", "x", "y", "balance_residual")}
    for ex in web_api.MIXTURE_EXAMPLES:
        r = web_api.mixture(ex)
        if not r["ok"]:
            continue
        for i, nm in enumerate(r["names"]):
            cols["mixture"].append(ex["id"]); cols["component"].append(nm)
            cols["z"].append(r["z"][i]); cols["x"].append(r["x"][i]); cols["y"].append(r["y"][i])
            cols["balance_residual"].append((1 - r["beta"]) * r["x"][i] + r["beta"] * r["y"][i] - r["z"][i])
    sc.add(write_csv(sc.file("material_balance.csv"), cols))
    sc.close()


def part_real_data(out):
    """G. the models against measured VLE and reference equations of state
    (scripts/validate_real_data.py); the data files are copied next to the results."""
    print("G. validation against measured and reference data")
    import importlib.util
    spec = importlib.util.spec_from_file_location("vrd", os.path.join(ROOT, "scripts", "validate_real_data.py"))
    vrd = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(vrd)
    val = os.path.join(ROOT, "validation")
    jobs = [("pure saturation vs reference EOS", vrd.pure_saturation,
             ["reference_saturation.csv", "reference_eos_sources.json"],
             "Vapor pressure and saturated liquid and vapor densities of ten fluids from vdW, SRK and PR "
             "against their reference equations of state (IAPWS-95, Span-Wagner, Lemmon et al., ...), "
             "ten temperatures from 0.5 Tc to 0.95 Tc each; points and average absolute deviations."),
            ("binary VLE NRTL vs measured", vrd.binary_activity,
             ["vle_ethanol_water_101kPa_Kamihama2012.csv", "vle_acetone_ethanol_101kPa_KuTu2005.csv",
              "vle_chloroform_methanol_101kPa_Li2014.csv"],
             "NRTL and the ideal solution against measured isobaric VLE at 101.3 kPa (ethanol + water, "
             "acetone + ethanol, chloroform + methanol; data from the NIST TRC ThermoML archive)."),
            ("methane ethane EOS vs measured", vrd.binary_eos, ["vle_methane_ethane_Han2012.csv"],
             "PR (E-PPR78, DECHEMA and zero k12), SRK and vdW bubble pressures against 59 measured "
             "methane + ethane points at 126-140 K (Han et al. 2012)."),
            ("field vs reference EOS", vrd.field, ["reference_field.csv"],
             "Field results (propane band width, interface height, water column, CO2 centrifuge) from PR "
             "and SRK against the reference equations of state integrated through the same field.")]
    for name, fn, inputs, desc in jobs:
        sc = Scenario(out, "validation_real_data", name, desc)
        for f in inputs:
            shutil.copy(os.path.join(val, f), sc.file(f))
            sc.add(f)
        vrd.OUT = sc.path
        before = set(os.listdir(sc.path))
        fn()
        for f in sorted(set(os.listdir(sc.path)) - before):
            sc.add(f)
        sc.close()


def write_index(out):
    lines = ["DATA FILES: phase-equilibrium calculator",
             "=" * 40, "",
             "Produced by scripts/export_data.py from the code in Codes/. Rerun it to regenerate",
             "every file here. Each folder is one scenario. Common files:",
             "  input.txt / *_input.txt         the input file the program reads",
             "  checker_report.txt              the input checker's verdict and notes",
             "  results.txt                     the solver's printed results",
             "  *.csv                           computed values; units are in the column names",
             "                                  (T_K, P_bar, position_m, density_kg_m3, f1_bar ...)",
             "",
             "Column names in the diagram CSV files:",
             "  x_<compound>, y_<compound>  liquid and vapor mole fraction of component 1",
             "  K1, K2                      K-values y_i/x_i",
             "  alpha                       relative volatility K1/K2",
             "  gamma1, gamma2              activity coefficients (activity route)",
             "  phiL1, phiL2, phiV1, phiV2  fugacity coefficients in the liquid and vapor",
             "  f1_bar, f2_bar              fugacities",
             "  GE_RT, Gmix_RT              excess and mixing Gibbs energy over RT of the liquid",
             "  stability                   x1 x2 d2(Gmix/RT)/dx1^2; negative = the liquid splits",
             "  Empty cells: no solution at that point (for example above a critical point).",
             "", "SCENARIOS", "---------"]
    for folder, desc, files in INDEX:
        lines.append(f"{folder}/")
        lines.append("    " + desc)
        lines.append("    files: " + ", ".join(files))
        lines.append("")
    write_text(os.path.join(out, "README.txt"), "\n".join(lines))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", default=os.path.join(ROOT, "submission", "Data"))
    a = ap.parse_args(argv)
    if os.path.isdir(a.out):
        shutil.rmtree(a.out)
    os.makedirs(a.out)
    part_a(a.out)
    part_pure(a.out)
    part_b(a.out)
    part_mix(a.out)
    part_validation(a.out)
    part_real_data(a.out)
    write_index(a.out)
    print(f"wrote {len(INDEX)} scenarios to {a.out}")


if __name__ == "__main__":
    main()
