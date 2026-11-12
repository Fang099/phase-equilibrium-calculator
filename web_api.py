"""
web_api.py -- the interface between the website's calculator and the solver.

The calculator page runs these same Python modules in the browser (through
Pyodide). It sends a small JSON request describing what the user chose and
gets back the numbers to draw. Nothing is calculated in JavaScript.

Every request goes through the normal pipeline:

    request -> input text -> input_checker -> thermo_models -> equilibrium_solver

so the page can show the exact input file and the checker's report for what
is on screen. The same functions work from ordinary Python:

    import web_api
    out = web_api.compute({"components": ["ethanol", "water"], "model": "nrtl",
                           "hold": "P", "P": 1.01325})
    out["columns"]["x1"], out["columns"]["T"], out["azeotrope"]

Units in requests and responses: T in K, P and fugacities in bar,
compositions as mole fractions of the first component.
"""

import json
import math
import time

import numpy as np

import compound_db
from input_checker import check_input
from thermo_models import NRTL, R, build_system
import equilibrium_solver as es

EOS_MODELS = ("PR", "SRK", "vdW")
ACTIVITY_MODELS = ("ideal", "margules", "vanlaar", "nrtl")
MODEL_LABELS = {"PR": "Peng-Robinson", "SRK": "Soave-Redlich-Kwong", "vdW": "van der Waals",
                "ideal": "Ideal solution (Raoult's law)", "margules": "Margules",
                "vanlaar": "van Laar", "nrtl": "NRTL"}


# --------------------------------------------------------------------------
# Catalog: what the page can offer
# --------------------------------------------------------------------------
def catalog():
    """Compounds, which pairs have database parameters, models and sources."""
    db = compound_db.load()
    comps = []
    for c in db["compounds"]:
        a = c["antoine"]
        comps.append({
            "id": c["id"], "name": c["name"], "formula": c["formula"], "cas": c["cas"],
            "family": c["family"], "M": c["M"]["value"], "Tc": c["Tc"]["value"],
            "Pc_bar": c["Pc"]["value"] / 1e5, "omega": c["omega"]["value"],
            "Tb": c["Tb"]["value"],
            "antoine": {"A": a["A"], "B": a["B"], "C": a["C"], "Tmin": a["Tmin"],
                        "Tmax": a["Tmax"], "form": a["form"]},
            "sources": {"M": c["M"]["source"], "Tc": c["Tc"]["source"],
                        "Pc": c["Pc"]["source"], "omega": c["omega"]["source"],
                        "antoine": a["source"], "Tb": c["Tb"]["source"]},
        })
    pairs = [{"a": p["a"], "b": p["b"], "nrtl": "nrtl" in p, "kij_pr": "kij_pr" in p}
             for p in db["pairs"]]
    return {"compounds": comps, "pairs": pairs, "sources": db["meta"]["sources"],
            "generated": db["meta"]["generated"],
            "models": [{"id": m, "label": MODEL_LABELS[m],
                        "family": "eos" if m in EOS_MODELS else "activity"}
                       for m in EOS_MODELS + ACTIVITY_MODELS]}


def pair_defaults(a, b):
    """Starting parameter values and sensible ranges for a pair of compounds.

    NRTL parameters and the Peng-Robinson k_ij come from the database when it
    has them. Margules and van Laar have no database; their A12 and A21 are
    ln(gamma) at infinite dilution, so they are matched to the NRTL model's
    infinite-dilution values at the mean boiling point when NRTL parameters
    exist, and are zero (an ideal solution) otherwise."""
    ca, cb = compound_db.find(a), compound_db.find(b)
    if ca is None or cb is None:
        raise ValueError("unknown compound")
    pr = compound_db.pair(a, b)
    T_ref = 0.5 * (ca["Tb"]["value"] + cb["Tb"]["value"])
    out = {"a": ca["id"], "b": cb["id"], "T_ref": round(T_ref, 2),
           "kij": {"value": pr.get("kij_pr", 0.0),
                   "source": pr.get("kij_source", "none")},
           "nrtl": None, "A": {"A12": 0.0, "A21": 0.0, "source": "none"}}
    if "dg_ab" in pr:
        out["nrtl"] = {"dg12": round(pr["dg_ab"], 2), "dg21": round(pr["dg_ba"], 2),
                       "alpha": pr["alpha"], "source": pr["nrtl_source"]}
        m = NRTL(np.array([[0, pr["dg_ab"]], [pr["dg_ba"], 0]]),
                 np.array([[0, pr["alpha"]], [pr["alpha"], 0]]))
        A12 = float(m.ln_gamma(T_ref, [0.0, 1.0])[0])
        A21 = float(m.ln_gamma(T_ref, [1.0, 0.0])[1])
        out["A"] = {"A12": round(A12, 4), "A21": round(A21, 4), "source": "nrtl_infinite_dilution"}
    Tb_lo, Tb_hi = sorted((ca["Tb"]["value"], cb["Tb"]["value"]))
    Tc_lo = min(ca["Tc"]["value"], cb["Tc"]["value"])
    out["ranges"] = {
        "T": [int(5 * math.floor((Tb_lo - 60) / 5)), int(5 * math.ceil(min(Tb_hi + 80, Tc_lo - 5) / 5))],
        "T_default": int(5 * round(T_ref / 5)),
        "P": [0.05, 30.0], "P_default": 1.01325,
        "antoine_T": [max(ca["antoine"]["Tmin"], cb["antoine"]["Tmin"]),
                      min(ca["antoine"]["Tmax"], cb["antoine"]["Tmax"])],
    }
    return out


# --------------------------------------------------------------------------
# Request -> input text
# --------------------------------------------------------------------------
def build_input(req):
    """The input file that describes the request, exactly as the command-line
    tool would read it."""
    a, b = req["components"]
    model = req.get("model", "PR")
    par = req.get("params") or {}
    d = pair_defaults(a, b)
    a, b = d["a"], d["b"]
    lines = [f"components = {a}, {b}", "database = builtin"]
    if model in EOS_MODELS:
        lines.append(f"model = {model}")
        kij = par.get("kij", d["kij"]["value"] if model == "PR" else 0.0)
        lines.append(f"kij({a}, {b}) = {float(kij):g}")
    elif model in ACTIVITY_MODELS:
        lines.append(f"activity = {model}")
        if model in ("margules", "vanlaar"):
            lines.append(f"A({a}, {b}) = {float(par.get('A12', d['A']['A12'])):g}")
            lines.append(f"A({b}, {a}) = {float(par.get('A21', d['A']['A21'])):g}")
        elif model == "nrtl":
            n = d["nrtl"] or {"dg12": 0.0, "dg21": 0.0, "alpha": 0.3}
            lines.append(f"dg({a}, {b}) = {float(par.get('dg12', n['dg12'])):g} J/mol")
            lines.append(f"dg({b}, {a}) = {float(par.get('dg21', n['dg21'])):g} J/mol")
            lines.append(f"alpha({a}, {b}) = {float(par.get('alpha', n['alpha'])):g}")
    else:
        raise ValueError(f"unknown model '{model}'")
    hold = req.get("hold", "T")
    if hold == "T":
        lines += [f"T = {float(req['T']):g} K", "sweep = x"]
    elif hold == "P":
        lines += [f"P = {float(req['P']):g} bar", "sweep = x"]
    elif hold == "x":
        lines += [f"x(L) = {float(req['x1']):g}", "sweep = T"]
    else:
        raise ValueError("hold must be 'T', 'P' or 'x'")
    return "\n".join(lines) + "\n"


# --------------------------------------------------------------------------
# Compute
# --------------------------------------------------------------------------
def _clean(arr, digits=7):
    """Arrays -> JSON-safe lists (NaN and inf become null)."""
    out = []
    for v in np.asarray(arr, float).ravel():
        out.append(None if not np.isfinite(v) else float(f"{v:.{digits}g}"))
    return out


def _columns(table):
    """Solver table (SI) -> response columns (K, bar, mole fractions)."""
    cols = {}
    for k, v in table.items():
        v = np.asarray(v, float)
        if k in ("P", "P_dew", "f1", "f2"):
            v = v / 1e5
        cols[k] = _clean(v)
    return cols


def _solve(req):
    text = build_input(req)
    rep = check_input(text)
    out = {"input_text": text, "report": str(rep), "status": rep.status}
    if rep.status != "SUFFICIENT":
        out["ok"] = False
        out["error"] = "The checker did not accept this input: " + rep.status
        return out
    system = build_system(text)
    eq = next(e for e in system.equilibria if e.phase_names[1] == "V")
    L, V = eq.phases
    npts = int(req.get("npts", 41))
    hold = req.get("hold", "T")
    notes = []
    out.update(route=eq.route, calculation=eq.calculation, names=list(system.names),
               hold=hold, azeotrope=None)
    if hold in ("T", "P"):
        d = (es.vle_diagram(L, V, T=float(req["T"]), npts=npts, spacing="ends") if hold == "T"
             else es.vle_diagram(L, V, P=float(req["P"]) * 1e5, npts=npts, spacing="ends"))
        notes += d.notes
        if len(d.x1) < 2:
            out["ok"] = False
            out["error"] = ("No vapor-liquid equilibrium could be solved at this condition "
                            "with this model. Try a lower temperature or pressure.")
            return out
        out["columns"] = _columns(es.diagram_properties(L, V, d))
        out["fixed"] = {"T": float(req["T"])} if hold == "T" else {"P": float(req["P"])}
        if d.azeotrope:
            xa, va = d.azeotrope
            out["azeotrope"] = {"x1": float(f"{xa:.5g}"),
                                "P" if hold == "T" else "T":
                                float(f"{(va / 1e5 if hold == 'T' else va):.6g}")}
    else:
        x1 = float(req["x1"])
        env = es.pt_envelope(L, V, [x1, 1 - x1], npts=npts)
        notes += env.notes
        out["columns"] = _columns(es.envelope_properties(L, V, env))
        out["fixed"] = {"x1": x1}
    stab = out["columns"].get("stability") or []
    bad = [x for x, v in zip(out["columns"]["x1"], stab) if v is not None and v < 0]
    out["unstable"] = None
    if bad and hold in ("T", "P"):
        lo_x, hi_x = min(bad), max(bad)
        out["unstable"] = [lo_x, hi_x]
        cols = out["columns"]
        for k, x in enumerate(cols["x1"]):          # leave that part of the curves out
            if lo_x <= x <= hi_x:
                for name in cols:
                    if name not in ("x1", "stability"):
                        cols[name][k] = None
        out["azeotrope"] = (None if out["azeotrope"] and lo_x <= out["azeotrope"]["x1"] <= hi_x
                            else out["azeotrope"])
        notes = [n for n in notes if not n.startswith("azeotrope at")]
        notes.append(f"With this model the liquid is unstable between x = {lo_x:.2f} and "
                     f"{hi_x:.2f}: it would split into two liquid phases, which this tool "
                     f"does not calculate yet. That part of the curves is left out.")
    if req.get("model") in ACTIVITY_MODELS:
        lo, hi = pair_defaults(*req["components"])["ranges"]["antoine_T"]
        Ts = [t for t in out["columns"]["T"] if t is not None]
        if Ts and (min(Ts) < lo - 10 or max(Ts) > hi + 10):
            notes.append(f"Vapor pressures are extrapolated: the Antoine constants of this "
                         f"pair are fitted from {lo:.0f} to {hi:.0f} K.")
    out["notes"] = notes
    out["ok"] = True
    return out


def compute(req):
    """Run one request. Returns a dict with 'ok', the input text, the checker
    report, and (when ok) 'columns' with one list per variable."""
    t0 = time.time()
    try:
        out = _solve(req)
        cmp_model = req.get("compare")
        if out.get("ok") and cmp_model and cmp_model != req.get("model"):
            req2 = dict(req, model=cmp_model, params=req.get("compare_params") or {},
                        compare=None)
            try:
                o2 = _solve(req2)
                if o2.get("ok"):
                    out["compare"] = {"model": cmp_model, "columns": o2["columns"],
                                      "azeotrope": o2.get("azeotrope")}
            except Exception as e:                      # the comparison is optional
                out["notes"].append(f"comparison with {cmp_model} failed: {e}")
    except Exception as e:
        out = {"ok": False, "error": str(e) or type(e).__name__}
        try:
            out["input_text"] = build_input(req)
        except Exception:
            pass
    out["ms"] = round(1000 * (time.time() - t0))
    return out


def handle(message):
    """JSON in, JSON out: {"op": "catalog"} | {"op": "pair", "a":, "b":} |
    {"op": "compute", "request": {...}}."""
    msg = json.loads(message) if isinstance(message, str) else message
    op = msg.get("op")
    try:
        if op == "catalog":
            res = catalog()
        elif op == "pair":
            res = pair_defaults(msg["a"], msg["b"])
        elif op == "compute":
            res = compute(msg["request"])
        else:
            res = {"ok": False, "error": f"unknown operation '{op}'"}
    except Exception as e:
        res = {"ok": False, "error": str(e) or type(e).__name__}
    return json.dumps(res)


# --------------------------------------------------------------------------
# Self-test
# --------------------------------------------------------------------------
if __name__ == "__main__":
    ok_all = True

    def report(name, ok, detail=""):
        global ok_all
        ok_all &= bool(ok)
        print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f": {detail}" if detail else ""))

    cat = catalog()
    report("catalog", len(cat["compounds"]) >= 20 and len(cat["models"]) == 7,
           f"{len(cat['compounds'])} compounds, {len(cat['pairs'])} pairs with parameters")

    r = compute({"components": ["ethanol", "water"], "model": "nrtl", "hold": "P", "P": 1.01325})
    az = r.get("azeotrope")
    report("ethanol-water NRTL azeotrope near the measured one (x = 0.894, 351.3 K)",
           r["ok"] and az and abs(az["x1"] - 0.894) < 0.03 and abs(az["T"] - 351.3) < 0.5,
           f"x1 = {az['x1']}, T = {az['T']} K, {r['ms']} ms" if az else str(r.get("error")))
    r = compute({"components": ["ethanol", "water"], "model": "nrtl", "hold": "T", "T": 350})
    c = r["columns"]
    k = len(c["x1"]) // 2
    lhs = c["y1"][k] * c["P"][k]                            # vapor side: y1 P
    rhs = c["x1"][k] * c["gamma1"][k] * c["f1"][-1]         # liquid side: x1 gamma1 f1(pure)
    report("equal fugacities: f1 = y1 P (vapor) = x1 gamma1 f1,pure (liquid)",
           abs(c["f1"][k] / lhs - 1) < 1e-6 and abs(lhs / rhs - 1) < 1e-6,
           f"f1 = {c['f1'][k]:.5f} bar, y1 P = {lhs:.5f} bar, x1 gamma1 f1,pure = {rhs:.5f} bar")

    r = compute({"components": ["acetone", "chloroform"], "model": "nrtl", "hold": "T", "T": 323.15})
    az = r.get("azeotrope")
    report("acetone-chloroform: negative deviation, maximum-boiling azeotrope",
           r["ok"] and az and max(r["columns"]["gamma1"]) <= 1.0 + 1e-9
           and az["P"] < min(r["columns"]["P"][0], r["columns"]["P"][-1]),
           f"x1 = {az['x1']}, P = {az['P']} bar" if az else str(r.get("error")))

    r = compute({"components": ["propane", "n-butane"], "model": "PR", "hold": "T", "T": 300,
                 "compare": "SRK"})
    report("PR diagram with SRK comparison", r["ok"] and "compare" in r,
           f"P ends {r['columns']['P'][0]:.3f} to {r['columns']['P'][-1]:.3f} bar")
    g = [v for v in r["columns"]["gamma1"] if v is not None]
    report("EOS route also reports activity coefficients (near 1 for this pair)",
           len(g) > 30 and all(0.9 < v < 1.2 for v in g), f"gamma1 from {min(g):.3f} to {max(g):.3f}")

    r = compute({"components": ["propane", "n-butane"], "model": "PR", "hold": "x", "x1": 0.5})
    c = r["columns"] if r["ok"] else {}
    report("P-T diagram at fixed composition: bubble P above dew P",
           r["ok"] and all(pb > pd for pb, pd in zip(c["P"], c["P_dew"]) if pb and pd),
           f"{len(c.get('T', []))} temperatures")

    d = pair_defaults("ethanol", "water")
    r1 = compute({"components": ["ethanol", "water"], "model": "vanlaar", "hold": "P", "P": 1.01325})
    report("van Laar defaults matched to NRTL at infinite dilution give an azeotrope",
           r1["ok"] and r1.get("azeotrope") is not None,
           f"A12 = {d['A']['A12']}, A21 = {d['A']['A21']}, azeotrope {r1.get('azeotrope')}")

    r = compute({"components": ["n-hexane", "n-heptane"], "model": "ideal", "hold": "T", "T": 340})
    c = r["columns"]
    report("ideal solution: gamma = 1 and no azeotrope",
           r["ok"] and r["azeotrope"] is None and all(abs(v - 1) < 1e-12 for v in c["gamma1"]))

    r = compute({"components": ["ethanol", "water"], "model": "vanlaar", "hold": "T", "T": 350,
                 "params": {"A12": 1.0, "A21": -1.0}})
    report("invalid parameters give a readable error, not a crash",
           not r["ok"] and "same sign" in r["error"], r.get("error", ""))

    rt = json.loads(handle(json.dumps({"op": "compute", "request": {
        "components": ["methanol", "water"], "model": "nrtl", "hold": "P", "P": 1.0}})))
    report("JSON round trip", rt["ok"] and "NaN" not in json.dumps(rt))

    print("\nALL TESTS PASSED" if ok_all else "\nSOME TESTS FAILED")
