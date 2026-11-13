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
import os
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
    pairs = [{"a": p["a"], "b": p["b"], "nrtl": "nrtl" in p, "kij_pr": "kij_pr" in p,
              "vanlaar": "vanlaar" in p}
             for p in db["pairs"]]
    return {"compounds": comps, "pairs": pairs, "sources": db["meta"]["sources"],
            "generated": db["meta"]["generated"],
            "models": [{"id": m, "label": MODEL_LABELS[m],
                        "family": "eos" if m in EOS_MODELS else "activity"}
                       for m in EOS_MODELS + ACTIVITY_MODELS]}


def pair_defaults(a, b):
    """Starting parameter values and sensible ranges for a pair of compounds.

    NRTL parameters, the Peng-Robinson k_ij and the temperature-dependent van
    Laar parameters A_ij(T) = a_ij + b_ij/T come from the database when it has
    them ('vanlaar': A12 and A21 are given at T_ref, with the slopes b12, b21).
    Margules has no database; its A12 and A21 are ln(gamma) at infinite
    dilution, matched to the NRTL model's values at the mean boiling point
    when NRTL parameters exist, and zero (an ideal solution) otherwise."""
    ca, cb = compound_db.find(a), compound_db.find(b)
    if ca is None or cb is None:
        raise ValueError("unknown compound")
    pr = compound_db.pair(a, b)
    T_ref = 0.5 * (ca["Tb"]["value"] + cb["Tb"]["value"])
    out = {"a": ca["id"], "b": cb["id"], "T_ref": round(T_ref, 2),
           "kij": {"value": pr.get("kij_pr", 0.0),
                   "source": pr.get("kij_source", "none")},
           "nrtl": None, "A": {"A12": 0.0, "A21": 0.0, "source": "none"}, "vanlaar": None}
    if "vanlaar" in pr:
        v = pr["vanlaar"]
        out["vanlaar"] = {"A12": round(v["a_ab"] + v["b_ab"] / T_ref, 4),
                          "A21": round(v["a_ba"] + v["b_ba"] / T_ref, 4),
                          "b12": v["b_ab"], "b21": v["b_ba"], "T_ref": round(T_ref, 2),
                          "Tmin": v["Tmin"], "Tmax": v["Tmax"],
                          "rms_lngamma": v["rms_lngamma"], "source": v["source"]}
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
        if model == "vanlaar":
            v = d["vanlaar"]
            if v is None and not ("A12" in par and "A21" in par):
                raise ValueError("the database has no van Laar parameters for this pair")
            A12 = float(par.get("A12", v["A12"] if v else 0.0))
            A21 = float(par.get("A21", v["A21"] if v else 0.0))
            if A12 * A21 <= 0:
                raise ValueError("van Laar needs A12 and A21 of the same sign")
            if v is not None and par.get("tdep", 1):
                # A_ij(T) = a_ij + b_ij/T, passing through the chosen value at T_ref
                Tr = v["T_ref"]
                lines.append(f"A({a}, {b}) = {compound_db.vanlaar_expr(A12 - v['b12'] / Tr, v['b12'])}")
                lines.append(f"A({b}, {a}) = {compound_db.vanlaar_expr(A21 - v['b21'] / Tr, v['b21'])}")
            else:
                lines.append(f"A({a}, {b}) = {A12:g}")
                lines.append(f"A({b}, {a}) = {A21:g}")
        elif model == "margules":
            lines.append(f"A({a}, {b}) = {float(par.get('A12', d['A']['A12'])):g}")
            lines.append(f"A({b}, {a}) = {float(par.get('A21', d['A']['A21'])):g}")
        elif model == "nrtl":
            if d["nrtl"] is None and not ("dg12" in par and "dg21" in par):
                raise ValueError("the database has no NRTL parameters for this pair")
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
    elif hold == "field":
        lines += _field_lines(req, a, b, model)
    else:
        raise ValueError("hold must be 'T', 'P', 'x' or 'field'")
    return "\n".join(lines) + "\n"


# --------------------------------------------------------------------------
# Compute
# --------------------------------------------------------------------------
def rackett_VL(name, T):
    """Saturated-liquid molar volume (m3/mol) from the Rackett equation,
    V = (R Tc/Pc) Z_RA^[1 + (1 - Tr)^(2/7)], with the Yamada-Gunn estimate
    Z_RA = 0.29056 - 0.08775 omega (needs only Tc, Pc, omega from the
    database). Used for the pressure dependence of the liquid in a field
    on the activity route."""
    rec = compound_db.find(name)
    Tc, Pc, w = rec["Tc"]["value"], rec["Pc"]["value"], rec["omega"]["value"]
    Tr = min(T / Tc, 0.99)
    return 8.314462618 * Tc / Pc * (0.29056 - 0.08775 * w) ** (1 + (1 - Tr) ** (2 / 7))


def _field_lines(req, a, b, model):
    """Input lines for a column in a field: the state (T, P, feed z) at the
    reference position (the bottom of a gravity column, the inner radius of
    a centrifuge) and the field itself."""
    f = req["field"]
    T, P, z1 = float(req["T"]), float(req["P"]), float(req["z1"])
    lines = [f"T = {T:g} K", f"P = {P:g} bar", f"z = {z1:g}"]
    if f["kind"] == "gravity":
        H = float(f["height"])
        lines += ["field = gravity", f"g = {float(f['g']):g} m/s2", f"height = 0, {H:g} m",
                  "reference = 0 m"]
    elif f["kind"] == "centrifugal":
        r1, r2 = float(f["r_in"]), float(f["r_out"])
        if not 0 <= r1 < r2:
            raise ValueError("the outer radius must be larger than the inner radius")
        lines += ["field = centrifugal", f"omega = {float(f['rpm']):g} rpm",
                  f"radius = {r1:g}, {r2:g} m", f"reference = {r1:g} m"]
    else:
        raise ValueError("field must be 'gravity' or 'centrifugal'")
    if model in ACTIVITY_MODELS:
        for c in (a, b):
            lines.append(f"{c}.V_L = {rackett_VL(c, T) * 1e6:.4g} cm3/mol   "
                         f"# Rackett equation (Yamada-Gunn), at {T:g} K")
    return lines


def _profile(req, text, out):
    """Equilibrium at the reference position (a flash at T, P, z), then the
    profile of pressure, composition and phase along the field."""
    system = build_system(text)
    eq = next(e for e in system.equilibria if e.phase_names[1] == "V")
    L, V = eq.phases
    z1, P = float(req["z1"]), float(req["P"]) * 1e5
    point = es.flash_TP(L, V, float(req["T"]), P, np.array([z1, 1 - z1]))
    prof = es.column_profile(system, L, V, point, npts=int(req.get("npts", 61)),
                             time_limit=float(req.get("time_limit", 8.0)))
    out.update(route=eq.route, calculation=eq.calculation, names=list(system.names),
               hold="field", azeotrope=None, three_phase=None, unstable=None,
               field=dict(req["field"]), fixed={"T": float(req["T"])})
    out["columns"] = {"position": _clean(prof.position), "P": _clean(prof.P / 1e5),
                      "x1": _clean(prof.x[:, 0]), "rho": _clean(prof.density),
                      "phase": list(prof.phase)}
    fin = lambda v: float(v) if np.isfinite(v) else None
    out["reference"] = {"phase": point.phase, "beta": fin(point.beta),
                        "x1": fin(point.x[0]), "y1": fin(point.y[0])}
    out["contacts"] = [{"position": float(c["position"]), "P": float(c["P"]) / 1e5,
                        "x_liquid": float(c["x_liquid"][0]), "x_vapor": float(c["x_vapor"][0])}
                       for c in prof.extra.get("contacts", [])]
    # saturation pressure of the local fluid (bubble P of a liquid, dew P of a
    # vapor) at every third point: where the pressure line meets it, the
    # other phase appears
    T = float(req["T"])
    psat = [None] * len(prof.position)
    for i in list(range(0, len(prof.position), 3)) + [len(prof.position) - 1]:
        try:
            psat[i] = es.saturation(L, V, prof.x[i], 0.0 if prof.phase[i] == "liquid" else 1.0, T=T).P / 1e5
        except Exception:
            pass
    out["columns"]["P_sat"] = [None if v is None else float(f"{v:.7g}") for v in psat]
    out["notes"] = [n for n in prof.notes if not n.startswith("single phase")]
    out["ok"] = True
    return out


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
    if req.get("hold") == "field":
        return _profile(req, text, out)
    system = build_system(text)
    eq = next(e for e in system.equilibria if e.phase_names[1] == "V")
    L, V = eq.phases
    npts = int(req.get("npts", 41))
    hold = req.get("hold", "T")
    notes = []
    out.update(route=eq.route, calculation=eq.calculation, names=list(system.names),
               hold=hold, azeotrope=None, three_phase=None)
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
        if d.three_phase:
            tp = d.three_phase
            val = tp["value"] / 1e5 if hold == "T" else tp["value"]
            out["three_phase"] = {"x_alpha": float(f"{tp['x_alpha']:.5g}"),
                                  "x_beta": float(f"{tp['x_beta']:.5g}"),
                                  "y": float(f"{tp['y']:.5g}"),
                                  "P" if hold == "T" else "T": float(f"{val:.6g}"),
                                  "heteroazeotrope": tp["heteroazeotrope"]}
            # a row in the middle of the gap: the bubble and dew curves pass
            # through it (the three-phase line), every other column is empty,
            # which breaks those curves across the two-liquid region
            cols = out["columns"]
            k = sum(1 for x in cols["x1"] if x <= tp["x_alpha"] * (1 + 1e-5) + 1e-12)
            mid = 0.5 * (tp["x_alpha"] + tp["x_beta"])
            for name in cols:
                keep = {"x1": mid, "y1": cols["y1"][k], "T": cols["T"][k], "P": cols["P"][k]}
                cols[name].insert(k, float(f"{keep[name]:.7g}") if name in keep else None)
            notes = [n for n in notes if not n.startswith("the liquid splits")]
            notes.append(f"This model predicts two liquid phases between x = {tp['x_alpha']:.3f} "
                         f"and {tp['x_beta']:.3f}. Along the flat part of the bubble curve, both "
                         f"liquids and a vapor of y = {tp['y']:.3f} coexist.")
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
    if hold == "x" and bad:
        notes.append("With this model a liquid of this composition is unstable and would "
                     "split into two liquid phases, so its bubble curve is not physical. "
                     "Choose a composition outside the two-liquid range.")
    if bad and hold in ("T", "P") and not out["three_phase"]:
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
    if req.get("model") == "vanlaar" and (req.get("params") or {}).get("tdep", 1):
        v = pair_defaults(*req["components"])["vanlaar"]
        Ts = [t for t in out["columns"]["T"] if t is not None]
        if v and Ts and (min(Ts) < v["Tmin"] - 5 or max(Ts) > v["Tmax"] + 5):
            notes.append(f"The van Laar parameters are extrapolated: their temperature "
                         f"dependence is fitted from {v['Tmin']:.0f} to {v['Tmax']:.0f} K.")
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


def surface(req):
    """The same diagram at several values of the variable held fixed (req
    'values': temperatures in K for hold = 'T', pressures in bar for 'P'):
    the layers of a 3D surface, e.g. Txy diagrams stacked along pressure.
    Layers that cannot be solved are left out."""
    t0 = time.time()
    hold = req.get("hold", "T")
    if hold not in ("T", "P"):
        raise ValueError("a surface stacks Pxy or Txy diagrams: hold must be 'T' or 'P'")
    keep = req.get("keep")
    layers, failed = [], 0
    for v in req["values"]:
        r = {k: val for k, val in req.items() if k not in ("values", "keep", "compare")}
        r[hold] = float(v)
        try:
            o = _solve(r)
        except Exception:
            o = {"ok": False}
        if not o.get("ok"):
            failed += 1
            continue
        cols = o["columns"] if not keep else {k: o["columns"][k] for k in keep if k in o["columns"]}
        layers.append({"value": float(v), "columns": cols, "azeotrope": o.get("azeotrope")})
    out = {"ok": bool(layers), "hold": hold, "layers": layers, "failed": failed,
           "ms": round(1000 * (time.time() - t0))}
    if not layers:
        out["error"] = "No layer of the surface could be calculated in this range."
    return out


def point(req):
    """Is there vapor-liquid equilibrium at the user's T, P and overall
    composition z1? A flash answers it; the bubble and dew pressures (at T)
    and temperatures (at P) of the same feed say where equilibrium occurs."""
    t0 = time.time()
    try:
        r = dict(req, hold="T", npts=3)
        lines = build_input(r).splitlines()
        lines = [l for l in lines if not l.startswith(("sweep", "T =", "P =", "x(L)"))]
        T, P, z1 = float(req["T"]), float(req["P"]), float(req["z1"])
        if not 0 < z1 < 1:
            raise ValueError("the overall mole fraction must be between 0 and 1")
        text = "\n".join(lines + [f"T = {T:g} K", f"P = {P:g} bar", f"z = {z1:g}"]) + "\n"
        rep = check_input(text)
        out = {"input_text": text, "report": str(rep), "status": rep.status}
        if rep.status != "SUFFICIENT":
            out.update(ok=False, error="The checker did not accept this input: " + rep.status)
            return out
        system = build_system(text)
        eq = next(e for e in system.equilibria if e.phase_names[1] == "V")
        L, V = eq.phases
        z = np.array([z1, 1 - z1])
        pt = es.flash_TP(L, V, T, P * 1e5, z)
        fin = lambda v: None if v is None or not np.isfinite(v) else float(f"{v:.6g}")
        out.update(ok=True, route=eq.route, phase=pt.phase, beta=fin(pt.beta),
                   x1=fin(pt.x[0]), y1=fin(pt.y[0]), T=T, P=P, z1=z1)
        if pt.phase == "two-phase":
            K = es.K_values(L, V, T, P * 1e5, pt.x, pt.y)
            g = L.ln_gamma(T, pt.x) if hasattr(L, "ln_gamma") else None
            out["K"] = [fin(k) for k in K]
            out["gamma"] = None if g is None else [fin(math.exp(v)) for v in g]
        bounds = {}
        for name, beta, kw in (("P_bubble", 0.0, {"T": T}), ("P_dew", 1.0, {"T": T}),
                               ("T_bubble", 0.0, {"P": P * 1e5}), ("T_dew", 1.0, {"P": P * 1e5})):
            try:
                s_ = es.saturation(L, V, z, beta, **kw)
                bounds[name] = fin(s_.P / 1e5) if name.startswith("P") else fin(s_.T)
            except Exception:
                bounds[name] = None
        out["bounds"] = bounds
    except Exception as e:
        out = {"ok": False, "error": str(e) or type(e).__name__}
    out["ms"] = round(1000 * (time.time() - t0))
    return out


def run_input(text, time_limit=20.0):
    """Check any input file and, if it is sufficient, solve it: exactly what
    `python equilibrium_solver.py file.txt` prints."""
    import contextlib
    import io
    t0 = time.time()
    rep = check_input(text)
    out = {"ok": True, "status": rep.status, "report": str(rep), "results": ""}
    if rep.status == "SUFFICIENT":
        buf = io.StringIO()
        es._deadline[0] = time.time() + time_limit
        try:
            with contextlib.redirect_stdout(buf):
                system = build_system(text)
                es.print_results(system, es.solve(system))
        except Exception as e:
            buf.write(f"\nThe calculation stopped: {e}\n")
        finally:
            es._deadline[0] = None
        out["results"] = buf.getvalue()
    out["ms"] = round(1000 * (time.time() - t0))
    return out


def handle(message):
    """JSON in, JSON out: {"op": "ping"} | {"op": "catalog"} |
    {"op": "pair", "a":, "b":} | {"op": "compute", "request": {...}} |
    {"op": "surface", "request": {..., "values": [...]}}."""
    msg = json.loads(message) if isinstance(message, str) else message
    op = msg.get("op")
    try:
        if op == "ping":
            res = {"ok": True, "engine": "python"}
        elif op == "catalog":
            res = catalog()
        elif op == "pair":
            res = pair_defaults(msg["a"], msg["b"])
        elif op == "compute":
            res = compute(msg["request"])
        elif op == "surface":
            res = surface(msg["request"])
        elif op == "point":
            res = point(msg["request"])
        elif op == "run_input":
            res = run_input(msg["text"])
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

    A_ = 2.6
    r = compute({"components": ["ethanol", "water"], "model": "margules", "hold": "T", "T": 350,
                 "params": {"A12": A_, "A21": A_}})
    tp = r.get("three_phase") or {}
    xa = tp.get("x_alpha", 0.5)
    exact = abs(math.log(xa / (1 - xa)) - A_ * (2 * xa - 1))       # symmetric Margules binodal
    report("a model that predicts two liquids gives the three-phase line (checked against "
           "the closed-form Margules result)",
           r["ok"] and tp and exact < 1e-3 and abs(tp["x_alpha"] + tp["x_beta"] - 1) < 1e-4
           and r["azeotrope"] is None and tp["heteroazeotrope"],
           f"x = {tp.get('x_alpha')} and {tp.get('x_beta')}, vapor y = {tp.get('y')}")

    r = compute({"components": ["n-hexane", "n-heptane"], "model": "ideal", "hold": "T", "T": 340})
    c = r["columns"]
    report("ideal solution: gamma = 1 and no azeotrope",
           r["ok"] and r["azeotrope"] is None and all(abs(v - 1) < 1e-12 for v in c["gamma1"]))

    r = compute({"components": ["ethanol", "water"], "model": "vanlaar", "hold": "T", "T": 350,
                 "params": {"A12": 1.0, "A21": -1.0}})
    report("invalid parameters give a readable error, not a crash",
           not r["ok"] and "same sign" in r["error"], r.get("error", ""))

    # field: the 500 m propane / n-butane column forms a gas cap near 306 m
    r = compute({"components": ["propane", "n-butane"], "model": "PR", "hold": "field", "T": 320,
                 "P": 25, "z1": 0.4, "npts": 41,
                 "field": {"kind": "gravity", "g": 9.80665, "height": 500}})
    report("gravity column in the calculator: gas-liquid contact", r["ok"] and len(r["contacts"]) == 1
           and abs(r["contacts"][0]["position"] - 306.1) < 1.5 and r["columns"]["phase"][-1] == "vapor",
           str(r.get("contacts")))
    r = compute({"components": ["ethanol", "water"], "model": "vanlaar", "hold": "field", "T": 350,
                 "P": 1.0, "z1": 0.3, "npts": 21,
                 "field": {"kind": "centrifugal", "rpm": 20000, "r_in": 0.05, "r_out": 0.10}})
    report("centrifuge with an activity model: Rackett liquid volumes, P rises outward",
           r["ok"] and "Rackett" in r["input_text"] and r["columns"]["P"][-1] > r["columns"]["P"][0] > 0.99,
           f"{r['columns']['P'][-1] if r.get('ok') else r.get('error')} bar at the rim")

    # "is there equilibrium here?": a state inside the two-phase region, and one above it
    r = point({"components": ["ethanol", "water"], "model": "vanlaar", "T": 360.0, "P": 1.01325, "z1": 0.4})
    r2 = point({"components": ["ethanol", "water"], "model": "vanlaar", "T": 360.0, "P": 2.0, "z1": 0.4})
    b_ = r.get("bounds", {})
    report("state check: two phases between the dew and bubble pressures, liquid above",
           r["ok"] and r["phase"] == "two-phase" and b_["P_dew"] < 1.01325 < b_["P_bubble"]
           and r["x1"] < 0.4 < r["y1"] and r2["phase"] == "liquid",
           f"x = {r.get('x1')}, y = {r.get('y1')}, beta = {r.get('beta')}")
    o = run_input(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "examples",
                                    "gravity_column_gas_cap.txt")).read())
    report("own input: checked, solved and printed", o["status"] == "SUFFICIENT"
           and "gas-liquid contact" in o["results"], f"{o['ms']} ms")

    # 3D: each layer of a surface is the ordinary diagram at that pressure
    sv = surface({"components": ["ethanol", "water"], "model": "nrtl", "hold": "P",
                  "values": [0.5, 1.01325, 2.0], "npts": 21})
    one = compute({"components": ["ethanol", "water"], "model": "nrtl", "hold": "P",
                   "P": 1.01325, "npts": 21})
    report("surface layers equal the single diagrams", sv["ok"] and len(sv["layers"]) == 3
           and sv["layers"][1]["columns"]["T"] == one["columns"]["T"]
           and sv["layers"][0]["columns"]["T"][-1] < sv["layers"][2]["columns"]["T"][-1],
           f"{sv['ms']} ms")

    # database van Laar: A_ij(T) = a + b/T, against the NRTL fit it came from
    pdv = pair_defaults("ethanol", "water")["vanlaar"]
    r = compute({"components": ["ethanol", "water"], "model": "vanlaar", "hold": "P", "P": 1.01325})
    rn = compute({"components": ["ethanol", "water"], "model": "nrtl", "hold": "P", "P": 1.01325})
    report("van Laar from the database: ethanol-water azeotrope near the NRTL one",
           r["ok"] and abs(r["azeotrope"]["x1"] - rn["azeotrope"]["x1"]) < 0.03
           and abs(r["azeotrope"]["T"] - rn["azeotrope"]["T"]) < 0.3
           and "/T" in r["input_text"] and pdv["rms_lngamma"] < 0.01,
           f"{r['azeotrope']} vs {rn['azeotrope']}")
    g = [compute({"components": ["ethanol", "water"], "model": "vanlaar", "hold": "T", "T": T,
                  "params": {"tdep": td}})["columns"]["gamma1"][0] for T in (340.0, 370.0) for td in (1, 0)]
    report("temperature dependence: gamma(inf) changes with T only when it is switched on",
           abs(g[0] - g[2]) > 0.1 and abs(g[1] - g[3]) < 1e-9, f"{g}")
    r = compute({"components": ["propane", "n-butane"], "model": "vanlaar", "hold": "T", "T": 300})
    report("a pair without database van Laar parameters is refused", not r["ok"]
           and "no van Laar" in r["error"], r.get("error", ""))

    rt = json.loads(handle(json.dumps({"op": "compute", "request": {
        "components": ["methanol", "water"], "model": "nrtl", "hold": "P", "P": 1.0}})))
    report("JSON round trip", rt["ok"] and "NaN" not in json.dumps(rt))

    print("\nALL TESTS PASSED" if ok_all else "\nSOME TESTS FAILED")
