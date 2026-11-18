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
import pure_fluid

EOS_MODELS = ("PR", "SRK", "vdW")
ACTIVITY_MODELS = ("ideal", "margules", "nrtl")
MODEL_LABELS = {"PR": "Peng-Robinson", "SRK": "Soave-Redlich-Kwong", "vdW": "van der Waals",
                "ideal": "Ideal solution (Raoult's law)", "margules": "Margules",
                "nrtl": "NRTL"}


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
            "Tb": c["Tb"]["value"] if c["Tb"] else None,
            "antoine": ({"A": a["A"], "B": a["B"], "C": a["C"], "Tmin": a["Tmin"],
                         "Tmax": a["Tmax"], "form": a["form"]} if a else None),
            "ppr78": (c["ppr78"]["groups"] if "ppr78" in c else None),
            "sources": {"M": c["M"]["source"], "Tc": c["Tc"]["source"],
                        "Pc": c["Pc"]["source"], "omega": c["omega"]["source"],
                        "antoine": a["source"] if a else None,
                        "Tb": c["Tb"]["source"] if c["Tb"] else None,
                        "ppr78": "eppr78" if "ppr78" in c else None},
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

    NRTL parameters come from the database when it has them. The Peng-Robinson
    k_ij comes from E-PPR78 (a function of T; 'kij' then reports its value at
    T_ref) or else the fixed DECHEMA value. Margules has no database; its A12 and A21 are ln(gamma) at infinite
    dilution, matched to the NRTL model's values at the mean boiling point
    when NRTL parameters exist, and zero (an ideal solution) otherwise."""
    ca, cb = compound_db.find(a), compound_db.find(b)
    if ca is None or cb is None:
        raise ValueError("unknown compound")
    pr = compound_db.pair(a, b)
    Tb_a, Tb_b = (_tb(c) for c in (ca, cb))
    T_ref = 0.5 * (Tb_a + Tb_b)
    if compound_db.ppr78_ready(a, b):
        kij = {"value": round(compound_db.PPR78Kij(a, b)(T_ref), 4), "source": "eppr78",
               "T_dependent": True}
    else:
        kij = {"value": pr.get("kij_pr", 0.0), "source": pr.get("kij_source", "none"),
               "T_dependent": False}
    out = {"a": ca["id"], "b": cb["id"], "T_ref": round(T_ref, 2), "kij": kij,
           "nrtl": None, "A": {"A12": 0.0, "A21": 0.0, "source": "none"}}
    if "dg_ab" in pr:
        out["nrtl"] = {"dg12": round(pr["dg_ab"], 2), "dg21": round(pr["dg_ba"], 2),
                       "alpha": pr["alpha"], "source": pr["nrtl_source"]}
        m = NRTL(np.array([[0, pr["dg_ab"]], [pr["dg_ba"], 0]]),
                 np.array([[0, pr["alpha"]], [pr["alpha"], 0]]))
        A12 = float(m.ln_gamma(T_ref, [0.0, 1.0])[0])
        A21 = float(m.ln_gamma(T_ref, [1.0, 0.0])[1])
        out["A"] = {"A12": round(A12, 4), "A21": round(A21, 4), "source": "nrtl_infinite_dilution"}
    Tb_lo, Tb_hi = sorted((Tb_a, Tb_b))
    Tc_lo = min(ca["Tc"]["value"], cb["Tc"]["value"])
    out["ranges"] = {
        "T": [int(5 * math.floor((Tb_lo - 60) / 5)), int(5 * math.ceil(min(Tb_hi + 80, Tc_lo - 5) / 5))],
        "T_default": int(5 * round(T_ref / 5)),
        "P": [0.05, 30.0], "P_default": 1.01325,
        "antoine_T": ([max(ca["antoine"]["Tmin"], cb["antoine"]["Tmin"]),
                       min(ca["antoine"]["Tmax"], cb["antoine"]["Tmax"])]
                      if ca["antoine"] and cb["antoine"] else None),
    }
    return out


def _tb(rec):
    """Normal boiling point, or for a compound without Antoine constants (CO2,
    which sublimes at 1 atm) a stand-in of 0.65 Tc for choosing ranges."""
    return rec["Tb"]["value"] if rec["Tb"] else round(0.65 * rec["Tc"]["value"], 2)


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
        if "kij" in par:
            lines.append(f"kij({a}, {b}) = {float(par['kij']):g}")
        elif model == "PR" and d["kij"]["source"] == "eppr78":
            lines.append(f"kij({a}, {b}) = PPR78")
        else:
            lines.append(f"kij({a}, {b}) = {float(d['kij']['value'] if model == 'PR' else 0.0):g}")
    elif model in ACTIVITY_MODELS:
        lines.append(f"activity = {model}")
        if model == "margules":
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


# --------------------------------------------------------------------------
# Pure fluids
# --------------------------------------------------------------------------
def _pure_field(f, start=None):
    """pure_fluid.PureField from the request's field: gravity over a height
    with the reference at the bottom, or a centrifuge between two radii with
    the reference at the inner radius."""
    if f["kind"] == "gravity":
        H = float(f["height"])
        if H <= 0:
            raise ValueError("the column height must be positive")
        return pure_fluid.PureField("gravity", float(f["g"]), 0.0, H, 0.0)
    if f["kind"] == "centrifugal":
        r1, r2 = float(f["r_in"]), float(f["r_out"])
        if not 0 <= r1 < r2:
            raise ValueError("the outer radius must be larger than the inner radius")
        return pure_fluid.PureField("centrifugal", float(f["rpm"]) * 2 * math.pi / 60, r1, r2, r1)
    raise ValueError("field must be 'gravity' or 'centrifugal'")


def pure(req):
    """One compound with a cubic equation of state: the vapor-pressure curve
    to the critical point, the coexisting densities, Clapeyron's enthalpy of
    vaporization, isotherms on reduced P-V axes with the equal-area line, the
    state at (T, P), and the Antoine vapor pressure for comparison. With a
    'field', also the column at (T, P at the reference position), the band
    of reference pressures over which the column holds both phases at four
    field strengths, and the interface position against field strength."""
    t0 = time.time()
    old_err = np.seterr(all="ignore")
    try:
        model = req.get("model", "PR")
        if model not in EOS_MODELS:
            raise ValueError("a pure fluid needs an equation of state (PR, SRK or vdW)")
        eos, M, rec = pure_fluid.make_eos(req["compound"], model)
        Tc, Pc, vc = pure_fluid.eos_critical_point(eos)
        cur = pure_fluid.saturation_curve(eos, npts=int(req.get("npts", 60)))
        fin = lambda v: None if v is None or not np.isfinite(v) else float(f"{v:.6g}")
        vec = lambda a: [fin(v) for v in a]
        a = rec["antoine"]
        P_ant = [fin(10 ** (a["A"] - a["B"] / (T + a["C"])) / 1e5) if a and a["Tmin"] - 1e-9 <= T <= a["Tmax"] + 1e-9
                 else None for T in cur.T] if a else [None] * len(cur.T)
        out = {"ok": True, "compound": rec["id"], "name": rec["name"], "model": model,
               "M": M * 1000, "Tc": Tc, "Pc": Pc / 1e5, "rhoc": fin(M / vc),
               "Tb_data": rec["Tb"]["value"] if rec["Tb"] else None,
               "antoine_range": [a["Tmin"], a["Tmax"]] if a else None,
               "curve": {"T": vec(cur.T), "P": vec(cur.P / 1e5), "rhoL": vec(M / cur.vL),
                         "rhoV": vec(M / cur.vV), "dHvap": vec(cur.dHvap / 1000), "P_antoine": P_ant}}
        try:
            out["Tb_eos"] = fin(pure_fluid.saturation_temperature(eos, 101325.0))
        except ValueError:
            out["Tb_eos"] = None
        # isotherms in reduced units
        iso = []
        for tr in (0.8, 0.85, 0.9, 0.95, 1.0, 1.1):
            it = pure_fluid.isotherm(eos, tr * Tc, npts=260, V_max_factor=6)
            keep = it["V"] / vc >= 0.33
            iso.append({"Tr": tr, "V": vec(it["V"][keep] / vc), "P": vec(it["P"][keep] / Pc),
                        "Psat": fin(it["Psat"] / Pc) if it["Psat"] else None,
                        "vL": fin(it["vL"] / vc) if it["vL"] else None,
                        "vV": fin(it["vV"] / vc) if it["vV"] else None,
                        "area": fin(it.get("maxwell_area"))})
        out["isotherms"] = iso
        out["dome"] = {"vL": vec(cur.vL / vc), "vV": vec(cur.vV / vc), "P": vec(cur.P / Pc)}
        # the state asked about
        T, P = float(req.get("T", 0.8 * Tc)), float(req.get("P", 1.01325)) * 1e5
        st = {"T": T, "P": P / 1e5}
        if T < Tc * (1 - 1e-4):
            Ps, vL, vV = pure_fluid.saturation_pressure(eos, T)
            st.update(Psat=fin(Ps / 1e5), rhoL_sat=fin(M / vL), rhoV_sat=fin(M / vV),
                      phase="liquid" if P > Ps else "vapor")
        else:
            st["phase"] = "supercritical fluid" if P > Pc else "gas (above Tc)"
        v = eos.molar_volume(T, P, pure_fluid.ONE, "stable")
        st["rho"] = fin(M / v)
        st["Z"] = fin(P * v / (8.314462618 * T))
        try:
            st["Tsat"] = fin(pure_fluid.saturation_temperature(eos, P)) if P < Pc else None
        except ValueError:
            st["Tsat"] = None
        out["state"] = st
        if req.get("field"):
            out["field"] = _pure_field_results(eos, M, T, P, req["field"], cur)
    except Exception as e:
        out = {"ok": False, "error": str(e) or type(e).__name__}
    finally:
        np.seterr(**old_err)
    out["ms"] = round(1000 * (time.time() - t0))
    return out


def _pure_field_results(eos, M, T, P_ref, f, cur):
    fin = lambda v: None if v is None or not np.isfinite(v) else float(f"{v:.6g}")
    vec = lambda a: [fin(v) for v in a]
    fld = _pure_field(f)
    col = pure_fluid.column(eos, M, T, P_ref, fld, npts=81)
    res = {"kind": fld.kind, "reference": fld.reference,
           "column": {"position": vec(col.position), "P": vec(col.P / 1e5), "rho": vec(col.rho),
                      "phase": col.phase, "Psat": fin(col.Psat / 1e5) if col.Psat else None,
                      "interface": fin(col.interface) if col.interface is not None else None,
                      "notes": col.notes}}
    # the band of reference pressures at four field strengths
    Ts = [t for t in cur.T[:-1] if t > cur.T[0]][::2]
    bands = []
    for frac in (0.25, 0.5, 1.0):
        fs = pure_fluid.PureField(fld.kind, fld.strength * (frac if fld.kind == "gravity" else math.sqrt(frac)),
                                  fld.start, fld.end, fld.reference)
        try:
            Tb, Ps, P1, P2 = pure_fluid.band(eos, M, fs, Ts)
            bands.append({"fraction": frac, "strength": fs.strength, "T": vec(Tb), "P_first": vec(P1 / 1e5),
                          "P_last": vec(P2 / 1e5), "Psat": vec(Ps / 1e5)})
        except Exception:
            pass
    res["bands"] = bands
    # interface position against field strength, at T and P_ref
    if T < eos.Tc[0] * (1 - 1e-4):
        sat = pure_fluid.saturation_pressure(eos, T)
        ln_fs = pure_fluid._ln_f(eos, T, sat[0], "vapor")
        ph = "liquid" if P_ref > sat[0] else "vapor"
        ln_fr = pure_fluid._ln_f(eos, T, P_ref, ph)
        xs, ys = [], []
        for k in range(1, 41):
            s_ = fld.strength * k / 40
            fs = pure_fluid.PureField(fld.kind, s_, fld.start, fld.end, fld.reference)
            phi_star = fs.phi(fs.reference) + 8.314462618 * T * (ln_fr - ln_fs) / M
            xs.append(s_ if fld.kind == "gravity" else s_ * 60 / (2 * math.pi))
            ys.append(fs.position_at(phi_star))
        res["interface_curve"] = {"strength": vec(xs), "position": [fin(v) if v is not None else None for v in ys]}
    return res


# --------------------------------------------------------------------------
# Mixtures of three or four components
# --------------------------------------------------------------------------
MIXTURE_EXAMPLES = [
    {"id": "nat-gas", "label": "Methane, ethane, propane, n-butane (Peng-Robinson)",
     "components": ["methane", "ethane", "propane", "n-butane"], "model": "PR",
     "T": 250.0, "P": 30.0, "z": [0.40, 0.25, 0.20, 0.15]},
    {"id": "co2-gas", "label": "Carbon dioxide, methane, n-butane (Peng-Robinson)",
     "components": ["carbon-dioxide", "methane", "n-butane"], "model": "PR",
     "T": 260.0, "P": 30.0, "z": [0.30, 0.30, 0.40]},
    {"id": "solvents", "label": "Acetone, methanol, water (NRTL)",
     "components": ["acetone", "methanol", "water"], "model": "nrtl",
     "T": 345.0, "P": 1.01325, "z": [0.30, 0.30, 0.40]},
    {"id": "alcohols", "label": "Methanol, ethanol, acetone, water (NRTL)",
     "components": ["methanol", "ethanol", "acetone", "water"], "model": "nrtl",
     "T": 345.0, "P": 1.01325, "z": [0.25, 0.25, 0.20, 0.30]},
]


def mixture_input(req):
    """Input file for a flash of 2 to 4 components at T, P and feed z."""
    comps = [compound_db.find(c) for c in req["components"]]
    if any(c is None for c in comps):
        raise ValueError("unknown compound")
    ids = [c["id"] for c in comps]
    if not 2 <= len(ids) <= 4 or len(set(ids)) != len(ids):
        raise ValueError("choose two to four different compounds")
    model = req.get("model", "PR")
    lines = [f"components = {', '.join(ids)}", "database = builtin"]
    if model in EOS_MODELS:
        lines.append(f"model = {model}")
        for i, a in enumerate(ids):
            for b in ids[i + 1:]:
                if model == "PR" and compound_db.ppr78_ready(a, b):
                    lines.append(f"kij({a}, {b}) = PPR78")
                elif model == "PR" and "kij_pr" in compound_db.pair(a, b):
                    lines.append(f"kij({a}, {b}) = {compound_db.pair(a, b)['kij_pr']:g}")
                else:
                    lines.append(f"kij({a}, {b}) = 0")
    elif model in ("ideal", "nrtl"):
        lines.append(f"activity = {model}")
        missing = [f"{a} + {b}" for i, a in enumerate(ids) for b in ids[i + 1:]
                   if model == "nrtl" and "dg_ab" not in compound_db.pair(a, b)]
        if missing:
            raise ValueError("the database has no NRTL parameters for " + ", ".join(missing))
        for c in comps:
            if c["antoine"] is None:
                raise ValueError(f"{c['name']} has no vapor-pressure data; use an equation of state")
    else:
        raise ValueError(f"unknown model '{model}'")
    z = [float(v) for v in req["z"]]
    if len(z) != len(ids) or any(v <= 0 for v in z):
        raise ValueError("give a positive amount of every component")
    tot = sum(z)
    z = [v / tot for v in z]
    lines += [f"T = {float(req['T']):g} K", f"P = {float(req['P']):g} bar",
              "z = " + ", ".join(f"{v:.6g}" for v in z)]
    return "\n".join(lines) + "\n", ids, np.array(z)


def mixture(req):
    """Flash of a mixture of two to four components at T and P: is it liquid,
    vapor or two-phase, the phase compositions and K-values, and the bubble
    and dew points of the same feed. With 'map': true and three components,
    also the bubble points over a grid of liquid compositions (bubble T at P
    for an activity model, bubble P at T for an equation of state), with the
    vapor in equilibrium with each, for a triangle diagram."""
    t0 = time.time()
    old_err = np.seterr(all="ignore")      # trial points near a critical point overflow harmlessly
    try:
        text, ids, z = mixture_input(req)
        rep = check_input(text)
        out = {"input_text": text, "report": str(rep), "status": rep.status, "names": ids}
        if rep.status != "SUFFICIENT":
            out.update(ok=False, error="The checker did not accept this input: " + rep.status)
            return out
        system = build_system(text)
        eq = next(e for e in system.equilibria if e.phase_names[1] == "V")
        L, V = eq.phases
        T, P = float(req["T"]), float(req["P"])
        fin = lambda v: None if v is None or not np.isfinite(v) else float(f"{v:.6g}")
        vec = lambda a: [fin(v) for v in a]
        pt = es.flash_TP(L, V, T, P * 1e5, z)
        out.update(ok=True, route=eq.route, phase=pt.phase, beta=fin(pt.beta), T=T, P=P,
                   z=vec(z), x=vec(pt.x), y=vec(pt.y))
        if pt.phase == "two-phase":
            out["K"] = vec(es.K_values(L, V, T, P * 1e5, pt.x, pt.y))
            if hasattr(L, "ln_gamma"):
                out["gamma"] = vec(np.exp(L.ln_gamma(T, pt.x)))
            else:
                out["phiL"] = vec(np.exp(L.ln_phi(T, P * 1e5, pt.x)))
                out["phiV"] = vec(np.exp(V.ln_phi(T, P * 1e5, pt.y)))
        kij = {}
        for i, a in enumerate(ids):
            for b in ids[i + 1:]:
                k = rep.inputs.pairs.get(frozenset((a, b)), {}).get("kij")
                if k is not None:
                    kij[f"{a}|{b}"] = fin(k(T) if callable(k) else k)
        if kij:
            out["kij"] = kij
        bounds = {}
        for name, beta, kw in (("P_bubble", 0.0, {"T": T}), ("P_dew", 1.0, {"T": T}),
                               ("T_bubble", 0.0, {"P": P * 1e5}), ("T_dew", 1.0, {"P": P * 1e5})):
            try:
                s_ = es.saturation(L, V, z, beta, **kw)
                bounds[name] = fin(s_.P / 1e5) if name.startswith("P") else fin(s_.T)
            except Exception:
                bounds[name] = None
        out["bounds"] = bounds
        if req.get("map") and len(ids) == 3:
            kind = "P" if req.get("model") in EOS_MODELS else "T"
            out["map"] = _ternary_map(L, V, T, P * 1e5, int(req.get("divisions", 10)), t0, kind)
    except Exception as e:
        out = {"ok": False, "error": str(e) or type(e).__name__}
    finally:
        np.seterr(**old_err)
    out["ms"] = round(1000 * (time.time() - t0))
    return out


def _ternary_map(L, V, T, P, n, t0, kind, limit=15.0):
    """Bubble points on a triangular grid of liquid compositions (n divisions
    per side), with the vapor each liquid boils into: the bubble temperature
    at P (kind 'T', for the activity route) or the bubble pressure at T
    (kind 'P', for an equation of state, where it is the more robust of the
    two). The grid is solved outward from the corner of the heaviest
    compound, each point starting from its solved neighbours."""
    from collections import deque
    pts = {}
    heavy = int(np.argmax(L.eos.Tc)) if hasattr(L, "eos") else None
    if heavy is None:
        heavy = int(np.argmin([math.exp(v) for v in L.ln_f0(T, P)]))
    start = tuple(n if k == heavy else 0 for k in range(3))
    queue, seen = deque([start]), {start}
    truncated = False
    while queue:
        ijk = queue.popleft()
        if time.time() - t0 > limit:
            truncated = True
            break
        x = np.clip(np.array(ijk, float) / n, 1e-6, None)
        x /= x.sum()
        nb = [(ijk[0] + a, ijk[1] + b, ijk[2] - a - b) for a, b in
              ((1, 0), (-1, 0), (0, 1), (0, -1), (1, -1), (-1, 1))]
        nb = [q for q in nb if min(q) >= 0]
        known = [pts[q]["v"] for q in nb if q in pts and pts[q]["v"] is not None]
        guess = (math.exp(sum(math.log(v) for v in known) / len(known)) if known else None)
        rec = {"x": [round(float(v), 4) for v in x], "v": None, "y": None}
        try:
            if kind == "T":
                s_ = es.bubble_T(L, V, P, x, guess=guess)
                val = float(s_.T)
            else:
                s_ = es.bubble_P(L, V, T, x, guess=guess)
                val = float(s_.P)
            if np.all(np.isfinite(s_.y)) and val > 0:
                rec.update(v=val, y=[round(float(v), 4) for v in s_.y])
        except Exception:
            pass
        pts[ijk] = rec
        for q in nb:
            if q not in seen:
                seen.add(q)
                queue.append(q)
    out = []
    for ijk, rec in pts.items():
        v = rec["v"]
        out.append({"x": rec["x"], "y": rec["y"],
                    "T": round(v, 3) if kind == "T" and v is not None else None,
                    "P": round(v / 1e5, 5) if kind == "P" and v is not None else None})
    return {"points": out, "kind": kind, "T": T, "P": P / 1e5, "n": n, "truncated": truncated}


def run_input(text, time_limit=20.0):
    """Check any input file and, if it is sufficient, solve it: exactly what
    `python equilibrium_solver.py file.txt` prints, plus the computed curves
    as chart data ('charts': phase diagrams, P-T envelopes, solubility
    curves of solids and the eutectic, profiles along a field)."""
    import contextlib
    import io
    t0 = time.time()
    rep = check_input(text)
    out = {"ok": True, "status": rep.status, "report": str(rep), "results": "", "charts": []}
    if rep.status == "SUFFICIENT":
        buf = io.StringIO()
        es._deadline[0] = time.time() + time_limit
        old_err = np.seterr(all="ignore")
        try:
            if len(rep.inputs.components) == 1:                 # a pure fluid
                txt, out["charts"] = pure_fluid.solve_input(rep)
                buf.write(txt + "\n")
            else:
                with contextlib.redirect_stdout(buf):
                    system = build_system(text)
                    results = es.solve(system)
                    es.print_results(system, results)
                out["charts"] = _result_charts(system, results)
        except Exception as e:
            buf.write(f"\nThe calculation stopped: {e}\n")
        finally:
            es._deadline[0] = None
            np.seterr(**old_err)
        out["results"] = buf.getvalue()
    out["ms"] = round(1000 * (time.time() - t0))
    return out


def _result_charts(system, results):
    """Chart data for the results of es.solve(): curves that share axes go on
    one chart (a Txy diagram, the solubility curves of both solids and the
    eutectic make one complete T-x diagram)."""
    fin = lambda v: None if v is None or not np.isfinite(v) else float(f"{v:.6g}")
    charts = {}

    def chart(key, title, xt, yt):
        if key not in charts:
            charts[key] = {"title": title, "xTitle": xt, "yTitle": yt, "series": [], "points": []}
        return charts[key]

    names = list(system.names) if hasattr(system, "names") else None
    for res in results:
        d = res.data
        nm = names or getattr(d, "names", None) or ["1", "2"]
        if isinstance(d, es.Diagram):
            pxy = d.kind == "Pxy"
            fixed = f"{d.fixed:.6g} K" if pxy else f"{d.fixed / 1e5:.6g} bar"
            v = (d.var / 1e5) if pxy else d.var
            c = chart(("x1", "P" if pxy else "T", fixed), f"{d.kind} diagram at {fixed}",
                      f"Mole fraction of {d.names[0]}", "Pressure (bar)" if pxy else "Temperature (K)")
            c["series"].append({"name": "Bubble curve (liquid)", "role": "liquid",
                                "pts": [[fin(a), fin(b)] for a, b in zip(d.x1, v)]})
            c["series"].append({"name": "Dew curve (vapor)", "role": "vapor",
                                "pts": [[fin(a), fin(b)] for a, b in zip(d.y1, v)]})
            if d.azeotrope:
                c["points"].append({"x": fin(d.azeotrope[0]), "y": fin(d.azeotrope[1] / 1e5 if pxy else d.azeotrope[1]),
                                    "text": "azeotrope"})
            if d.three_phase:
                tp = d.three_phase
                yv = tp.get("P", tp.get("T"))
                yv = yv / 1e5 if pxy and yv else yv
                c["series"].append({"name": "Three-phase line", "role": "third",
                                    "pts": [[fin(tp["x_alpha"]), fin(yv)], [fin(tp["x_beta"]), fin(yv)]]})
        elif isinstance(d, es.Envelope):
            c = chart(("T", "P", "env"), "P–T envelope at fixed composition", "Temperature (K)", "Pressure (bar)")
            c["series"].append({"name": "Bubble pressure", "role": "liquid",
                                "pts": [[fin(a), fin(b / 1e5)] for a, b in zip(d.T, d.P_bubble)]})
            c["series"].append({"name": "Dew pressure", "role": "vapor",
                                "pts": [[fin(a), fin(b / 1e5)] for a, b in zip(d.T, d.P_dew)]})
        elif isinstance(d, es.SLECurve):
            first = nm[0]
            xs = d.x_solid if d.solid == first else 1 - d.x_solid
            c = chart(("x1", "T", f"{d.P / 1e5:.6g} bar"), f"T–x diagram at {d.P / 1e5:.6g} bar",
                      f"Mole fraction of {first}", "Temperature (K)")
            c["title"] = f"Complete T–x diagram at {d.P / 1e5:.6g} bar"
            c["series"].append({"name": f"Solid {d.solid} appears (solubility curve)",
                                "role": "third" if d.solid == first else "solid2",
                                "pts": [[fin(a), fin(b)] for a, b in zip(xs, d.T)]})
        elif isinstance(d, dict) and res.equilibrium == "Eutectic":
            x = d["x"] if d.get("component") == nm[0] else 1 - d["x"]
            for k, c in charts.items():
                if k[1] == "T" and k[0] == "x1":
                    c["points"].append({"x": fin(x), "y": fin(d["T"]), "text": "eutectic"})
                    for se in c["series"]:            # below the eutectic a solubility curve is metastable
                        if se["name"].startswith("Solid"):
                            se["pts"] = [q for q in se["pts"] if q[1] is not None and q[1] >= d["T"] - 1e-6]
                    c["series"].append({"name": "Eutectic line", "role": "muted",
                                        "pts": [[0, fin(d["T"])], [1, fin(d["T"])]]})
        for prof in res.profiles:
            where = "Height (m)" if prof.kind == "gravity" else "Radius (m)"
            c = chart(("pos", "P", prof.kind), f"Pressure along the {prof.kind} field ({res.equilibrium})",
                      where, "Pressure (bar)")
            c["series"].append({"name": "Pressure", "role": "liquid",
                                "pts": [[fin(a), fin(b / 1e5)] for a, b in zip(prof.position, prof.P)]})
            c2 = chart(("pos", "x", prof.kind), f"Composition along the {prof.kind} field ({res.equilibrium})",
                       where, f"Mole fraction of {prof.names[0]}")
            c2["series"].append({"name": "Local phase", "role": "vapor",
                                 "pts": [[fin(a), fin(b)] for a, b in zip(prof.position, prof.x[:, 0])]})
    return list(charts.values())


def example_files():
    """The example input files, with the checker's verdict on each."""
    folder = os.path.join(os.path.dirname(os.path.abspath(__file__)), "examples")
    out = []
    if not os.path.isdir(folder):
        return {"ok": True, "examples": []}
    for fn in sorted(os.listdir(folder)):
        if not fn.endswith(".txt"):
            continue
        text = open(os.path.join(folder, fn)).read()
        title = next((ln.lstrip("# ").strip() for ln in text.splitlines() if ln.startswith("#")), fn)
        out.append({"file": fn, "title": title, "text": text, "status": check_input(text).status})
    return {"ok": True, "examples": out}


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
        elif op == "examples":
            res = example_files()
        elif op == "pure":
            res = pure(msg["request"])
        elif op == "mixture":
            res = mixture(msg["request"])
        elif op == "mixture_examples":
            res = {"ok": True, "examples": MIXTURE_EXAMPLES}
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
    report("catalog", len(cat["compounds"]) >= 20 and len(cat["models"]) == 6,
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

    r = compute({"components": ["ethanol", "water"], "model": "nrtl", "hold": "T", "T": 350,
                 "params": {"dg12": 100.0, "dg21": 100.0, "alpha": -1.0}})
    report("invalid parameters give a readable error, not a crash",
           not r["ok"] and r.get("error"), r.get("error", ""))

    # field: the 500 m propane / n-butane column forms a gas cap near 306 m
    r = compute({"components": ["propane", "n-butane"], "model": "PR", "hold": "field", "T": 320,
                 "P": 25, "z1": 0.4, "npts": 41,
                 "field": {"kind": "gravity", "g": 9.80665, "height": 500}})
    report("gravity column in the calculator: gas-liquid contact", r["ok"] and len(r["contacts"]) == 1
           and abs(r["contacts"][0]["position"] - 306.1) < 1.5 and r["columns"]["phase"][-1] == "vapor",
           str(r.get("contacts")))
    r = compute({"components": ["ethanol", "water"], "model": "nrtl", "hold": "field", "T": 350,
                 "P": 1.0, "z1": 0.3, "npts": 21,
                 "field": {"kind": "centrifugal", "rpm": 20000, "r_in": 0.05, "r_out": 0.10}})
    report("centrifuge with an activity model: Rackett liquid volumes, P rises outward",
           r["ok"] and "Rackett" in r["input_text"] and r["columns"]["P"][-1] > r["columns"]["P"][0] > 0.99,
           f"{r['columns']['P'][-1] if r.get('ok') else r.get('error')} bar at the rim")

    # "is there equilibrium here?": a state inside the two-phase region, and one above it
    r = point({"components": ["ethanol", "water"], "model": "nrtl", "T": 360.0, "P": 1.01325, "z1": 0.4})
    r2 = point({"components": ["ethanol", "water"], "model": "nrtl", "T": 360.0, "P": 2.0, "z1": 0.4})
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

    # E-PPR78 k_ij(T): the CO2 / ethane azeotrope at 263.15 K (bubble-pressure maximum
    # near 30.8 bar at 30-35 % ethane with k12 = 0.13, CL 202 mixtures slides)
    kc = compound_db.PPR78Kij("carbon-dioxide", "ethane")
    r = compute({"components": ["carbon-dioxide", "ethane"], "model": "PR", "hold": "T",
                 "T": 263.15, "npts": 41})
    az = r.get("azeotrope") or {}
    report("E-PPR78: CO2-ethane k12 = 0.13 and its azeotrope at 263.15 K",
           r["ok"] and abs(kc(263.15) - 0.13) < 0.01 and 0.62 < az.get("x1", 0) < 0.72
           and abs(az.get("P", 0) - 30.8) < 0.6 and "PPR78" in r["input_text"],
           f"k12 = {kc(263.15):.4f}, azeotrope {az}")
    report("E-PPR78 k_ij varies with temperature", abs(kc(230.0) - kc(320.0)) > 1e-3,
           f"{kc(230.0):.4f} at 230 K, {kc(320.0):.4f} at 320 K")
    r = compute({"components": ["carbon-dioxide", "ethane"], "model": "PR", "hold": "T",
                 "T": 263.15, "npts": 41, "params": {"kij": 0.0}})
    report("with k12 = 0 the azeotrope disappears (the slide's point)",
           r["ok"] and r["azeotrope"] is None, str(r.get("azeotrope")))

    # pure fluids
    r = pure({"compound": "propane", "model": "PR", "T": 300.0, "P": 10.0,
              "field": {"kind": "gravity", "g": 9.81, "height": 20.0}})
    i300 = None
    report("pure propane: PR vapor pressure at 300 K near the Antoine value",
           r["ok"] and abs(r["state"]["Psat"] / 9.98 - 1) < 0.02, f"{r.get('state', {}).get('Psat')} bar")
    report("pure propane under gravity: an interface, liquid below and vapor above",
           r["ok"] and r["field"]["column"]["interface"] is not None
           and r["field"]["column"]["phase"][0] == "liquid" and r["field"]["column"]["phase"][-1] == "vapor",
           f"interface at {r['field']['column']['interface']} m")
    b = r["field"]["bands"]
    k = len(b[-1]["T"]) // 2
    report("the band of reference pressures widens with the field",
           b[0]["P_last"][k] - b[0]["P_first"][k] < b[-1]["P_last"][k] - b[-1]["P_first"][k],
           f"{b[0]['P_first'][k]}-{b[0]['P_last'][k]} vs {b[-1]['P_first'][k]}-{b[-1]['P_last'][k]} bar")
    r = pure({"compound": "nitrogen", "model": "SRK", "T": 300.0, "P": 10.0,
              "field": {"kind": "centrifugal", "rpm": 30000, "r_in": 0.01, "r_out": 0.1}})
    report("nitrogen above Tc in a centrifuge: one phase, pressure rising outward",
           r["ok"] and "interface_curve" not in r["field"] and r["field"]["column"]["P"][-1] > r["field"]["column"]["P"][0],
           f"{r['field']['column']['P'][0]} to {r['field']['column']['P'][-1]} bar")

    # mixtures of three and four components
    for ex in MIXTURE_EXAMPLES:
        m = mixture(dict(ex, map=len(ex["components"]) == 3))
        ok = m["ok"] and m["phase"] == "two-phase"
        if ok:
            z_, x_, y_, b_ = np.array(m["z"]), np.array(m["x"]), np.array(m["y"]), m["beta"]
            ok = (np.max(np.abs((1 - b_) * x_ + b_ * y_ - z_)) < 1e-5
                  and abs(sum(x_) - 1) < 1e-5 and abs(sum(y_) - 1) < 1e-5
                  and np.max(np.abs(np.array(m["K"]) - y_ / x_)) < 1e-3 * max(m["K"])
                  and m["bounds"]["T_bubble"] < ex["T"] < m["bounds"]["T_dew"])
            if "map" in m:
                key = m["map"]["kind"]
                ok = ok and sum(p[key] is not None for p in m["map"]["points"]) >= 30
        report(f"flash, {ex['label']}: mass balance, K = y/x, between bubble and dew T", ok,
               f"beta = {m.get('beta')}, {m.get('ms')} ms" if m.get("ok") else m.get("error"))
    m = mixture({"components": ["acetone", "methanol", "n-heptane"], "model": "nrtl",
                 "T": 340, "P": 1.0, "z": [1, 1, 1]})
    report("a mixture with a pair missing from the NRTL database is refused",
           not m["ok"] and "no NRTL" in m["error"], m.get("error", ""))

    rt = json.loads(handle(json.dumps({"op": "compute", "request": {
        "components": ["methanol", "water"], "model": "nrtl", "hold": "P", "P": 1.0}})))
    report("JSON round trip", rt["ok"] and "NaN" not in json.dumps(rt))

    print("\nALL TESTS PASSED" if ok_all else "\nSOME TESTS FAILED")
