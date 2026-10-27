"""
input_checker.py  --  Stage 1 of the Phase Equilibrium Calculator

Decides whether the information a user gives is sufficient to calculate
phase equilibrium, which equilibria can be calculated, and what is missing.

FREE-FORM INPUT
The user lists whatever they know, one quantity per line:

    components = benzene, naphthalene
    P = 1 atm
    sweep = x
    benzene.Tm = 278.68 K
    benzene.dHfus = 9.87 kJ/mol
    activity = vanlaar
    A(benzene, naphthalene) = 0.2
    kij(propane, n-butane) = 0.003
    x(L) = 0.3                 (n or n-1 mole fractions, in component order)

Every quantity is optional. The checker works out what can be calculated
from whatever is given, so the same input format serves any mixture.

Run `python input_checker.py` for demos, or
    `python input_checker.py my_input.txt` to check your own file.
Type `python input_checker.py --vocab` to list every recognized quantity.

FRAMEWORK
Calculations use the Helmholtz form of the fundamental equation,
    dA = -S dT - P dV + sum_i mu_i dN_i              (+ Phi dm with a field)
Phi = g z (gravity) or -omega^2 r^2 / 2 (centrifugal); m = sum_i M_i N_i.
  * Cubic EOS (vdW, PR, SRK) are pressure-explicit, so they give A directly:
        A_res(T, V, N) = -integral_inf^V (P - N R T / V) dV
    and mu_i, fugacity coefficients, densities follow by differentiation.
  * Vapor pressures, latent heats and melting points give the chemical
    potential of a pure component in each phase (Gibbs-Helmholtz); an
    optional activity model (Margules, van Laar, NRTL) corrects the liquid
    for non-ideal mixing. This route needs no critical constants, so it also
    covers components that never vaporize.
  * Each phase has n + 2 + m independent variables; equality of T, P and
    mu_i (or mu_i + M_i Phi with a field) gives
        F = n + 2 + m - pi            (open specification)
        2 + m                         (closed system, feed z known: Duhem)

This version handles two-phase equilibria: vapor-liquid (VLE) and
solid-liquid with a pure solid (SLE). External fields: gravity, centrifugal.
"""

import ast
import difflib
import math
import re
import sys
from dataclasses import dataclass, field
from itertools import combinations
from typing import Optional

R = 8.314462618          # J/(mol K)
COMP_TOL = 1e-6
G_STANDARD = 9.80665     # m/s2


# ==========================================================================
# 1. Units  (everything is converted to SI internally; molar mass -> kg/mol)
# ==========================================================================
# unit key (lower-case, no spaces) -> (factor, offset):  SI = value*factor + offset
UNITS = {
    "temperature": {"k": (1, 0), "c": (1, 273.15), "degc": (1, 273.15),
                    "°c": (1, 273.15)},
    "pressure": {"pa": (1, 0), "kpa": (1e3, 0), "mpa": (1e6, 0),
                 "bar": (1e5, 0), "atm": (101325, 0), "mmhg": (133.322368, 0),
                 "torr": (133.322368, 0), "psi": (6894.757, 0)},
    "energy_mol": {"j/mol": (1, 0), "kj/mol": (1e3, 0), "cal/mol": (4.184, 0),
                   "kcal/mol": (4184, 0)},
    "heat_cap": {"j/mol/k": (1, 0), "j/(molk)": (1, 0), "j/mol-k": (1, 0),
                 "cal/mol/k": (4.184, 0)},
    "molar_mass": {"g/mol": (1e-3, 0), "kg/mol": (1, 0), "kg/kmol": (1e-3, 0)},
    "molar_volume": {"m3/mol": (1, 0), "cm3/mol": (1e-6, 0), "l/mol": (1e-3, 0),
                     "ml/mol": (1e-6, 0), "m3/kmol": (1e-3, 0)},
    "a_vdw": {"pam6/mol2": (1, 0), "barl2/mol2": (0.1, 0),
              "atml2/mol2": (0.101325, 0)},
    "length": {"m": (1, 0), "cm": (1e-2, 0), "mm": (1e-3, 0), "km": (1e3, 0)},
    "accel": {"m/s2": (1, 0), "m/s^2": (1, 0)},
    "angvel": {"rad/s": (1, 0), "rpm": (2 * math.pi / 60, 0)},
    "none": {"": (1, 0), "-": (1, 0)},
}
# Unit assumed when the user gives none (a warning is always printed)
DEFAULT_UNIT = {"temperature": "K", "pressure": "Pa", "energy_mol": "J/mol",
                "heat_cap": "J/mol/K", "molar_mass": "g/mol",
                "molar_volume": "cm3/mol", "a_vdw": "Pa m6/mol2",
                "length": "m", "accel": "m/s2", "angvel": "rad/s", "none": ""}


def _unit_key(u):
    return u.lower().replace(" ", "").replace("*", "").replace("·", "")


def to_si(value, unit, dim):
    """Convert a number in `unit` to SI. Raises ValueError on unknown unit."""
    key = _unit_key(unit)
    table = UNITS[dim]
    if key not in table:
        raise ValueError(f"unknown unit '{unit}' for {dim.replace('_', ' ')} "
                         f"(allowed: {', '.join(sorted(k for k in table if k))})")
    f, off = table[key]
    return value * f + off


def from_si(value, unit, dim):
    f, off = UNITS[dim][_unit_key(unit)]
    return (value - off) / f


# ==========================================================================
# 2. Vocabulary of recognized quantities
# ==========================================================================
# (scope, name) -> (kind, dimension, description)
#   scope: system | component | pair
#   kind : scalar | list | names | word | bool | composition | antoine | relation
VOCAB = {
    # ---- system definition ----
    ("system", "components"): ("names", None, "component names, comma separated"),
    ("system", "model"): ("word", None, "cubic EOS: vdW, PR or SRK (optional)"),
    ("system", "activity"): ("word", None, "liquid activity model: ideal (default), margules, vanlaar or nrtl"),
    ("system", "phases"): ("word", None, "restrict to one equilibrium: 'L, V' or 'L, S(name)'"),
    ("system", "sweep"): ("word", None, "variable swept to draw a diagram: x, T or P"),
    # ---- state ----
    ("system", "T"): ("scalar", "temperature", "temperature"),
    ("system", "P"): ("scalar", "pressure", "pressure"),
    ("system", "x"): ("composition", "none", "phase composition: x(L) liquid, x(V) vapor"),
    ("system", "y"): ("composition", "none", "vapor composition (same as x(V))"),
    ("system", "z"): ("composition", "none", "overall (feed) composition"),
    ("system", "beta"): ("scalar", "none", "vapor fraction (moles vapor / total)"),
    # ---- external field ----
    ("system", "field"): ("word", None, "external field: gravity or centrifugal"),
    ("system", "g"): ("scalar", "accel", "gravitational acceleration"),
    ("system", "omega"): ("scalar", "angvel", "angular velocity (centrifugal field)"),
    ("system", "height"): ("list", "length", "height range of the system: z_start, z_end"),
    ("system", "radius"): ("list", "length", "radius range of the system: r_start, r_end"),
    ("system", "reference"): ("scalar", "length", "position where the given T, P, x apply"),
    # ---- pure-component data ----
    ("component", "Tc"): ("scalar", "temperature", "critical temperature"),
    ("component", "Pc"): ("scalar", "pressure", "critical pressure"),
    ("component", "omega"): ("scalar", "none", "acentric factor"),
    ("component", "a"): ("scalar", "a_vdw", "van der Waals a (instead of Tc, Pc)"),
    ("component", "b"): ("scalar", "molar_volume", "van der Waals b (instead of Tc, Pc)"),
    ("component", "M"): ("scalar", "molar_mass", "molar mass"),
    ("component", "Tb"): ("scalar", "temperature", "normal boiling point (at 1 atm)"),
    ("component", "dHvap"): ("scalar", "energy_mol", "enthalpy of vaporization"),
    ("component", "Tm"): ("scalar", "temperature", "melting point"),
    ("component", "dHfus"): ("scalar", "energy_mol", "enthalpy of fusion"),
    ("component", "dCp_fus"): ("scalar", "heat_cap", "Cp(liquid) - Cp(solid) (optional)"),
    ("component", "V_L"): ("scalar", "molar_volume", "liquid molar volume"),
    ("component", "V_S"): ("scalar", "molar_volume", "solid molar volume"),
    ("component", "antoine"): ("antoine", None, "A, B, C (P unit, T unit[, ln|log10])"),
    ("component", "antoine_range"): ("list", "temperature", "Tmin, Tmax of Antoine validity"),
    ("component", "Psat"): ("relation", None, "vapor pressure relation in T, e.g. exp(...) [kPa, K]"),
    ("component", "nonvolatile"): ("bool", None, "yes if the component never enters the vapor"),
    # ---- pair data ----
    ("pair", "kij"): ("scalar", "none", "binary interaction parameter for the EOS"),
    ("pair", "A"): ("scalar", "none", "Margules / van Laar: A(i, j) = ln gamma_i at infinite dilution in j"),
    ("pair", "dg"): ("scalar", "energy_mol", "NRTL: dg(i, j) = g_ij - g_jj, so tau_ij = dg(i, j)/RT"),
    ("pair", "alpha"): ("scalar", "none", "NRTL non-randomness alpha(i, j) (default 0.3)"),
}

MODEL_ALIASES = {"vdw": "vdW", "vanderwaals": "vdW", "van der waals": "vdW",
                 "pr": "PR", "pengrobinson": "PR", "peng-robinson": "PR",
                 "srk": "SRK", "soave": "SRK", "soaveredlichkwong": "SRK",
                 "soave-redlich-kwong": "SRK"}
ACTIVITY_ALIASES = {"ideal": "ideal", "none": "ideal", "margules": "margules",
                    "vanlaar": "vanlaar", "van laar": "vanlaar", "nrtl": "nrtl"}
ORDERED_PAIR = {"A", "dg"}           # value depends on the order (i, j)
PAIR_QTYS = {"kij", "A", "dg", "alpha"}
FIELD_ALIASES = {"gravity": "gravity", "gravitational": "gravity",
                 "centrifugal": "centrifugal", "centrifuge": "centrifugal"}


# ==========================================================================
# 3. Safe parsing of user relations (e.g. Psat as a function of T)
# ==========================================================================
_SAFE_FUNCS = {"exp": math.exp, "log": math.log, "ln": math.log,
               "log10": math.log10, "sqrt": math.sqrt}
_SAFE_NODES = (ast.Expression, ast.BinOp, ast.UnaryOp, ast.Constant, ast.Name,
               ast.Call, ast.Load, ast.Add, ast.Sub, ast.Mult, ast.Div,
               ast.Pow, ast.USub, ast.UAdd)


def compile_relation(expr, variables):
    """Turn a text expression into a Python function, allowing only arithmetic,
    the listed variables and exp/log/ln/log10/sqrt. Raises ValueError."""
    expr = expr.replace("^", "**")
    try:
        tree = ast.parse(expr, mode="eval")
    except SyntaxError as e:
        raise ValueError(f"cannot read expression '{expr}': {e.msg}")
    for node in ast.walk(tree):
        if not isinstance(node, _SAFE_NODES):
            raise ValueError(f"'{type(node).__name__}' is not allowed in a relation")
        if isinstance(node, ast.Name) and node.id not in variables \
                and node.id not in _SAFE_FUNCS:
            raise ValueError(f"unknown symbol '{node.id}' (allowed variables: "
                             f"{', '.join(variables)}; functions: "
                             f"{', '.join(_SAFE_FUNCS)})")
        if isinstance(node, ast.Call) and (not isinstance(node.func, ast.Name)
                                           or node.func.id not in _SAFE_FUNCS
                                           or node.keywords):
            raise ValueError("only exp, log, ln, log10 and sqrt may be called")
    code = compile(tree, "<relation>", "eval")
    return lambda **kw: eval(code, {"__builtins__": {}}, {**_SAFE_FUNCS, **kw})


# ==========================================================================
# 4. Parsed input
# ==========================================================================
@dataclass
class Inputs:
    components: list = field(default_factory=list)
    model: Optional[str] = None
    phases: Optional[tuple] = None
    sweep: Optional[str] = None
    T: Optional[float] = None
    P: Optional[float] = None
    beta: Optional[float] = None
    comp: dict = field(default_factory=dict)        # 'L'/'V'/'z' -> list
    field_kind: Optional[str] = None
    g: Optional[float] = None
    omega: Optional[float] = None
    height: Optional[list] = None
    radius: Optional[list] = None
    reference: Optional[float] = None
    props: dict = field(default_factory=dict)       # comp -> {name: value}
    pairs: dict = field(default_factory=dict)       # frozenset -> {name: value}
    opairs: dict = field(default_factory=dict)      # (i, j) -> {name: value}, ordered
    activity: str = "ideal"


LINE_RE = re.compile(
    r"^\s*([A-Za-z][\w\-]*)(?:\.([A-Za-z]\w*))?\s*(?:\(([^)]*)\))?\s*=\s*(.*?)\s*$")


def _split_number_unit(text):
    """'0, 100 m' -> ([0.0, 100.0], 'm');  '1 atm' -> ([1.0], 'atm')"""
    tokens = text.replace(",", " ").split()
    nums, i = [], 0
    while i < len(tokens):
        try:
            nums.append(float(tokens[i]))
            i += 1
        except ValueError:
            break
    return nums, " ".join(tokens[i:])


def parse_input(text):
    """Parse free-form text. Returns (Inputs, errors, warnings)."""
    inp, errors, warnings = Inputs(), [], []
    pending = []    # property lines, checked once the component list is known

    def convert(qty, raw, dim, lineno):
        nums, unit = _split_number_unit(raw)
        if not nums:
            raise ValueError("no number found")
        if unit == "" and dim != "none":
            unit = DEFAULT_UNIT[dim]
            warnings.append(f"line {lineno}: no unit for {qty}; assumed {unit}")
        return [to_si(v, unit, dim) for v in nums]

    for lineno, line in enumerate(text.splitlines(), 1):
        line = line.split("#", 1)[0].strip()
        if not line:
            continue
        m = LINE_RE.match(line)
        if not m:
            errors.append(f"line {lineno}: cannot read '{line}' "
                          f"(expected 'quantity = value [unit]')")
            continue
        first, second, arg, raw = m.groups()
        if second:                          # component.property
            scope, qty, target = "component", second, first
        elif arg and first in PAIR_QTYS:
            scope, qty, target = "pair", first, arg
        else:
            scope, qty, target = "system", first, arg
        key = (scope, qty)
        if key not in VOCAB:
            names = [n for s, n in VOCAB if s == scope]
            close = difflib.get_close_matches(qty, names, n=1)
            hint = f" Did you mean '{close[0]}'?" if close else ""
            errors.append(f"line {lineno}: unknown {scope} quantity '{qty}'.{hint}")
            continue
        if scope != "system":
            pending.append((lineno, scope, qty, target, raw))
            continue

        kind, dim, _ = VOCAB[key]
        try:
            if qty == "components":
                names = [s.strip() for s in raw.split(",") if s.strip()]
                if len({n.lower() for n in names}) != len(names):
                    raise ValueError("a component is listed twice")
                inp.components = names
            elif qty == "model":
                k = raw.lower().replace(" ", "")
                if k not in MODEL_ALIASES and raw.lower() not in MODEL_ALIASES:
                    raise ValueError(f"unknown model '{raw}' (use vdW, PR or SRK)")
                inp.model = MODEL_ALIASES.get(k) or MODEL_ALIASES[raw.lower()]
            elif qty == "activity":
                k = raw.lower().strip()
                if k not in ACTIVITY_ALIASES:
                    raise ValueError(f"unknown activity model '{raw}' "
                                     f"(use ideal, margules, vanlaar or nrtl)")
                inp.activity = ACTIVITY_ALIASES[k]
            elif qty == "field":
                k = raw.lower().strip()
                if k in ("electric", "magnetic"):
                    raise ValueError(f"{k} fields are not supported in this version")
                if k not in FIELD_ALIASES:
                    raise ValueError(f"unknown field '{raw}' (use gravity or centrifugal)")
                inp.field_kind = FIELD_ALIASES[k]
            elif qty == "sweep":
                k = raw.strip()
                if k not in ("x", "T", "P"):
                    raise ValueError("sweep must be x, T or P")
                inp.sweep = k
            elif qty == "phases":
                inp.phases = tuple(p.strip() for p in raw.split(",") if p.strip())
            elif kind == "composition":
                phase = "z" if qty == "z" else "V" if qty == "y" else (arg or "").strip()
                if phase not in ("L", "V", "z"):
                    raise ValueError("write x(L) for liquid or x(V) for vapor")
                vals, unit = _split_number_unit(raw)
                if unit:
                    raise ValueError(f"unexpected text '{unit}' in mole fractions")
                inp.comp[phase] = vals
            elif kind == "list":
                vals = convert(qty, raw, dim, lineno)
                if len(vals) != 2:
                    raise ValueError("give two values: start, end")
                setattr(inp, qty, vals)
            else:                                   # scalar
                vals = convert(qty, raw, dim, lineno)
                if len(vals) != 1:
                    raise ValueError("expected a single value")
                setattr(inp, qty, vals[0])
        except ValueError as e:
            errors.append(f"line {lineno}: {qty}: {e}")

    # ---- component and pair properties ----
    known = {c.lower(): c for c in inp.components}
    if not inp.components:
        errors.append("no 'components = ...' line given")
    for c in inp.components:
        inp.props[c] = {}

    for lineno, scope, qty, target, raw in pending:
        kind, dim, _ = VOCAB[(scope, qty)]
        try:
            if scope == "pair":
                names = [s.strip() for s in target.split(",")]
                if len(names) != 2 or names[0].lower() == names[1].lower():
                    raise ValueError(f"write {qty}(component1, component2)")
                unknown = [n for n in names if n.lower() not in known]
                if unknown:
                    raise ValueError(f"'{unknown[0]}' is not in the component list")
                vals = convert(qty, raw, dim, lineno)
                if len(vals) != 1:
                    raise ValueError("expected a single value")
                if qty in ORDERED_PAIR:
                    key2 = (known[names[0].lower()], known[names[1].lower()])
                    inp.opairs.setdefault(key2, {})[qty] = vals[0]
                else:
                    pair = frozenset(known[n.lower()] for n in names)
                    inp.pairs.setdefault(pair, {})[qty] = vals[0]
                continue

            if target.lower() not in known:
                raise ValueError(f"'{target}' is not in the component list")
            c = known[target.lower()]
            if kind == "bool":
                v = raw.lower()
                if v not in ("yes", "no", "true", "false", "1", "0"):
                    raise ValueError("use yes or no")
                inp.props[c][qty] = v in ("yes", "true", "1")
            elif kind == "antoine":
                m = re.match(r"^(.*?)\(([^)]*)\)\s*$", raw)
                body, units = (m.group(1), m.group(2)) if m else (raw, "")
                nums, rest = _split_number_unit(body)
                if len(nums) != 3 or rest:
                    raise ValueError("give A, B, C followed by (P unit, T unit)")
                parts = [u.strip() for u in units.split(",") if u.strip()]
                if len(parts) < 2:
                    parts = ["mmHg", "C"]
                    warnings.append(f"line {lineno}: no units for {c} Antoine "
                                    f"constants; assumed log10(P/mmHg) = A - B/(T/C + C)")
                base = parts[2].lower() if len(parts) > 2 else "log10"
                if base not in ("log10", "ln"):
                    raise ValueError("third unit entry must be log10 or ln")
                to_si(1, parts[0], "pressure"), to_si(1, parts[1], "temperature")
                inp.props[c][qty] = dict(A=nums[0], B=nums[1], C=nums[2],
                                         P_unit=parts[0], T_unit=parts[1], base=base)
            elif kind == "relation":
                m = re.match(r"^(.*?)\[([^\]]*)\]\s*$", raw)
                body, units = (m.group(1), m.group(2)) if m else (raw, "")
                parts = [u.strip() for u in units.split(",") if u.strip()]
                if len(parts) != 2:
                    parts = ["Pa", "K"]
                    warnings.append(f"line {lineno}: no units for {c} Psat relation; "
                                    f"assumed result in Pa with T in K "
                                    f"(write e.g. '... [kPa, K]')")
                to_si(1, parts[0], "pressure"), to_si(1, parts[1], "temperature")
                fn = compile_relation(body.strip(), ["T"])
                inp.props[c][qty] = dict(expr=body.strip(), fn=fn,
                                         P_unit=parts[0], T_unit=parts[1])
            elif kind == "list":
                vals = convert(f"{c}.{qty}", raw, dim, lineno)
                if len(vals) != 2:
                    raise ValueError("give two values: min, max")
                inp.props[c][qty] = vals
            else:
                vals = convert(f"{c}.{qty}", raw, dim, lineno)
                if len(vals) != 1:
                    raise ValueError("expected a single value")
                inp.props[c][qty] = vals[0]
        except ValueError as e:
            errors.append(f"line {lineno}: {target}.{qty}: {e}")

    return inp, errors, warnings


# ==========================================================================
# 5. Vapor pressure from whatever the user gave (used for checks here,
#    and later by the solver)
# ==========================================================================
def psat_source(p):
    """Name of the best available vapor-pressure source, or None."""
    if "antoine" in p:
        return "Antoine equation"
    if "Psat" in p:
        return "user Psat(T) relation"
    if "dHvap" in p and "Tb" in p:
        return "Clausius-Clapeyron (dHvap, Tb)"
    return None


def psat(p, T):
    """Vapor pressure in Pa at T (K) from the first available source."""
    if "antoine" in p:
        a = p["antoine"]
        t = from_si(T, a["T_unit"], "temperature")
        arg = a["A"] - a["B"] / (t + a["C"])
        val = 10 ** arg if a["base"] == "log10" else math.exp(arg)
        return to_si(val, a["P_unit"], "pressure")
    if "Psat" in p:
        r = p["Psat"]
        val = r["fn"](T=from_si(T, r["T_unit"], "temperature"))
        return to_si(val, r["P_unit"], "pressure")
    if "dHvap" in p and "Tb" in p:
        return 101325 * math.exp(-p["dHvap"] / R * (1 / T - 1 / p["Tb"]))
    return None


# ==========================================================================
# 6. Report
# ==========================================================================
@dataclass
class EquilibriumResult:
    label: str
    phases: tuple
    status: str = "UNKNOWN"
    route: Optional[str] = None
    route_kind: Optional[str] = None     # 'eos' | 'vp' | 'thermal' (for the next module)
    criteria: Optional[str] = None
    dof_line: Optional[str] = None
    calc: Optional[str] = None
    missing: list = field(default_factory=list)
    errors: list = field(default_factory=list)
    notes: list = field(default_factory=list)
    options: list = field(default_factory=list)


@dataclass
class Report:
    status: str = "UNKNOWN"
    summary: str = ""
    fundamental: list = field(default_factory=list)
    equilibria: list = field(default_factory=list)
    errors: list = field(default_factory=list)
    warnings: list = field(default_factory=list)
    inputs: Optional[Inputs] = None      # parsed, SI-converted input (for the next module)

    def __str__(self):
        out = [f"OVERALL: {self.status}" + (f" -- {self.summary}" if self.summary else "")]
        if self.errors:
            out.append("Errors:")
            out += [f"  - {e}" for e in self.errors]
        if self.fundamental:
            out.append("Fundamental relation:")
            out += [f"  {s}" for s in self.fundamental]
        for k, eq in enumerate(self.equilibria, 1):
            out.append(f"[{k}] {eq.label}  ->  {eq.status}")
            for title, val in (("Route", eq.route), ("Criteria", eq.criteria),
                               ("Phase rule", eq.dof_line), ("Calculation", eq.calc)):
                if val:
                    out.append(f"    {title}: {val}")
            for title, items in (("Errors", eq.errors), ("Missing", eq.missing),
                                 ("You could specify", eq.options),
                                 ("Notes", eq.notes)):
                if items:
                    out.append(f"    {title}:")
                    out += [f"      - {s}" for s in items]
        if self.warnings:
            out.append("Input warnings:")
            out += [f"  - {w}" for w in self.warnings]
        return "\n".join(out)


# ==========================================================================
# 7. Checks
# ==========================================================================
EOS_NEEDS = {"vdW": "Tc and Pc (or a and b)", "PR": "Tc, Pc and omega",
             "SRK": "Tc, Pc and omega"}


def eos_route(inp):
    """Can the chosen cubic EOS describe every component in both phases?
    Returns (ok, missing, notes); ok is None when no model was chosen."""
    if inp.model is None:
        return None, [], []
    missing, notes = [], []
    for c in inp.components:
        p = inp.props[c]
        if inp.model == "vdW":
            if not (("Tc" in p and "Pc" in p) or ("a" in p and "b" in p)):
                missing.append(f"{c}: Tc and Pc (or vdW a and b)")
        else:
            absent = [q for q in ("Tc", "Pc", "omega") if q not in p]
            if absent:
                missing.append(f"{c}: {', '.join(absent)}")
    for ci, cj in combinations(inp.components, 2):
        if "kij" not in inp.pairs.get(frozenset((ci, cj)), {}):
            notes.append(f"kij({ci}, {cj}) not given; 0 will be used")
    if not missing:
        notes.insert(0, f"{inp.model} with van der Waals one-fluid mixing rules")
    return (not missing), missing, notes


ACTIVITY_LABEL = {"ideal": "ideal solution (gamma = 1)",
                  "margules": "two-parameter Margules activity model",
                  "vanlaar": "van Laar activity model",
                  "nrtl": "NRTL activity model"}


def activity_check(inp):
    """Parameters the chosen liquid activity model still needs.
    Returns (missing, notes)."""
    model, missing, notes = inp.activity, [], []
    if model == "ideal":
        return missing, notes
    for ci, cj in combinations(inp.components, 2):
        if model in ("margules", "vanlaar"):
            for a, b in ((ci, cj), (cj, ci)):
                if "A" not in inp.opairs.get((a, b), {}):
                    missing.append(f"A({a}, {b}) for the {model} model "
                                   f"(ln gamma of {a} at infinite dilution in {b})")
        else:
            for a, b in ((ci, cj), (cj, ci)):
                if "dg" not in inp.opairs.get((a, b), {}):
                    missing.append(f"dg({a}, {b}) for the NRTL model (energy, e.g. J/mol)")
            if "alpha" not in inp.pairs.get(frozenset((ci, cj)), {}):
                notes.append(f"alpha({ci}, {cj}) not given; 0.3 will be used")
    return missing, notes


def vapor_pressure_route(inp):
    """Pure-component vapor pressures + liquid activity model + ideal gas.
    Returns (ok, missing, notes, volatile_components)."""
    volatile = [c for c in inp.components if not inp.props[c].get("nonvolatile")]
    missing, notes = [], []
    if not volatile:
        missing.append("at least one volatile component (all are marked nonvolatile)")
    for c in volatile:
        src = psat_source(inp.props[c])
        if src is None:
            missing.append(f"{c}: vapor pressure -- give antoine, a Psat(T) relation, "
                           f"or dHvap and Tb; or mark {c}.nonvolatile = yes")
        elif src.startswith("Clausius"):
            notes.append(f"{c}: Psat from Clausius-Clapeyron with constant dHvap "
                         f"(less accurate far from Tb)")
    nonvol = [c for c in inp.components if c not in volatile]
    if nonvol and volatile:
        notes.append(f"nonvolatile: {', '.join(nonvol)} -- absent from the vapor")
    a_miss, a_notes = activity_check(inp)
    missing += a_miss
    notes += a_notes
    if not missing:
        lead = f"liquid: {ACTIVITY_LABEL[inp.activity]}; vapor: ideal gas"
        if inp.activity == "ideal":
            lead += "; no azeotropes can be predicted"
        notes.insert(0, lead)
    return (not missing), missing, notes, volatile


def describe_fundamental(inp):
    names = inp.components
    mu = " + ".join(f"mu_{c} dN_{c}" for c in names)
    out = [f"dA = -S dT - P dV + {mu}"]
    if inp.field_kind == "gravity":
        out.append("     + Phi dm,  Phi = g z,  m = sum M_i N_i   [gravity]")
    elif inp.field_kind == "centrifugal":
        out.append("     + Phi dm,  Phi = -omega^2 r^2 / 2,  m = sum M_i N_i   [centrifugal]")
    return out


def _criteria(inp, phases, present_in_vapor):
    t = ("T uniform; P hydrostatic (varies with position); "
         if inp.field_kind else "T, P uniform; ")
    pot = "mu_i + M_i Phi" if inp.field_kind else "mu_i"
    if phases[1] == "V":
        if len(present_in_vapor) == len(inp.components):
            return t + f"{pot} equal in L and V for every component"
        return t + f"{pot} equal in L and V for {', '.join(present_in_vapor)}"
    solid = phases[1][2:-1]
    return t + f"{pot} of {solid} equal in the liquid and the pure solid"


def _calc_name(inp, is_sle, solid, has):
    T, P, xL, xV, z, b = (has[k] for k in ("T", "P", "L", "V", "z", "beta"))
    sw = inp.sweep
    if is_sle:
        if sw == "x" and P:
            return f"T-x solid-liquid curve for solid {solid} at fixed P"
        if sw == "T" and P:
            return f"saturation composition vs T for solid {solid} at fixed P"
        if P and xL:
            return f"temperature at which solid {solid} forms (P, x -> T)"
        if P and T:
            return f"liquid composition in equilibrium with solid {solid} (P, T -> x)"
        return None
    if z:
        if T and P:
            return "isothermal flash (T, P, z -> x, y, beta)"
        if P and b:
            return "flash at fixed P and vapor fraction"
        if T and b:
            return "flash at fixed T and vapor fraction"
        return None
    if sw == "x":
        if T:
            return "Pxy diagram at fixed T (bubble and dew curves)"
        if P:
            return "Txy diagram at fixed P (bubble and dew curves)"
    if sw == "T" and P:
        return "Txy diagram at fixed P (bubble and dew curves)"
    if sw == "P" and T:
        return "Pxy diagram at fixed T (bubble and dew curves)"
    if xL and T:
        return "bubble-point pressure (T, x -> P, y)"
    if xL and P:
        return "bubble-point temperature (P, x -> T, y)"
    if xV and T:
        return "dew-point pressure (T, y -> P, x)"
    if xV and P:
        return "dew-point temperature (P, y -> T, x)"
    return None


def evaluate(inp, phases, global_errors):
    n = len(inp.components)
    is_sle = phases[1].startswith("S(")
    solid = phases[1][2:-1] if is_sle else None
    label = (f"Solid-liquid: liquid + pure solid {solid}" if is_sle
             else "Vapor-liquid: liquid + vapor")
    res = EquilibriumResult(label=label, phases=phases)

    # ---------- how the chemical potentials will be calculated ----------
    route_kind = None
    vapor_comps = list(inp.components)
    if is_sle:
        p = inp.props[solid]
        absent = [f"{solid}.{q} ({d})" for q, d in
                  (("Tm", "melting point"), ("dHfus", "enthalpy of fusion")) if q not in p]
        res.missing += absent
        a_miss, a_notes = activity_check(inp)
        res.missing += a_miss
        res.notes += a_notes
        res.route = (f"solid {solid}: dHfus and Tm (Gibbs-Helmholtz); "
                     f"liquid: {ACTIVITY_LABEL[inp.activity]}")
        route_kind = "thermal"
        if "dCp_fus" not in p:
            res.notes.append(f"{solid}.dCp_fus not given; heat-capacity term neglected")
        if inp.model:
            res.notes.append(f"the EOS is not used for solid-liquid equilibrium; the "
                             f"liquid uses the {ACTIVITY_LABEL[inp.activity]}")
    else:
        e_ok, e_miss, e_notes = eos_route(inp)
        v_ok, v_miss, v_notes, volatile = vapor_pressure_route(inp)
        if e_ok:
            route_kind, res.route = "eos", f"cubic EOS ({inp.model}) for both phases"
            res.notes += e_notes
            if inp.activity != "ideal":
                res.notes.append(f"the {inp.activity} activity model is not used for "
                                 f"vapor-liquid equilibrium on the EOS route")
            if any(inp.props[c].get("nonvolatile") for c in inp.components):
                res.notes.append("'nonvolatile' is ignored on the EOS route; the EOS "
                                 "decides each component's volatility")
        elif v_ok:
            route_kind = "vp"
            res.route = ("vapor pressures + ideal solution (Raoult's law)"
                         if inp.activity == "ideal" else
                         f"vapor pressures + {ACTIVITY_LABEL[inp.activity]} "
                         f"(modified Raoult's law)")
            vapor_comps = volatile
            res.notes += v_notes
            if inp.model:
                res.notes.append(f"{inp.model} route incomplete ({'; '.join(e_miss)}); "
                                 f"vapor-pressure route used instead")
        else:
            vapor_comps = volatile
            if inp.model:
                res.missing.append(f"EITHER for the {inp.model} route: " + "; ".join(e_miss))
                res.missing.append("OR for the vapor-pressure route: " + "; ".join(v_miss))
            else:
                res.missing += v_miss
                res.missing.append("(alternatively choose model = PR, SRK or vdW and "
                                   "give its constants for each component)")
        # vapor-pressure sanity: T above Tc on the vapor-pressure route
        if route_kind == "vp" and inp.T is not None:
            for c in vapor_comps:
                if "Tc" in inp.props[c] and inp.T > inp.props[c]["Tc"]:
                    res.notes.append(f"T is above Tc of {c}; Raoult's law is not "
                                     f"valid -- use an EOS")
        res.notes += _antoine_range_notes(inp, vapor_comps)
    res.criteria = _criteria(inp, phases, vapor_comps)
    res.route_kind = route_kind

    # ---------- external field ----------
    m = 1 if inp.field_kind else 0
    field_given = 0
    if inp.field_kind:
        for c in inp.components:
            if "M" not in inp.props[c]:
                res.missing.append(f"{c}.M (molar mass; couples {c} to the field)")
        if route_kind != "eos":
            for c in inp.components:
                if "V_L" not in inp.props[c]:
                    res.missing.append(f"{c}.V_L (liquid molar volume; needed for the "
                                       f"pressure dependence of mu in the field)")
            if is_sle and "V_S" not in inp.props[solid]:
                res.missing.append(f"{solid}.V_S (solid molar volume)")
        if inp.field_kind == "gravity":
            if inp.g is None:
                res.notes.append(f"g not given; standard value {G_STANDARD} m/s2 used")
            field_given = 1
            if inp.height is None:
                res.missing.append("height = z_start, z_end (extent of the system)")
            ext = inp.height
        else:
            if inp.omega is None:
                res.missing.append("omega (angular velocity, rad/s or rpm)")
            else:
                field_given = 1
            if inp.radius is None:
                res.missing.append("radius = r_start, r_end (extent of the system)")
            ext = inp.radius
        if ext is not None:
            if inp.reference is None:
                res.notes.append(f"reference position not given; the state is taken "
                                 f"at the start, {ext[0]} m")
            elif not min(ext) <= inp.reference <= max(ext):
                res.errors.append("reference position lies outside the system extent")

    # ---------- degrees of freedom ----------
    has = {"T": inp.T is not None, "P": inp.P is not None,
           "L": "L" in inp.comp, "V": "V" in inp.comp, "z": "z" in inp.comp,
           "beta": inp.beta is not None}
    sweep_n = 0
    if inp.sweep:
        sw = inp.sweep
        if (sw == "x" and (has["L"] or has["V"])) or (sw in ("T", "P") and has[sw]):
            res.errors.append(f"{sw} is both swept and given a value")
        elif sw == "x" and n != 2:
            res.errors.append("sweep = x is only defined for binary mixtures here")
        elif has["z"]:
            res.errors.append("sweep cannot be combined with a feed composition z")
        else:
            sweep_n = 1

    n_vap = len(vapor_comps)
    if has["z"]:
        if is_sle:
            res.errors.append("flash with a feed z is only implemented for VLE")
        dof = 2 + m
        state = has["T"] + has["P"] + has["beta"]
        res.dof_line = (f"closed system, feed known (Duhem): 2 + m = {dof}; "
                        f"specified {state + field_given}")
        if has["L"] or has["V"]:
            res.notes.append("phase compositions given with z; they are normally "
                             "results of a flash, not inputs")
    else:
        dof = n + 2 + m - 2
        state = has["T"] + has["P"] + sweep_n + (n - 1) * has["L"]
        if has["V"]:
            if is_sle:
                res.notes.append("vapor composition ignored for solid-liquid equilibrium")
            elif n_vap > 1:
                state += n_vap - 1
            else:
                res.notes.append(f"vapor is pure {vapor_comps[0]}; its composition "
                                 f"is not an independent variable")
        if has["beta"]:
            res.notes.append("vapor fraction is ignored without a feed composition z")
        res.dof_line = (f"F = n + 2 + m - pi = {n} + 2 + {m} - 2 = {dof}; "
                        f"specified {state + field_given}"
                        + (" (incl. swept variable)" if sweep_n else ""))

    state_dof = dof - m
    if res.errors:
        res.status = "INVALID"
    elif state < state_dof:
        res.status = "INSUFFICIENT"
        res.missing.insert(0, f"{state_dof - state} more state variable(s)")
        if is_sle:
            res.options = ["P and x(L) -> temperature at which the solid forms",
                           "P and T -> saturated liquid composition",
                           "P and sweep = x -> full solid-liquid curve"]
        elif has["z"]:
            res.options = ["any two of T, P, beta (with z) -> flash"]
        else:
            res.options = ["T and x(L) -> bubble P", "P and x(L) -> bubble T",
                           "T and x(V) -> dew P", "P and x(V) -> dew T",
                           "T (or P) and sweep = x -> Pxy (Txy) diagram",
                           "T, P and z -> flash"]
    elif state > state_dof:
        res.status = "OVERSPECIFIED"
        res.errors.append(f"{state} state variables given but only {state_dof} can be "
                          f"fixed; remove {state - state_dof} (the rest are results)")
    else:
        res.calc = _calc_name(inp, is_sle, solid, has)
        if res.calc is None:
            res.notes.append("variable count is right but this combination is not a "
                             "standard calculation")
        if inp.field_kind:
            res.calc = (res.calc or "equilibrium") + (
                f", with profile along {'height' if inp.field_kind == 'gravity' else 'radius'}")
        res.status = "INSUFFICIENT" if res.missing else "SUFFICIENT"
    return res


def _antoine_range_notes(inp, comps):
    notes = []
    if inp.T is None:
        return notes
    for c in comps:
        p = inp.props[c]
        if "antoine" in p and "antoine_range" in p:
            lo, hi = p["antoine_range"]
            if not lo <= inp.T <= hi:
                notes.append(f"T outside the Antoine range of {c} "
                             f"({lo:.1f}-{hi:.1f} K); Psat is extrapolated")
    return notes


def _global_checks(inp, errors, warnings):
    n = len(inp.components)
    if 0 < n < 2:
        errors.append("at least two components are needed for a mixture")
    if inp.activity in ("margules", "vanlaar") and n > 2:
        errors.append(f"the {inp.activity} model is for binary mixtures; use "
                      f"activity = nrtl for {n} components")
    if inp.activity == "ideal" and (inp.opairs or any(
            q in v for v in inp.pairs.values() for q in ("alpha",))):
        warnings.append("activity-model parameters given but activity = ideal; "
                        "they are ignored (set activity = margules, vanlaar or nrtl)")
    if inp.T is not None and inp.T <= 0:
        errors.append("temperature must be above 0 K")
    if inp.P is not None and inp.P <= 0:
        errors.append("pressure must be positive")
    if inp.beta is not None and not 0 <= inp.beta <= 1:
        errors.append("beta must lie between 0 and 1")
    if inp.field_kind is None and any(v is not None for v in
                                      (inp.g, inp.omega, inp.height, inp.radius,
                                       inp.reference)):
        warnings.append("field data given but no 'field = ...' line; ignored")
    if inp.field_kind == "gravity" and inp.radius is not None:
        warnings.append("radius is given but the field is gravity; ignored")
    if inp.field_kind == "centrifugal" and inp.height is not None:
        warnings.append("height is given but the field is centrifugal; ignored")
    if inp.radius is not None and min(inp.radius) < 0:
        errors.append("radius values must not be negative")

    # compositions: allow n-1 values, complete the last one
    for ph, vals in list(inp.comp.items()):
        name = {"L": "x(L)", "V": "x(V)", "z": "z"}[ph]
        if n < 2:
            break
        if len(vals) == n - 1:
            last = 1 - sum(vals)
            if last < -COMP_TOL:
                errors.append(f"{name}: fractions add up to more than 1")
                continue
            vals = vals + [max(last, 0.0)]
        elif len(vals) != n:
            errors.append(f"{name}: expected {n} (or {n - 1}) values, got {len(vals)}")
            continue
        if any(v < 0 or v > 1 for v in vals):
            errors.append(f"{name}: mole fractions must lie between 0 and 1")
        elif abs(sum(vals) - 1) > COMP_TOL:
            errors.append(f"{name}: mole fractions add up to {sum(vals):.6f}, not 1")
        inp.comp[ph] = vals

    # property sanity
    for c, p in inp.props.items():
        for q in ("Tc", "Pc", "M", "Tb", "Tm", "dHvap", "dHfus", "V_L", "V_S", "a", "b"):
            if q in p and p[q] <= 0:
                errors.append(f"{c}.{q} must be positive")
        if "omega" in p and not -0.5 <= p["omega"] <= 1.5:
            warnings.append(f"{c}.omega = {p['omega']} is unusual for an acentric factor")
        if "Tc" in p and "Tb" in p and p["Tb"] >= p["Tc"]:
            errors.append(f"{c}: Tb must be below Tc")
        if "Psat" in p:
            T_test = p.get("Tb", 298.15)
            try:
                val = psat(p, T_test)
                if not (val > 0 and math.isfinite(val)):
                    raise ValueError
                if "Tb" in p and not 0.5 < val / 101325 < 2:
                    warnings.append(f"{c}: Psat relation gives {val / 1e3:.1f} kPa at "
                                    f"Tb; expected about 101.3 kPa -- check its units")
            except (ValueError, ZeroDivisionError, OverflowError):
                errors.append(f"{c}: Psat relation does not give a positive finite "
                              f"value at T = {T_test:.2f} K")


def check_input(text):
    """Main entry point: free-form text -> Report."""
    inp, errors, warnings = parse_input(text)
    rep = Report(errors=errors, warnings=warnings, inputs=inp)
    _global_checks(inp, rep.errors, rep.warnings)
    if rep.errors or len(inp.components) < 2:
        rep.status = "INVALID"
        rep.summary = "fix the errors below first"
        return rep

    rep.fundamental = describe_fundamental(inp)

    # which equilibria to examine
    if inp.phases:
        ph = inp.phases
        valid = (len(ph) == 2 and ph[0] == "L" and
                 (ph[1] == "V" or re.fullmatch(r"S\(.+\)", ph[1])))
        if not valid:
            rep.status, rep.summary = "INVALID", "unsupported phases"
            rep.errors.append("phases must be 'L, V' or 'L, S(component)' in this version"
                              " (liquid-liquid and three-phase come later)")
            return rep
        if ph[1].startswith("S("):
            s = ph[1][2:-1]
            match = [c for c in inp.components if c.lower() == s.lower()]
            if not match:
                rep.status, rep.summary = "INVALID", "unknown solid"
                rep.errors.append(f"S({s}): '{s}' is not in the component list")
                return rep
            ph = ("L", f"S({match[0]})")
        candidates = [ph]
    else:
        candidates = [("L", "V")]
        for c in inp.components:
            p = inp.props[c]
            if "Tm" in p or "dHfus" in p:
                candidates.append(("L", f"S({c})"))

    rep.equilibria = [evaluate(inp, ph, rep.errors) for ph in candidates]

    ok = [e for e in rep.equilibria if e.status == "SUFFICIENT"]
    if ok:
        rep.status = "SUFFICIENT"
        rep.summary = f"{len(ok)} of {len(rep.equilibria)} equilibria can be calculated"
    elif any(e.status == "OVERSPECIFIED" for e in rep.equilibria):
        rep.status, rep.summary = "OVERSPECIFIED", "too many variables fixed"
    elif all(e.status == "INVALID" for e in rep.equilibria):
        rep.status, rep.summary = "INVALID", "see errors"
    else:
        rep.status, rep.summary = "INSUFFICIENT", "see what is missing"
    if not inp.phases and len(candidates) == 1:
        rep.warnings.append("only VLE examined; give Tm and dHfus of a component to "
                            "also check solid-liquid equilibrium")
    return rep


def print_vocabulary():
    for scope in ("system", "component", "pair"):
        print(f"\n{scope.upper()} quantities"
              + {"system": "", "component": "  (write  name.quantity = ...)",
                 "pair": "  (write  kij(name1, name2) = ...; A and dg depend on the order)"}[scope])
        for (s, q), (kind, dim, desc) in VOCAB.items():
            if s == scope:
                unit = f"[default {DEFAULT_UNIT[dim]}]" if dim and dim != "none" else ""
                print(f"  {q:14s} {desc} {unit}")


# ==========================================================================
# 8. Demo cases  (property values are typical literature values, for
#    demonstration only -- check them against a data source before real use)
# ==========================================================================
DEMOS = {
    "1. Propane / n-butane, PR, Pxy diagram at 300 K": """
components = propane, n-butane
model = PR
T = 300 K
sweep = x
propane.Tc = 369.8 K
propane.Pc = 42.48 bar
propane.omega = 0.152
n-butane.Tc = 425.1 K
n-butane.Pc = 37.96 bar
n-butane.omega = 0.200
kij(propane, n-butane) = 0.003
""",
    "2. Same system, SRK, acentric factor and state missing": """
components = propane, n-butane
model = SRK
propane.Tc = 369.8 K
propane.Pc = 42.48 bar
propane.omega = 0.152
n-butane.Tc = 425.1 K
n-butane.Pc = 37.96 bar
P = 5 bar
""",
    "3. Benzene / naphthalene at 1 atm from latent heats only": """
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
""",
    "4. Propane / n-butane flash in a 500 m gravity column, one molar mass missing": """
components = propane, n-butane
model = PR
T = 320 K
P = 12 bar
z = 0.4
field = gravity
height = 0, 500 m
propane.Tc = 369.8 K
propane.Pc = 42.48 bar
propane.omega = 0.152
propane.M = 44.10 g/mol
n-butane.Tc = 425.1 K
n-butane.Pc = 37.96 bar
n-butane.omega = 0.200
""",
    "5. Ethanol / water, user Psat relations, centrifuge": """
components = ethanol, water
T = 350 K
sweep = x
field = centrifugal
omega = 20000 rpm
radius = 0.05, 0.10 m
reference = 0.05 m
ethanol.Psat = exp(16.8958 - 3795.17/(T - 273.15 + 230.918)) [kPa, K]
water.Psat   = exp(16.3872 - 3885.70/(T - 273.15 + 230.170)) [kPa, K]
ethanol.M = 46.07 g/mol
water.M = 18.015 g/mol
ethanol.V_L = 58.7 cm3/mol
water.V_L = 18.07 cm3/mol
""",
    "8. Ethanol / water at 1 atm, van Laar activity model": """
components = ethanol, water
P = 1 atm
sweep = x
activity = vanlaar
A(ethanol, water) = 1.6798
A(water, ethanol) = 0.9227
ethanol.Psat = exp(16.8958 - 3795.17/(T - 273.15 + 230.918)) [kPa, K]
water.Psat   = exp(16.3872 - 3885.70/(T - 273.15 + 230.170)) [kPa, K]
""",
    "6. Input mistakes: typo and bad composition": """
components = methane, ethane
model = vdW
T = 200 K
x(V) = 0.9, 0.2
methane.Tcc = 190.6 K
""",
    "7. Overspecified: T, P, x and y all given (vdW with a, b for ethane)": """
components = methane, ethane
model = vdW
T = 200 K
P = 20 bar
x(L) = 0.3
x(V) = 0.8
methane.Tc = 190.6 K
methane.Pc = 45.99 bar
ethane.a = 5.562 bar L2/mol2
ethane.b = 0.0638 L/mol
""",
}

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--vocab":
        print_vocabulary()
    elif len(sys.argv) > 1:
        with open(sys.argv[1]) as f:
            print(check_input(f.read()))
    else:
        for title, text in DEMOS.items():
            print("=" * 78)
            print(title)
            print("-" * 78)
            print(check_input(text))
