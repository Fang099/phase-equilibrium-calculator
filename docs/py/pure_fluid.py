"""Phase equilibrium of a pure fluid from a cubic equation of state, with and
without an external field.

A pure fluid has one component, so the phase rule leaves F = 1 - 2 + 2 = 1
for two phases: at each temperature below Tc there is one pressure, the
vapor pressure Psat(T), at which liquid and vapor coexist. The condition is
equal fugacity,

    f_L(T, Psat) = f_V(T, Psat),   with  ln f = ln P + ln phi(T, P),

and ln phi comes from the residual Helmholtz energy of the equation of state
(thermo_models.CubicEOS.ln_phi_TV), the same derivative used for mixtures.

In a gravitational or centrifugal field the equilibrium condition becomes

    mu + M Phi = constant  along the fluid, at uniform T,

where Phi is the potential energy per unit mass (g h for gravity,
-omega^2 r^2 / 2 for a centrifuge) and M the molar mass. Since
d mu = RT d ln f at constant T, this gives the fugacity at every position:

    ln f(pos) = ln f(ref) - M [Phi(pos) - Phi(ref)] / RT.

At each position the pressure is the one at which the fluid has that
fugacity (liquid if f > f_sat, vapor if f < f_sat), so the column holds
liquid where the potential is low and vapor where it is high, and the
interface sits where f = f_sat. Seen from a reference point, the single
vapor-pressure curve P = Psat(T) widens into a band: between "first liquid
appears somewhere" and "the last vapor disappears", whose width grows with
the field strength.

Functions
---------
make_eos(compound, model)               the cubic EOS of one compound from the database
spinodal(eos, T)                        pressures where the isotherm's loop turns (T < Tc)
saturation_pressure(eos, T)             Psat and the liquid and vapor molar volumes
saturation_temperature(eos, P)          Tsat at a pressure
saturation_curve(eos, ...)              Psat, densities, enthalpy of vaporization over T
isotherm(eos, T, ...)                   P against V, with the equal-area (Maxwell) line
pressure_at_fugacity(eos, T, ln_f)      inverse of ln f(T, P) on the stable branch
column(eos, M, T, P_ref, field, ...)    P, density and phase along a field at fixed T
band(eos, M, field, T_values)           reference pressures of first liquid / last vapor
"""
import math
from dataclasses import dataclass, field as dcfield

import numpy as np
try:
    from scipy.optimize import brentq
except ImportError:                     # in the browser: the project's own root finder
    from numerics import brentq

from thermo_models import CubicEOS, R

ONE = np.array([1.0])


# --------------------------------------------------------------------------
# the equation of state of one compound
# --------------------------------------------------------------------------
def make_eos(compound=None, model="PR", Tc=None, Pc=None, omega=None):
    """CubicEOS for one compound: from the database by name, or from Tc (K),
    Pc (Pa) and omega. Returns (eos, M in kg/mol or None, record)."""
    rec = None
    M = None
    if compound is not None:
        import compound_db
        rec = compound_db.find(compound)
        if rec is None:
            raise ValueError(f"'{compound}' is not in the database")
        Tc, Pc, omega = rec["Tc"]["value"], rec["Pc"]["value"], rec["omega"]["value"]
        M = rec["M"]["value"] / 1000.0
    eos = CubicEOS(model, Tc=[Tc], Pc=[Pc], omega=[omega])
    return eos, M, rec


def _ln_f(eos, T, P, phase):
    """ln(f / Pa) of the pure fluid on one branch ('liquid' or 'vapor')."""
    return math.log(P) + float(eos.ln_phi(T, P, ONE, phase)[0])


def _Z(eos, T, P, phase):
    v = eos.molar_volume(T, P, ONE, phase)
    return P * v / (R * T), v


# --------------------------------------------------------------------------
# saturation
# --------------------------------------------------------------------------
def spinodal(eos, T):
    """The pressures at the two turning points of the isotherm P(V) below Tc
    (the liquid spinodal, often negative, and the vapor spinodal). Between
    them the cubic has three roots. Returns (P_low, P_high) or None above
    the critical temperature of the equation of state."""
    b, a = eos.b[0], eos.a_pure(T)[0]
    d1, d2 = eos.d1, eos.d2
    Pv = lambda V: R * T / (V - b) - a / ((V + d1 * b) * (V + d2 * b))
    V = np.geomspace(b * (1 + 1e-6), b * 1e4, 600)
    P = Pv(V)
    dP = np.diff(P)
    turns = np.where(np.sign(dP[:-1]) != np.sign(dP[1:]))[0]
    if len(turns) < 2:
        return None

    def dPdV(v):
        h = v * 1e-7
        return Pv(v + h) - Pv(v - h)

    ext = []
    for i in turns[:2]:
        try:
            v = brentq(dPdV, V[i], V[i + 2], xtol=1e-14 * V[i])
        except ValueError:
            v = V[i + 1]
        ext.append(Pv(v))
    return min(ext), max(ext)


def saturation_pressure(eos, T):
    """Vapor pressure (Pa) of the pure fluid at T (K), with the liquid and
    vapor molar volumes (m3/mol): the pressure between the spinodals at which
    ln f_L = ln f_V. Raises ValueError at or above the critical point."""
    sp = spinodal(eos, T)
    if sp is None:
        raise ValueError(f"T = {T:.2f} K is not below the critical temperature")
    hi = sp[1] * (1 - 1e-9)
    # lower end: the liquid spinodal, or (when that is negative, at low T)
    # a pressure well below the Lee-Kesler/Wilson estimate of Psat; at
    # vanishing pressure the cubic cannot be resolved into three roots
    Tc, Pc, w = eos.Tc[0], eos.Pc[0], eos.omega[0]
    P_est = Pc * 10 ** (7 / 3 * (1 + w) * (1 - Tc / T))
    lo = max(sp[0] * (1 + 1e-9), min(P_est, hi) * 1e-3)
    if not hi > lo:
        raise ValueError("no two-phase region at this temperature")

    def g(lnP):
        P = math.exp(lnP)
        return _ln_f(eos, T, P, "liquid") - _ln_f(eos, T, P, "vapor")

    # g > 0 below Psat (the vapor is the stable phase, lower fugacity) and
    # g < 0 above it
    a_, b_ = math.log(lo), math.log(hi)
    while g(a_) <= 0 and a_ > math.log(max(sp[0], 1e-30) + 1e-30) - 30:
        a_ -= 2.0
    lnP = brentq(g, a_, b_, xtol=1e-13)
    P = math.exp(lnP)
    vL = eos.molar_volume(T, P, ONE, "liquid")
    vV = eos.molar_volume(T, P, ONE, "vapor")
    return P, vL, vV


def saturation_temperature(eos, P, T_lo=None):
    """Temperature (K) at which the vapor pressure equals P (Pa)."""
    Tc = eos.Tc[0]
    T_lo = T_lo or 0.3 * Tc
    T_hi = Tc * (1 - 1e-3)
    f = lambda T: math.log(saturation_pressure(eos, T)[0]) - math.log(P)
    if f(T_hi) < 0:
        raise ValueError("P is above the critical pressure of the equation of state")
    if f(T_lo) > 0:
        raise ValueError("P is below the vapor pressure at the lowest temperature tried")
    return brentq(f, T_lo, T_hi, xtol=1e-9)


@dataclass
class SatCurve:
    T: np.ndarray                 # K
    P: np.ndarray                 # Pa
    vL: np.ndarray                # m3/mol
    vV: np.ndarray
    dHvap: np.ndarray             # J/mol, from Clapeyron: T (vV - vL) dPsat/dT
    Tc: float
    Pc: float                     # of the equation of state (its own critical point)
    vc: float


def eos_critical_point(eos):
    """The critical point of the cubic itself: Tc and Pc are the inputs, and
    the critical volume is Zc RT/P (vdW 3/8, SRK 1/3, PR 0.3074)."""
    Tc, Pc = eos.Tc[0], eos.Pc[0]
    Zc = {"vdW": 0.375, "SRK": 1 / 3, "PR": 0.30740}[eos.model]
    return Tc, Pc, Zc * R * Tc / Pc


def saturation_curve(eos, T_min=None, npts=60):
    """Vapor pressure and the coexisting densities from T_min (default
    0.5 Tc) up to the critical point, points crowded near Tc."""
    Tc, Pc, vc = eos_critical_point(eos)
    T_min = T_min or 0.5 * Tc
    u = np.linspace(0, 1, npts)
    Ts = Tc - (Tc - T_min) * (1 - u) ** 2          # finer near Tc
    Ts = Ts[Ts < Tc * (1 - 2e-4)]
    rows = []
    for T in Ts:
        try:
            rows.append((T,) + saturation_pressure(eos, T))
        except (ValueError, RuntimeError):
            pass
    T, P, vL, vV = (np.array(c) for c in zip(*rows))
    lnP = np.log(P)
    dlnP = np.gradient(lnP, T)
    dH = T * (vV - vL) * P * dlnP                   # Clapeyron
    T = np.append(T, Tc); P = np.append(P, Pc)
    vL = np.append(vL, vc); vV = np.append(vV, vc); dH = np.append(dH, 0.0)
    return SatCurve(T, P, vL, vV, dH, Tc, Pc, vc)


def isotherm(eos, T, npts=300, V_max_factor=None):
    """P against V on one isotherm (one mole), and below Tc the equal-area
    (Maxwell) line at Psat between vL and vV. Returns a dict."""
    b = eos.b[0]
    Tc, Pc, vc = eos_critical_point(eos)
    V_hi = (V_max_factor or 30) * vc
    V = np.geomspace(b * 1.05, V_hi, npts)
    P = np.array([eos.pressure(T, v, ONE) for v in V])
    out = {"T": T, "V": V, "P": P, "Psat": None, "vL": None, "vV": None}
    if T < Tc * (1 - 1e-4):
        try:
            Ps, vL, vV = saturation_pressure(eos, T)
            out.update(Psat=Ps, vL=vL, vV=vV)
            # equal areas: the integral of (P - Psat) dV from vL to vV is zero
            Vs = np.linspace(vL, vV, 4001)
            dP = np.array([eos.pressure(T, v, ONE) - Ps for v in Vs])
            area = float(np.sum(0.5 * (dP[1:] + dP[:-1]) * np.diff(Vs)))
            out["maxwell_area"] = area / (Ps * (vV - vL))
        except ValueError:
            pass
    return out


# --------------------------------------------------------------------------
# fugacity -> pressure, and the fluid in a field
# --------------------------------------------------------------------------
def pressure_at_fugacity(eos, T, ln_f, sat=None):
    """The pressure (Pa) at which the pure fluid at T has fugacity exp(ln_f),
    on its stable branch, and that branch: liquid if f > f_sat, vapor if
    f < f_sat; above Tc there is a single branch ('fluid')."""
    Tc = eos.Tc[0]
    if T < Tc * (1 - 1e-4):
        Ps = sat[0] if sat else saturation_pressure(eos, T)[0]
        ln_fs = _ln_f(eos, T, Ps, "vapor")
        if abs(ln_f - ln_fs) < 1e-11:
            return Ps, "saturated"
        if ln_f >= ln_fs:
            phase, lo, hi = "liquid", math.log(Ps), math.log(Ps) + 1.0
            while _ln_f(eos, T, math.exp(hi), "liquid") < ln_f:
                lo, hi = hi, hi + 2.0
                if hi > math.log(1e12):
                    raise ValueError("fugacity out of range")
        else:
            phase, lo, hi = "vapor", ln_f - 2.0, math.log(Ps)
            while _ln_f(eos, T, math.exp(lo), "vapor") > ln_f:
                lo -= 3.0
        g = lambda lnP: _ln_f(eos, T, math.exp(lnP), phase) - ln_f
        return math.exp(brentq(g, lo, hi, xtol=1e-12)), phase
    phase = "fluid"
    lo, hi = ln_f - 3.0, ln_f + 3.0
    br = "vapor"
    g = lambda lnP: _ln_f(eos, T, math.exp(lnP), br) - ln_f
    while g(hi) < 0:
        hi += 2.0
    while g(lo) > 0:
        lo -= 2.0
    return math.exp(brentq(g, lo, hi, xtol=1e-12)), phase


@dataclass
class PureField:
    """A gravitational or centrifugal field over a length of fluid."""
    kind: str                     # 'gravity' | 'centrifugal'
    strength: float               # g in m/s2, or omega in rad/s
    start: float                  # m (height from the bottom, or radius)
    end: float                    # m
    reference: float              # m, where the given pressure applies

    def phi(self, pos):
        """Potential energy per unit mass, J/kg (zero at pos = 0)."""
        if self.kind == "gravity":
            return self.strength * pos
        return -0.5 * self.strength ** 2 * pos ** 2

    def position_at(self, phi):
        """The position (within the field) with potential phi, or None."""
        if self.kind == "gravity":
            pos = phi / self.strength if self.strength else None
        else:
            if self.strength == 0 or phi > 0:
                return None
            pos = math.sqrt(-2 * phi) / self.strength
        if pos is None or not (min(self.start, self.end) - 1e-12 <= pos <= max(self.start, self.end) + 1e-12):
            return None
        return pos


@dataclass
class Column:
    position: np.ndarray          # m
    P: np.ndarray                 # Pa
    ln_f: np.ndarray              # ln(f / Pa)
    rho: np.ndarray               # kg/m3
    phase: list
    interface: object             # position of the liquid-vapor interface (m) or None
    Psat: object                  # Pa, at T (None above Tc)
    notes: list = dcfield(default_factory=list)


def column(eos, M, T, P_ref, fld, npts=101):
    """The pure fluid at uniform T in a field: from the pressure P_ref at the
    reference position, the fugacity everywhere,
        ln f(pos) = ln f(ref) - M [Phi(pos) - Phi(ref)] / RT,
    then the pressure, phase and density at each position."""
    Tc = eos.Tc[0]
    sat = None
    if T < Tc * (1 - 1e-4):
        sat = saturation_pressure(eos, T)
    ref_phase = "fluid" if sat is None else ("liquid" if P_ref > sat[0] else "vapor")
    ln_f_ref = _ln_f(eos, T, P_ref, "vapor" if ref_phase != "liquid" else "liquid")
    pos = np.linspace(fld.start, fld.end, npts)
    P, lnf, rho, ph = [], [], [], []
    for x in pos:
        lf = ln_f_ref - M * (fld.phi(x) - fld.phi(fld.reference)) / (R * T)
        p, phase = pressure_at_fugacity(eos, T, lf, sat)
        v = eos.molar_volume(T, p, ONE, "liquid" if phase == "liquid" else "vapor")
        P.append(p); lnf.append(lf); rho.append(M / v); ph.append(phase)
    interface, notes = None, []
    if sat is not None:
        ln_fs = _ln_f(eos, T, sat[0], "vapor")
        phi_star = fld.phi(fld.reference) + R * T * (ln_f_ref - ln_fs) / M
        interface = fld.position_at(phi_star)
        if interface is None:
            notes.append(f"the whole {'column' if fld.kind == 'gravity' else 'rotor'} is "
                         f"{ph[0]}: the fugacity never reaches f_sat inside it")
    else:
        notes.append("above the critical temperature: one fluid phase, its density "
                     "changing continuously along the field")
    return Column(pos, np.array(P), np.array(lnf), np.array(rho), ph, interface,
                  sat[0] if sat else None, notes)


def band(eos, M, fld, T_values):
    """For each T below Tc, the reference pressures between which the column
    holds both phases:
        P_first  liquid first appears, where the potential is lowest
        P_last   the last vapor disappears, where the potential is highest
    With no field both equal Psat(T). Returns arrays (T, Psat, P_first, P_last)."""
    phis = [fld.phi(fld.start), fld.phi(fld.end)]
    if fld.kind == "centrifugal" and fld.start < 0 < fld.end:
        phis.append(0.0)
    phi_min, phi_max, phi_ref = min(phis), max(phis), fld.phi(fld.reference)
    out = []
    for T in T_values:
        try:
            sat = saturation_pressure(eos, T)
        except ValueError:
            continue
        ln_fs = _ln_f(eos, T, sat[0], "vapor")
        # first liquid: f at phi_min reaches f_sat; last vapor: f at phi_max does
        lf_first = ln_fs + M * (phi_min - phi_ref) / (R * T)
        lf_last = ln_fs + M * (phi_max - phi_ref) / (R * T)
        p1, _ = pressure_at_fugacity(eos, T, lf_first, sat)
        p2, _ = pressure_at_fugacity(eos, T, lf_last, sat)
        out.append((T, sat[0], p1, p2))
    T, Ps, P1, P2 = (np.array(c) for c in zip(*out))
    return T, Ps, P1, P2


# --------------------------------------------------------------------------
# an input file with one component (after the input checker)
# --------------------------------------------------------------------------
def solve_input(rep):
    """Solve a pure-fluid input that the checker found SUFFICIENT. Returns
    (text, charts): the printed results, and chart data in the format of
    web_api.run_input."""
    from input_checker import psat as psat_relation
    inp = rep.inputs
    c = inp.components[0]
    p = inp.props[c]
    T, P = inp.T, inp.P
    lines, charts = [], []
    fin = lambda v: None if v is None or not np.isfinite(v) else float(f"{v:.6g}")
    if inp.model is None:                                  # vapor-pressure relation only
        f = lambda t: psat_relation(p, t)
        if T is not None and P is not None:
            Ps = f(T)
            lines.append(f"At T = {T:.2f} K the vapor pressure is {Ps / 1e5:.6g} bar, so at "
                         f"P = {P / 1e5:.6g} bar the stable phase is {'liquid' if P > Ps else 'vapor'}.")
        elif T is not None:
            lines.append(f"Vapor pressure of {c} at {T:.2f} K: {f(T) / 1e5:.6g} bar")
        elif P is not None:
            lo, hi = 50.0, p.get("Tc", 2000.0) * 0.999
            Tb = brentq(lambda t: math.log(f(t)) - math.log(P), lo, hi)
            lines.append(f"Boiling temperature of {c} at {P / 1e5:.6g} bar: {Tb:.3f} K")
        lo, hi = p.get("antoine_range", [None, None]) if "antoine_range" in p else (None, None)
        Tc = p.get("Tc")
        lo = lo or (0.5 * Tc if Tc else 250.0)
        hi = hi or (0.95 * Tc if Tc else 500.0)
        Ts = np.linspace(lo, hi, 60)
        charts.append({"title": f"Vapor pressure of {c} (vapor-pressure relation)", "xTitle": "Temperature (K)",
                       "yTitle": "Pressure (bar)", "points": [],
                       "series": [{"name": "P_sat(T)", "role": "liquid",
                                   "pts": [[fin(t), fin(f(t) / 1e5)] for t in Ts]}]})
        return "\n".join(lines), charts

    eos = CubicEOS(inp.model, Tc=[p.get("Tc")], Pc=[p.get("Pc")], omega=[p.get("omega", 0.0)],
                   a_vdw=[p.get("a")] if "a" in p else None, b_vdw=[p.get("b")] if "b" in p else None)
    M = p.get("M")
    Tc, Pc, vc = eos_critical_point(eos)
    lines.append(f"{inp.model} critical point of {c}: Tc = {Tc:.2f} K, Pc = {Pc / 1e5:.4g} bar"
                 + (f", rho_c = {M / vc:.4g} kg/m3" if M else ""))
    rho = (lambda v: f"{M / v:.5g} kg/m3") if M else (lambda v: f"{v * 1e6:.5g} cm3/mol")
    if inp.field_kind:
        if inp.field_kind == "gravity":
            fld = PureField("gravity", inp.g, inp.height[0], inp.height[1],
                            inp.reference if inp.reference is not None else inp.height[0])
        else:
            fld = PureField("centrifugal", inp.omega, inp.radius[0], inp.radius[1],
                            inp.reference if inp.reference is not None else inp.radius[0])
        col = column(eos, M, T, P, fld, npts=41)
        lines.append(f"Fluid at {T:.2f} K in a {fld.kind} field, P = {P / 1e5:.6g} bar at {fld.reference:g} m:")
        if col.Psat:
            lines.append(f"  vapor pressure at {T:.2f} K: {col.Psat / 1e5:.6g} bar")
        lines.append("  position (m)   P (bar)      density (kg/m3)   phase")
        for k in range(0, len(col.position), 4):
            lines.append(f"  {col.position[k]:11.4g}   {col.P[k] / 1e5:10.6g}   {col.rho[k]:12.5g}      {col.phase[k]}")
        lines.append(f"  interface: {col.interface:.5g} m" if col.interface is not None else "  " + "; ".join(col.notes))
        where = "Height (m)" if fld.kind == "gravity" else "Radius (m)"
        charts.append({"title": f"Pressure along the {fld.kind} field", "xTitle": where, "yTitle": "Pressure (bar)",
                       "points": [], "series": [{"name": "Pressure", "role": "liquid",
                                                 "pts": [[fin(a), fin(b / 1e5)] for a, b in zip(col.position, col.P)]}]})
        charts.append({"title": f"Density along the {fld.kind} field", "xTitle": where, "yTitle": "Density (kg/m3)",
                       "points": [], "series": [{"name": "Density", "role": "vapor",
                                                 "pts": [[fin(a), fin(b)] for a, b in zip(col.position, col.rho)]}]})
        return "\n".join(lines), charts
    if T is not None and P is not None:
        v = eos.molar_volume(T, P, ONE, "stable")
        msg = f"At T = {T:.2f} K and P = {P / 1e5:.6g} bar: "
        if T < Tc * (1 - 1e-4):
            Ps = saturation_pressure(eos, T)[0]
            msg += f"{'liquid' if P > Ps else 'vapor'} (vapor pressure {Ps / 1e5:.6g} bar)"
        else:
            msg += "above Tc: one fluid phase"
        lines.append(msg + f", density {rho(v)}, Z = {P * v / (R * T):.5g}")
    elif T is not None:
        Ps, vL, vV = saturation_pressure(eos, T)
        lines.append(f"Vapor pressure at {T:.2f} K: {Ps / 1e5:.6g} bar; saturated liquid {rho(vL)}, "
                     f"saturated vapor {rho(vV)}")
    elif P is not None:
        lines.append(f"Boiling temperature at {P / 1e5:.6g} bar: {saturation_temperature(eos, P):.3f} K")
    cur = saturation_curve(eos, npts=50)
    lines.append("Vapor-pressure curve:  T (K)   Psat (bar)")
    for k in range(0, len(cur.T), 7):
        lines.append(f"  {cur.T[k]:9.2f}   {cur.P[k] / 1e5:10.5g}")
    charts.append({"title": f"Vapor pressure of {c} ({inp.model})", "xTitle": "Temperature (K)",
                   "yTitle": "Pressure (bar)", "points": [{"x": fin(Tc), "y": fin(Pc / 1e5), "text": "critical point"}],
                   "series": [{"name": "P_sat(T)", "role": "liquid",
                               "pts": [[fin(a), fin(b / 1e5)] for a, b in zip(cur.T, cur.P)]}]})
    if M:
        charts.append({"title": "Coexisting densities", "xTitle": "Density (kg/m3)", "yTitle": "Temperature (K)",
                       "points": [], "series": [
                           {"name": "Saturated liquid", "role": "liquid", "pts": [[fin(M / a), fin(b)] for a, b in zip(cur.vL, cur.T)]},
                           {"name": "Saturated vapor", "role": "vapor", "pts": [[fin(M / a), fin(b)] for a, b in zip(cur.vV, cur.T)]}]})
    return "\n".join(lines), charts


# --------------------------------------------------------------------------
# tests
# --------------------------------------------------------------------------
if __name__ == "__main__":
    import compound_db
    ok_all = True

    def report(name, ok, detail):
        global ok_all
        ok_all &= bool(ok)
        print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}")

    # 1. van der Waals in reduced units: Psat/Pc at Tr = 0.9 is 0.6470 (equal-area rule)
    eos = CubicEOS("vdW", Tc=[300.0], Pc=[50e5], omega=[0.0])
    P, vL, vV = saturation_pressure(eos, 270.0)
    report("vdW reduced vapor pressure at Tr = 0.9", abs(P / 50e5 - 0.6470) < 2e-3, f"Pr = {P / 50e5:.4f}")
    iso = isotherm(eos, 270.0)
    report("equal areas under the vdW loop", abs(iso["maxwell_area"]) < 1e-3,
           f"(area above - area below)/(Psat (vV - vL)) = {iso['maxwell_area']:.2e}")

    # 2. PR vapor pressures against Antoine (Poling) at Tr = 0.7 and the normal boiling point
    errs = []
    for c in ("water", "carbon-dioxide", "propane", "n-hexane", "benzene", "methanol", "ethanol",
              "nitrogen", "methane"):
        e, M, rec = make_eos(c, "PR")
        a = rec["antoine"]
        if a is None:
            continue
        T = rec["Tb"]["value"]
        Pa = 10 ** (a["A"] - a["B"] / (T + a["C"]))
        Pe = saturation_pressure(e, T)[0]
        errs.append((c, 100 * (Pe / Pa - 1)))
    worst = max(abs(x[1]) for x in errs)
    report("PR Psat at the normal boiling point vs Antoine", worst < 8,
           ", ".join(f"{c} {d:+.1f}%" for c, d in errs))

    # 3. Tsat inverts Psat
    e, M, rec = make_eos("water", "PR")
    Tb = saturation_temperature(e, 101325.0)
    report("PR normal boiling point of water", abs(Tb - 373.15) < 3, f"{Tb:.2f} K")

    # 4. the critical point closes the curve
    cur = saturation_curve(e, npts=40)
    report("liquid and vapor densities meet at the critical point",
           cur.vV[-2] / cur.vL[-2] < 1.5 and abs(cur.P[-1] - rec["Pc"]["value"]) < 1,
           f"vV/vL = {cur.vV[-2] / cur.vL[-2]:.3f} at T = {cur.T[-2]:.2f} K")
    i = np.argmin(abs(cur.T - 373.15))
    report("Clapeyron enthalpy of vaporization of water near 373 K", 37e3 < cur.dHvap[i] < 45e3,
           f"{cur.dHvap[i] / 1000:.1f} kJ/mol at {cur.T[i]:.1f} K (measured 40.7)")

    # 5. a gas column under gravity: ideal-gas limit ln(P/P0) = -M g h / RT
    e, M, rec = make_eos("nitrogen", "PR")
    fld = PureField("gravity", 9.81, 0.0, 1000.0, 0.0)
    col = column(e, M, 300.0, 1e5, fld, npts=11)
    analytic = 1e5 * math.exp(-M * 9.81 * 1000 / (R * 300.0))
    report("nitrogen at 1 bar, 1000 m column: barometric formula",
           abs(col.P[-1] / analytic - 1) < 1e-3, f"{col.P[-1] / 1e5:.5f} bar vs {analytic / 1e5:.5f}")

    # 6. a liquid column: dP/dh = -rho g
    e, M, rec = make_eos("water", "PR")
    col = column(e, M, 300.0, 10e5, PureField("gravity", 9.81, 0.0, 50.0, 0.0), npts=11)
    hyd = 10e5 - col.rho[0] * 9.81 * 50.0
    report("water at 10 bar, 50 m column: hydrostatics", abs(col.P[-1] - hyd) / 10e5 < 2e-3,
           f"{col.P[-1] / 1e5:.4f} bar vs rho g h {hyd / 1e5:.4f} bar")

    # 7. a liquid-vapor interface: propane at 300 K, bottom pressure just above Psat
    e, M, rec = make_eos("propane", "PR")
    Ps = saturation_pressure(e, 300.0)[0]
    col = column(e, M, 300.0, Ps + 0.05e5, PureField("gravity", 9.81, 0.0, 20.0, 0.0), npts=81)
    rhoL = col.rho[0]
    h_est = 0.05e5 / (rhoL * 9.81)
    report("propane interface height ~ (P0 - Psat)/(rho_L g)",
           col.interface is not None and abs(col.interface / h_est - 1) < 0.05,
           f"{col.interface:.3f} m vs {h_est:.3f} m; phases {col.phase[0]} -> {col.phase[-1]}")
    ok_phase = all(p == ("liquid" if x < col.interface else "vapor") for x, p in zip(col.position, col.phase)
                   if abs(x - col.interface) > 0.3)
    report("liquid below the interface, vapor above", ok_phase, "")

    # 8. the band: with no field it collapses to Psat; it widens with the field
    T = np.array([260.0, 300.0, 340.0])
    Tb_, Ps_, P1, P2 = band(e, M, PureField("gravity", 0.0, 0.0, 100.0, 0.0), T)
    report("no field: P_first = P_last = Psat", np.allclose(P1, Ps_, rtol=1e-8) and np.allclose(P2, Ps_, rtol=1e-8), "")
    _, Ps_, P1, P2 = band(e, M, PureField("gravity", 9.81, 0.0, 100.0, 0.0), T)
    _, _, P1b, P2b = band(e, M, PureField("gravity", 98.1, 0.0, 100.0, 0.0), T)
    report("gravity, reference at the bottom: band from Psat upward, wider for larger g",
           np.allclose(P1, Ps_, rtol=1e-8) and np.all(P2 > Ps_) and np.all(P2b - P1b > P2 - P1),
           f"at 300 K: {P1[1] / 1e5:.3f} to {P2[1] / 1e5:.3f} bar (g), to {P2b[1] / 1e5:.3f} bar (10 g)")
    w = 2 * math.pi * 30000 / 60
    _, Ps_, P1, P2 = band(e, M, PureField("centrifugal", w, 0.01, 0.1, 0.01), T)
    report("centrifuge, reference at the inner radius: band below Psat",
           np.all(P1 < Ps_) and np.allclose(P2, Ps_, rtol=1e-8),
           f"at 300 K: {P1[1] / 1e5:.3f} to {P2[1] / 1e5:.3f} bar")

    print("\nALL TESTS PASSED" if ok_all else "\nSOME TESTS FAILED")
