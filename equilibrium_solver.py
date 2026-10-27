"""
equilibrium_solver.py  --  Stage 3 of the Phase Equilibrium Calculator

Solves the equilibrium criteria for the phase models built by
thermo_models.py. Everything rests on one condition from the fundamental
equation: the chemical potential of each component is the same in every
phase, i.e. equal fugacities  ln f_i(phase 1) = ln f_i(phase 2).

VAPOR-LIQUID
  With K_i = y_i / x_i = phi_i(L) / phi_i(V), a material balance on a feed z
  split into a vapor fraction beta gives (Rachford-Rice)
      sum_i z_i (K_i - 1) / (1 + beta (K_i - 1)) = 0
  beta = 0 is a bubble point, beta = 1 a dew point. One routine,
  `saturation`, finds the unknown T or P for any beta; `flash_TP` finds beta
  when T and P are both fixed. K depends on x and y for an EOS, so it is
  updated by successive substitution until it stops changing.

SOLID-LIQUID (pure solid i)
      ln f_i(liquid; T, P, x) = ln f_i(pure solid; T, P)

EXTERNAL FIELD (gravity, centrifugal)
  At equilibrium  ln f_i + M_i Phi / RT  is the same at every position.
  Phases that coexist at one position share the same Phi, so the field
  cancels between them: the LOCAL phase diagram is unchanged. What the field
  does is make pressure and composition vary with position. `column_profile`
  solves for P(pos) and x(pos), and inserts a gas-liquid contact where the
  local fluid reaches saturation; the hydrostatic pressure gradient
  dP/dz = -rho g is not imposed but follows from the condition above
  (via Gibbs-Duhem), which the tests check.

Run `python equilibrium_solver.py` for the verification tests and demos, or
`python equilibrium_solver.py examples/sugar_water.txt` to solve an input file.
"""

import math
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
from scipy.optimize import brentq, fsolve

from thermo_models import (R, P0, EOSPhase, IdealGasPhase, IdealSolutionLiquid,
                           PureSolidPhase, build_system)


class SolverError(Exception):
    pass


# ==========================================================================
# 1. K-values
# ==========================================================================
def K_values(L, V, T, P, x, y):
    """K_i = phi_i(L)/phi_i(V); 0 for components that cannot enter the vapor."""
    K = np.exp(L.ln_phi(T, P, x) - V.ln_phi(T, P, y))
    K[~V.present] = 0.0
    return K


def initial_K(L, V, T, P):
    """Starting K-values: Wilson's correlation for an EOS (needs only Tc, Pc,
    omega), the exact ideal-solution value otherwise."""
    n = len(L.names)
    if isinstance(L, EOSPhase):
        e = L.eos
        K = e.Pc / P * np.exp(5.373 * (1 + e.omega) * (1 - e.Tc / T))
    else:
        x = np.full(n, 1 / n)
        K = np.exp(L.ln_phi(T, P, x) - V.ln_phi(T, P, x))
    K[~V.present] = 0.0
    return K


def _rr_residual(z, K, beta):
    """Rachford-Rice function."""
    return float(np.sum(z * (K - 1) / (1 + beta * (K - 1))))


def _compositions(z, K, beta):
    denom = 1 + beta * (K - 1)
    x = np.divide(z, denom, out=np.zeros_like(z), where=z > 0)
    y = K * x
    if y.sum() <= 0:
        raise SolverError("no volatile component is present")
    return x / x.sum(), y / y.sum()


# ==========================================================================
# 2. Saturation points: bubble (beta = 0), dew (beta = 1), or any fixed beta
# ==========================================================================
@dataclass
class VLEPoint:
    T: float
    P: float
    x: np.ndarray            # liquid
    y: np.ndarray            # vapor
    beta: float              # vapor fraction
    phase: str = "two-phase" # or 'liquid', 'vapor'


def _sat_eval(L, V, T, P, z, beta, Kbox, tol=1e-12, maxit=300):
    """At fixed T and P, iterate K to convergence for the given beta and
    return the Rachford-Rice residual. Kbox holds K between calls (warm start)."""
    K = Kbox[0].copy()
    pres = V.present
    for _ in range(maxit):
        x, y = _compositions(z, K, beta)
        Knew = K_values(L, V, T, P, x, y)
        diff = np.max(np.abs(np.log(Knew[pres]) - np.log(K[pres])))
        K = Knew
        if diff < tol:
            break
    Kbox[0] = K
    x, y = _compositions(z, K, beta)
    # Trivial solution (K = 1 for every volatile component): the EOS found only
    # one root, so both 'phases' are the same fluid. Steer the root finder:
    # a liquid-like fluid means P too high / T too low, and vice versa.
    if pres.sum() > 1 and np.max(np.abs(np.log(K[pres]))) < 1e-6:
        Zc = P * L.molar_volume(T, P, x) / (R * T)
        return (-1.0 if Zc < 0.3 else 1.0), x, y, True
    return _rr_residual(z, K, beta), x, y, False


def _ideal_guess(L, V, z, beta, T=None, P=None):
    """Guess the unknown T or P from K-values that ignore composition."""
    if P is None:                        # T known: K ~ Ki0 * (P_ref / P)
        Pref = 1e5
        K0 = initial_K(L, V, T, Pref)
        f = lambda lnP: _rr_residual(z, K0 * Pref / math.exp(lnP), beta)
        lo, hi = math.log(1e2), math.log(1e9)
        if f(lo) * f(hi) > 0:
            return 1e5
        return math.exp(brentq(f, lo, hi))
    Ts = np.linspace(60, 1500, 300)       # P known: scan T for a sign change
    vals = []
    for Tt in Ts:
        try:
            vals.append(_rr_residual(z, initial_K(L, V, Tt, P), beta))
        except (ValueError, OverflowError, FloatingPointError):
            vals.append(np.nan)
    for k in range(len(Ts) - 1):
        if np.isfinite(vals[k]) and np.isfinite(vals[k + 1]) and vals[k] * vals[k + 1] < 0:
            return float(Ts[k])
    return 300.0


def _bracket(f, u0, step, umin, umax, maxit=60):
    """Step outward from u0 (alternating sides, growing steps) until f changes
    sign. Points where f cannot be evaluated are skipped."""
    def safe(u):
        try:
            v = f(u)
            return v if np.isfinite(v) else None
        except (SolverError, ValueError, OverflowError, ZeroDivisionError,
                FloatingPointError, np.linalg.LinAlgError):
            return None
    f0 = safe(u0)
    left, right = (u0, f0), (u0, f0)
    s = step
    for _ in range(maxit):
        for side in (-1, 1):
            base = left if side < 0 else right
            u = min(max(base[0] + side * s, umin), umax)
            fu = safe(u)
            if fu is None:
                continue
            if base[1] is not None and base[1] * fu <= 0:
                return (u, base[0]) if side < 0 else (base[0], u)
            if side < 0:
                left = (u, fu)
            else:
                right = (u, fu)
            if left[1] is not None and right[1] is not None and left[1] * right[1] <= 0:
                return left[0], right[0]
        s *= 1.6
    raise SolverError("could not bracket a solution")


def saturation(L, V, z, beta=0.0, T=None, P=None, guess=None):
    """Find P (if T given) or T (if P given) at which feed z has vapor
    fraction beta. beta = 0: bubble point; beta = 1: dew point."""
    z = np.asarray(z, float) / np.sum(z)
    if (T is None) == (P is None):
        raise ValueError("give exactly one of T and P")
    if beta >= 1 and np.any((z > 0) & ~V.present):
        raise SolverError("a dew point is impossible: the feed contains a "
                          "nonvolatile component")
    if guess is None:
        guess = _ideal_guess(L, V, z, beta, T=T, P=P)
    Kbox = [initial_K(L, V, T if T is not None else guess,
                      P if P is not None else guess)]

    if P is None:                                  # unknown: ln P
        f = lambda u: _sat_eval(L, V, T, math.exp(u), z, beta, Kbox)[0]
        a, b = _bracket(f, math.log(guess), 0.05, math.log(1.0), math.log(1e9))
        u = brentq(f, a, b, xtol=1e-12, rtol=1e-12)
        Tsol, Psol = T, math.exp(u)
    else:                                          # unknown: T
        f = lambda u: _sat_eval(L, V, u, P, z, beta, Kbox)[0]
        a, b = _bracket(f, guess, 1.0, 20.0, 3000.0)
        u = brentq(f, a, b, xtol=1e-10, rtol=1e-12)
        Tsol, Psol = u, P
    res, x, y, trivial = _sat_eval(L, V, Tsol, Psol, z, beta, Kbox)
    if trivial:
        raise SolverError("only a trivial solution (x = y) was found; the state "
                          "is probably at or beyond the mixture critical point")
    if beta == 0:
        x = z.copy()
    if beta == 1:
        y = z.copy()
    return VLEPoint(Tsol, Psol, x, y, beta)


def bubble_P(L, V, T, x, guess=None):
    return saturation(L, V, x, 0.0, T=T, guess=guess)


def bubble_T(L, V, P, x, guess=None):
    return saturation(L, V, x, 0.0, P=P, guess=guess)


def dew_P(L, V, T, y, guess=None):
    return saturation(L, V, y, 1.0, T=T, guess=guess)


def dew_T(L, V, P, y, guess=None):
    return saturation(L, V, y, 1.0, P=P, guess=guess)


# ==========================================================================
# 3. Isothermal flash
# ==========================================================================
def flash_TP(L, V, T, P, z, tol=1e-12, maxit=500):
    """Split feed z at fixed T and P. Returns VLEPoint (phase 'liquid' or
    'vapor' if only one phase is present)."""
    z = np.asarray(z, float) / np.sum(z)
    K = initial_K(L, V, T, P)
    pres = V.present
    for _ in range(maxit):
        f0, f1 = _rr_residual(z, K, 0.0), _rr_residual(z, K, 1.0 - 1e-12)
        if f0 <= 0:
            beta = 0.0
        elif f1 >= 0:
            beta = 1.0
        else:
            beta = brentq(lambda b: _rr_residual(z, K, b), 0.0, 1.0 - 1e-12,
                          xtol=1e-14)
        if beta in (0.0, 1.0):
            break
        x, y = _compositions(z, K, beta)
        Knew = K_values(L, V, T, P, x, y)
        diff = np.max(np.abs(np.log(Knew[pres]) - np.log(K[pres])))
        K = Knew
        if diff < tol:
            break
    if beta == 0.0:
        return VLEPoint(T, P, z.copy(), np.full_like(z, np.nan), 0.0, "liquid")
    if beta == 1.0:
        return VLEPoint(T, P, np.full_like(z, np.nan), z.copy(), 1.0, "vapor")
    x, y = _compositions(z, K, beta)
    if pres.sum() > 1 and np.max(np.abs(np.log(K[pres]))) < 1e-6:
        Zc = P * L.molar_volume(T, P, z) / (R * T)
        ph = "liquid" if Zc < 0.3 else "vapor"
        return VLEPoint(T, P, z.copy() if ph == "liquid" else np.full_like(z, np.nan),
                        z.copy() if ph == "vapor" else np.full_like(z, np.nan),
                        0.0 if ph == "liquid" else 1.0, ph)
    return VLEPoint(T, P, x, y, beta)


# ==========================================================================
# 4. Binary phase diagrams
# ==========================================================================
@dataclass
class Diagram:
    kind: str                  # 'Pxy' or 'Txy'
    fixed: float               # T (K) or P (Pa)
    x1: np.ndarray             # liquid mole fraction of component 1 (bubble curve)
    y1: np.ndarray             # vapor mole fraction of component 1 (dew curve)
    var: np.ndarray            # P (Pa) or T (K) at each point
    names: list
    notes: list = field(default_factory=list)


def vle_diagram(L, V, T=None, P=None, npts=41, x_min=0.0):
    """Pxy at fixed T, or Txy at fixed P, for a binary. The bubble curve is
    (x1, var) and the dew curve is (y1, var): each bubble-point calculation
    gives one point of each, so no separate dew calculations are needed."""
    names = L.names
    if len(names) != 2:
        raise SolverError("diagrams are for binary mixtures")
    xs = np.linspace(1.0, x_min, npts)          # start from component 1 pure
    out_x, out_y, out_v, notes, fails = [], [], [], [], 0
    guess = None
    for x1 in xs:
        x = np.array([x1, 1 - x1])
        if not np.any((x > 0) & V.present):
            continue                             # nothing volatile in the liquid
        try:
            pt = saturation(L, V, x, 0.0, T=T, P=P, guess=guess)
        except SolverError:
            fails += 1
            continue
        guess = pt.P if P is None else pt.T
        out_x.append(x1)
        out_y.append(pt.y[0])
        out_v.append(pt.P if P is None else pt.T)
    if fails:
        notes.append(f"{fails} point(s) could not be solved (near a mixture "
                     f"critical point, or outside the model's range)")
    if not V.present.all():
        nv = [n for n, p in zip(names, V.present) if not p]
        notes.append(f"{', '.join(nv)} is nonvolatile: the vapor is pure, so the "
                     f"dew curve is the line y1 = {1.0 if V.present[0] else 0.0}")
    order = np.argsort(out_x)
    return Diagram("Pxy" if P is None else "Txy", T if P is None else P,
                   np.array(out_x)[order], np.array(out_y)[order],
                   np.array(out_v)[order], names, notes)


# ==========================================================================
# 5. Solid-liquid equilibrium (pure solid)
# ==========================================================================
def _sle_residual(L, S, T, P, x):
    i = S.index
    return L.ln_f(T, P, x)[i] - S.ln_f(T, P)[i]


def sle_temperature(L, S, P, x, guess=None):
    """Temperature at which pure solid i first forms from liquid x."""
    x = np.asarray(x, float) / np.sum(x)
    Tm = S.p["Tm"]
    f = lambda T: _sle_residual(L, S, T, P, x)
    a, b = _bracket(f, guess if guess else Tm * 0.98, 1.0, 5.0, Tm * 1.5)
    return brentq(f, a, b, xtol=1e-10, rtol=1e-12)


def sle_composition(L, S, T, P):
    """Binary: liquid mole fraction of i in equilibrium with pure solid i.
    Returns None if solid i cannot exist at this T (liquid i is stable)."""
    i = S.index
    j = 1 - i

    def h(xi):
        x = np.zeros(2)
        x[i], x[j] = xi, 1 - xi
        return _sle_residual(L, S, T, P, x)
    if h(1.0) < 0:
        return None
    return brentq(h, 1e-300, 1.0, xtol=1e-16, rtol=1e-14)


@dataclass
class SLECurve:
    solid: str
    P: float
    x_solid: np.ndarray        # liquid mole fraction of the solid-forming component
    T: np.ndarray


def sle_curve(L, S, P, npts=41, x_min=0.02):
    xs = np.linspace(1.0, x_min, npts)
    out_x, out_T, guess = [], [], None
    i, j = S.index, 1 - S.index
    for xi in xs:
        x = np.zeros(2)
        x[i], x[j] = xi, 1 - xi
        try:
            T = sle_temperature(L, S, P, x, guess)
        except SolverError:
            continue
        guess = T
        out_x.append(xi)
        out_T.append(T)
    return SLECurve(S.names[i], P, np.array(out_x), np.array(out_T))


def eutectic(L, S1, S2, P):
    """Binary eutectic: liquid in equilibrium with both pure solids.
    Returns (x of S1's component, T) or None."""
    i = S1.index

    def F(xi):
        x = np.zeros(2)
        x[i], x[1 - i] = xi, 1 - xi
        return sle_temperature(L, S1, P, x) - sle_temperature(L, S2, P, x)
    try:
        a, b = 1e-6, 1 - 1e-6
        if F(a) * F(b) > 0:
            return None
        xi = brentq(F, a, b, xtol=1e-14)
    except SolverError:
        return None
    x = np.zeros(2)
    x[i], x[1 - i] = xi, 1 - xi
    return xi, sle_temperature(L, S1, P, x)


# ==========================================================================
# 6. External-field profiles
# ==========================================================================
@dataclass
class Profile:
    position: np.ndarray       # m
    P: np.ndarray              # Pa
    x: np.ndarray              # (npos, n) mole fractions of the local phase
    phase: list                # phase label at each position
    density: np.ndarray        # kg/m3 (nan if unknown)
    names: list
    kind: str                  # 'gravity' | 'centrifugal'
    notes: list = field(default_factory=list)
    extra: dict = field(default_factory=dict)


def _solve_at(phase, T, idx, c, n, pos, u0, M, fld):
    """At one position, solve  ln f_i(T, P, x) = c_i - M_i Phi(pos)/RT  for
    ln P and the free mole fractions of the components present (idx)."""
    target = c - M[idx] * fld.potential(pos) / (R * T)

    def full_x(free):
        x = np.zeros(n)
        if len(idx) > 1:
            x[idx[:-1]] = free
            x[idx[-1]] = 1 - np.sum(free)
        else:
            x[idx[0]] = 1.0
        return x

    def eqs(v):
        x = full_x(v[1:])
        if np.any(x[idx] <= 0):
            return np.full(len(idx), 1e3)
        try:
            return phase.ln_f(T, math.exp(v[0]), x)[idx] - target
        except (ValueError, OverflowError):
            return np.full(len(idx), 1e3)
    sol, info, ier, msg = fsolve(eqs, u0, full_output=True, xtol=1e-13)
    if ier != 1 or np.max(np.abs(info["fvec"])) > 1e-8:
        raise SolverError(f"no solution at {pos:.4g} m")
    return math.exp(sol[0]), full_x(sol[1:]), sol


def _step_to(phase, T, idx, c, n, prev, pos, u, M, fld, depth=0):
    """Solve at `pos` starting from the solution `u` at `prev`. The pressure
    guess is the hydrostatic estimate P - rho dPhi; if the solve still fails
    the step is halved (up to 2^12 sub-steps) before giving up, so a large
    step is never mistaken for the phase ceasing to exist."""
    u0 = np.array(u, float)
    P_prev = math.exp(u0[0])
    x_prev = np.zeros(n)
    if len(idx) > 1:
        x_prev[idx[:-1]] = u0[1:]
        x_prev[idx[-1]] = 1 - np.sum(u0[1:])
    else:
        x_prev[idx[0]] = 1.0
    rho = _density(phase, T, P_prev, x_prev, M)
    if np.isfinite(rho):
        P_guess = P_prev - rho * (fld.potential(pos) - fld.potential(prev))
        if P_guess > 0:
            u0[0] = math.log(P_guess)
    try:
        return _solve_at(phase, T, idx, c, n, pos, u0, M, fld)
    except SolverError:
        if depth >= 12:
            raise
        mid = 0.5 * (prev + pos)
        _, _, u_mid = _step_to(phase, T, idx, c, n, prev, mid, u, M, fld, depth + 1)
        return _step_to(phase, T, idx, c, n, mid, pos, u_mid, M, fld, depth + 1)


def _density(phase, T, P, x, M):
    try:
        v = phase.molar_volume(T, P, x)
    except (ValueError, SolverError):
        return np.nan
    return np.nan if v is None else float(np.dot(x, M) / v)


def _sat_pressure(L, V, T, x, kind):
    """Bubble P of a liquid (kind 'L') or dew P of a vapor (kind 'V'),
    used to test whether the local fluid is still a single phase."""
    try:
        return saturation(L, V, x, 0.0 if kind == "L" else 1.0, T=T).P
    except (SolverError, ValueError):
        return None


def _unstable(kind, P, Psat):
    """Liquid below its bubble P, or vapor above its dew P, splits."""
    if Psat is None:
        return False
    return P < Psat if kind == "L" else P > Psat


def _march_side(system, L, V, T, start, P_start, x_start, kind, positions, out, notes):
    """March one phase from `start` through `positions` (ordered away from
    start). If the phase reaches saturation, locate the contact by bisection
    and continue with the other phase from there."""
    fld, M, n = system.field, system.M, len(system.names)
    ph = L if kind == "L" else V
    x_start = np.asarray(x_start, float)
    idx = np.where(ph.present & (x_start > 0))[0]
    # the conserved quantity  ln f_i + M_i Phi / RT, evaluated at the start
    c = ph.ln_f(T, P_start, x_start)[idx] + M[idx] * fld.potential(start) / (R * T)
    u = np.concatenate([[math.log(P_start)], x_start[idx[:-1]]])
    prev = start
    for k, pos in enumerate(positions):
        crossed, P, x, u_new = False, None, None, None
        try:
            P, x, u_new = _step_to(ph, T, idx, c, n, prev, pos, u, M, fld)
            crossed = V is not None and _unstable(kind, P, _sat_pressure(L, V, T, x, kind))
        except SolverError:
            crossed = V is not None        # the phase cannot exist here: assume a contact
            if not crossed:
                raise
        if not crossed:
            out[pos] = (P, x, "liquid" if kind == "L" else "vapor",
                        _density(ph, T, P, x, M))
            prev, u = pos, u_new
            continue

        # ---- bisection for the contact between prev and pos ----
        a, b = prev, pos
        for _ in range(60):
            mid = 0.5 * (a + b)
            try:
                Pm, xm, _u = _step_to(ph, T, idx, c, n, prev, mid, u, M, fld)
                bad = _unstable(kind, Pm, _sat_pressure(L, V, T, xm, kind))
            except SolverError:
                bad = True
            if bad:
                b = mid
            else:
                a = mid
            if abs(b - a) < 1e-9 * max(1.0, abs(b)):
                break
        pc = a
        Pc, xc, _u = _step_to(ph, T, idx, c, n, prev, pc, u, M, fld)
        sat = saturation(L, V, xc, 0.0 if kind == "L" else 1.0, T=T)
        other_kind = "V" if kind == "L" else "L"
        x_other = sat.y if kind == "L" else sat.x
        notes.append(f"gas-liquid contact at {pc:.4g} m (P = {Pc / 1e5:.4g} bar): "
                     f"{'vapor' if other_kind == 'V' else 'liquid'} beyond this point")
        return _march_side(system, L, V, T, pc, sat.P, x_other, other_kind,
                           positions[k:], out, notes)
    return out


def column_profile(system, L, V, point, npts=21):
    """P, composition and phase along the field direction, from an
    equilibrium state given at the reference position.

    * Two phases at the reference: the interface is placed there; the denser
      liquid takes the side of lower potential Phi (below for gravity,
      outward for a centrifuge), the vapor the other side.
    * One phase at the reference: that phase is followed in both directions,
      and a second phase appears where the local fluid reaches saturation
      (e.g. a gas cap forming above a liquid as pressure falls)."""
    fld = system.field
    positions = np.linspace(min(fld.extent), max(fld.extent), npts)
    ref = fld.reference
    side_hi = [p for p in positions if p >= ref]                     # ascending
    side_lo = sorted((p for p in positions if p < ref), reverse=True)
    out, notes = {}, []
    T = point.T

    if point.phase == "two-phase":
        probe = ref + (max(fld.extent) - min(fld.extent)) * 1e-3
        hi_is_low_phi = fld.potential(probe) < 0
        sides = [(side_hi, "L" if hi_is_low_phi else "V"),
                 (side_lo, "V" if hi_is_low_phi else "L")]
        notes.append(f"interface at the reference position ({ref} m); "
                     f"liquid on the low-potential side")
        for side, kind in sides:
            if side:
                comp = point.x if kind == "L" else point.y
                _march_side(system, L, V, T, ref, point.P, comp, kind, side, out, notes)
    else:
        kind = "L" if point.phase == "liquid" else "V"
        comp = point.x if kind == "L" else point.y
        for side in (side_hi, side_lo):
            if side:
                _march_side(system, L, V, T, ref, point.P, comp, kind, side, out, notes)
        if not any("contact" in s for s in notes):
            notes.append("single phase throughout")

    pos_sorted = sorted(out)
    return Profile(np.array(pos_sorted),
                   np.array([out[p][0] for p in pos_sorted]),
                   np.array([out[p][1] for p in pos_sorted]),
                   [out[p][2] for p in pos_sorted],
                   np.array([out[p][3] for p in pos_sorted]),
                   system.names, fld.kind, notes, {"T": T})


def liquid_profile(system, L, T, P_ref, x_ref, npts=21):
    """Profile of a liquid with no vapor check (used for solid-liquid cases)."""
    fld = system.field
    positions = np.linspace(min(fld.extent), max(fld.extent), npts)
    ref = fld.reference
    out, notes = {}, []
    for side in ([p for p in positions if p >= ref],
                 sorted((p for p in positions if p < ref), reverse=True)):
        if side:
            _march_side(system, L, None, T, ref, P_ref, x_ref, "L", side, out, notes)
    pos_sorted = sorted(out)
    return Profile(np.array(pos_sorted),
                   np.array([out[p][0] for p in pos_sorted]),
                   np.array([out[p][1] for p in pos_sorted]),
                   [out[p][2] for p in pos_sorted],
                   np.array([out[p][3] for p in pos_sorted]),
                   system.names, fld.kind, notes, {"T": T})

# ==========================================================================
# 7. Dispatcher: solve whatever the checked input asks for
# ==========================================================================
@dataclass
class Result:
    equilibrium: str
    calculation: str
    data: object
    notes: list = field(default_factory=list)
    profiles: list = field(default_factory=list)


def solve(system):
    """Run the calculation implied by the input for every equilibrium that
    the checker found sufficient. Returns a list of Result."""
    st = system.state
    T, P, sweep = st["T"], st["P"], st["sweep"]
    xL = None if st["xL"] is None else np.array(st["xL"])
    xV = None if st["xV"] is None else np.array(st["xV"])
    z = None if st["z"] is None else np.array(st["z"])
    results = []
    solids = []

    for eq in system.equilibria:
        A, B = eq.phases
        res = Result(eq.label, eq.calculation or "", None)
        try:
            if eq.phase_names[1] == "V":
                _solve_vle(system, A, B, T, P, xL, xV, z, st["beta"], sweep, res)
            else:
                solids.append(B)
                _solve_sle(system, A, B, T, P, xL, sweep, res)
        except SolverError as e:
            res.notes.append(f"FAILED: {e}")
        results.append(res)

    if len(solids) == 2 and sweep and P is not None:
        L = system.equilibria[0].phases[0]
        eu = eutectic(L, solids[0], solids[1], P)
        if eu:
            xi, Te = eu
            nm = solids[0].names[solids[0].index]
            results.append(Result("Eutectic", "liquid + both pure solids",
                                  {"x": xi, "component": nm, "T": Te, "P": P}))
    return results


def _solve_vle(system, L, V, T, P, xL, xV, z, beta, sweep, res):
    fld = system.field
    point = None
    if z is not None:
        if T is not None and P is not None:
            point = flash_TP(L, V, T, P, z)
        else:
            point = saturation(L, V, z, beta, T=T, P=P)
        res.data = point
    elif sweep:
        if sweep == "x":
            res.data = vle_diagram(L, V, T=T, P=P)
        elif sweep == "T":
            res.data = vle_diagram(L, V, P=P)
        else:
            res.data = vle_diagram(L, V, T=T)
        res.notes += res.data.notes
        if fld:
            res.notes.append("the field does not change the local phase diagram "
                             "(Phi is the same for coexisting phases at one "
                             "position); profiles shown for three bubble points")
            for x1 in (0.25, 0.5, 0.75):
                x = np.array([x1, 1 - x1])
                pt = saturation(L, V, x, 0.0, T=T, P=P)
                pt.phase = "two-phase"     # incipient vapor at the bubble point
                prof = column_profile(system, L, V, pt)
                prof.extra.update(label=f"bubble point x1 = {x1}")
                res.profiles.append(prof)
        return
    elif xL is not None:
        point = bubble_P(L, V, T, xL) if T is not None else bubble_T(L, V, P, xL)
        res.data = point
    elif xV is not None:
        point = dew_P(L, V, T, xV) if T is not None else dew_T(L, V, P, xV)
        res.data = point

    if fld and point is not None:
        res.profiles.append(column_profile(system, L, V, point))


def _solve_sle(system, L, S, T, P, xL, sweep, res):
    if sweep:
        res.data = sle_curve(L, S, P)
        if isinstance(L, IdealSolutionLiquid):
            res.notes.append("ideal solution: real solutions (e.g. sugar-water) "
                             "can deviate strongly; an activity model is needed "
                             "for quantitative results")
        return
    if xL is not None:
        Tsle = sle_temperature(L, S, P, xL)
        res.data = {"T": Tsle, "P": P, "x": xL}
    else:
        xi = sle_composition(L, S, T, P)
        if xi is None:
            res.notes.append(f"T is above the melting point: no solid "
                             f"{S.names[S.index]} forms")
            res.data = None
            return
        x = np.zeros(2)
        x[S.index], x[1 - S.index] = xi, 1 - xi
        res.data = {"T": T, "P": P, "x": x}
    if system.field and res.data:
        d = res.data
        prof = liquid_profile(system, L, d["T"], d["P"], d["x"])
        i = S.index
        prof.extra["supersaturation"] = np.array(
            [L.ln_f(d["T"], Pp, xx)[i] - S.ln_f(d["T"], Pp)[i]
             for Pp, xx in zip(prof.P, prof.x)])
        prof.notes.append("supersaturation = ln(f_liquid / f_solid) of "
                          f"{S.names[i]}: > 0 means solid forms there")
        res.profiles.append(prof)


# ==========================================================================
# 8. Text output
# ==========================================================================
def _fmt_x(x, names):
    return ", ".join(f"{n} {v:.4f}" for n, v in zip(names, x) if np.isfinite(v))


def print_results(system, results):
    names = system.names
    for r in results:
        print(f"* {r.equilibrium}" + (f"  [{r.calculation}]" if r.calculation else ""))
        d = r.data
        if isinstance(d, VLEPoint):
            print(f"    T = {d.T:.3f} K, P = {d.P / 1e5:.5f} bar, state: {d.phase}"
                  + (f", vapor fraction {d.beta:.5f}" if d.phase == "two-phase" else ""))
            if np.all(np.isfinite(d.x)):
                print(f"    liquid: {_fmt_x(d.x, names)}")
            if np.all(np.isfinite(d.y)):
                print(f"    vapor : {_fmt_x(d.y, names)}")
        elif isinstance(d, Diagram):
            unit = "bar" if d.kind == "Pxy" else "K"
            fixed = (f"T = {d.fixed:.2f} K" if d.kind == "Pxy"
                     else f"P = {d.fixed / 1e5:.4f} bar")
            print(f"    {d.kind} diagram at {fixed}, {len(d.x1)} points")
            print(f"      x_{names[0]:<10s} y_{names[0]:<10s} {'P' if d.kind == 'Pxy' else 'T'} [{unit}]")
            step = max(1, len(d.x1) // 8)
            for k in list(range(0, len(d.x1), step)) + ([len(d.x1) - 1] if (len(d.x1) - 1) % step else []):
                v = d.var[k] / 1e5 if d.kind == "Pxy" else d.var[k]
                print(f"      {d.x1[k]:<12.4f} {d.y1[k]:<12.4f} {v:.4f}")
        elif isinstance(d, SLECurve):
            print(f"    solid {d.solid}: liquidus at P = {d.P / 1e5:.4f} bar, {len(d.T)} points")
            step = max(1, len(d.T) // 6)
            for k in range(0, len(d.T), step):
                print(f"      x_{d.solid} = {d.x_solid[k]:.3f}   T = {d.T[k]:.2f} K")
        elif isinstance(d, dict) and "component" in d:
            print(f"    x_{d['component']} = {d['x']:.5g}, T = {d['T']:.3f} K")
        elif isinstance(d, dict):
            print(f"    T = {d['T']:.3f} K, liquid: {_fmt_x(d['x'], names)}")
        for n in r.notes:
            print(f"    note: {n}")
        for prof in r.profiles:
            lab = prof.extra.get("label", "")
            print(f"    {prof.kind} profile {lab}:")
            print(f"      {'pos [m]':>9s} {'phase':>7s} {'P [bar]':>11s} "
                  f"{'x_' + names[0]:>12s} {'rho':>9s}"
                  + ("  supersat." if "supersaturation" in prof.extra else ""))
            step = max(1, len(prof.position) // 6)
            ks = list(range(0, len(prof.position), step))
            if ks[-1] != len(prof.position) - 1:
                ks.append(len(prof.position) - 1)
            for k in ks:
                line = (f"      {prof.position[k]:9.4f} {prof.phase[k]:>7s} "
                        f"{prof.P[k] / 1e5:11.5f} {prof.x[k][0]:12.6f} "
                        f"{prof.density[k]:9.2f}")
                if "supersaturation" in prof.extra:
                    line += f"  {prof.extra['supersaturation'][k]:+.2e}"
                print(line)
            for n in prof.notes:
                print(f"      note: {n}")


# ==========================================================================
# 9. Verification tests and demos
# ==========================================================================
def run_file(path):
    """Check an input file and, if it is sufficient, solve and print it."""
    from input_checker import check_input
    with open(path) as fh:
        text = fh.read()
    rep = check_input(text)
    print(rep)
    if rep.status != "SUFFICIENT":
        return
    print("\n" + "=" * 78 + "\nRESULTS\n" + "=" * 78)
    system = build_system(text)
    print_results(system, solve(system))


if __name__ == "__main__" and len(__import__("sys").argv) > 1:
    run_file(__import__("sys").argv[1])
elif __name__ == "__main__":
    from input_checker import DEMOS
    ok_all = True

    def report(name, ok, detail):
        global ok_all
        ok_all &= bool(ok)
        print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}")

    PR_TEXT = DEMOS["1. Propane / n-butane, PR, Pxy diagram at 300 K"]
    sysPR = build_system(PR_TEXT)
    L, V = sysPR.equilibria[0].phases
    T = 300.0

    print("1. Pure-component limits of the PR Pxy diagram (300 K)")
    d = vle_diagram(L, V, T=T)
    report("x_propane = 1 end", abs(d.var[-1] / 9.986e5 - 1) < 1e-3,
           f"{d.var[-1] / 1e5:.3f} bar (pure propane Psat from stage 2: 9.986 bar)")
    report("y = x at both ends", abs(d.y1[0] - d.x1[0]) < 1e-9 and abs(d.y1[-1] - d.x1[-1]) < 1e-9,
           "dew and bubble curves meet")
    report("bubble curve above dew curve", np.all(d.y1[1:-1] > d.x1[1:-1]),
           "propane (lighter) enriched in the vapor")

    print("\n2. Bubble P -> dew P round trip (x_propane = 0.4)")
    b = bubble_P(L, V, T, [0.4, 0.6])
    dw = dew_P(L, V, T, b.y, guess=b.P)
    report("same P and x back", abs(dw.P / b.P - 1) < 1e-8 and np.max(np.abs(dw.x - b.x)) < 1e-8,
           f"P = {b.P / 1e5:.5f} bar, y_propane = {b.y[0]:.5f}")
    bt = bubble_T(L, V, b.P, [0.4, 0.6])
    report("bubble T at that P returns 300 K", abs(bt.T - T) < 1e-6, f"{bt.T:.6f} K")

    print("\n3. Isothermal flash between the bubble and dew pressures")
    Pd = dew_P(L, V, T, [0.4, 0.6]).P
    Pm = 0.5 * (b.P + Pd)
    fl = flash_TP(L, V, T, Pm, [0.4, 0.6])
    mb = np.max(np.abs((1 - fl.beta) * fl.x + fl.beta * fl.y - np.array([0.4, 0.6])))
    lf = np.max(np.abs(L.ln_f(T, Pm, fl.x) - V.ln_f(T, Pm, fl.y)))
    report("material balance", mb < 1e-12, f"error {mb:.1e}, vapor fraction {fl.beta:.4f}")
    report("equal fugacities", lf < 1e-9, f"max |ln fL - ln fV| = {lf:.1e}")
    report("below dew P -> vapor", flash_TP(L, V, T, 0.9 * Pd, [0.4, 0.6]).phase == "vapor", "")
    report("above bubble P -> liquid", flash_TP(L, V, T, 1.1 * b.P, [0.4, 0.6]).phase == "liquid", "")

    print("\n4. Ideal solution: bubble P equals sum x_i Psat_i (ethanol-water, 350 K)")
    sysE = build_system(DEMOS["5. Ethanol / water, user Psat relations, centrifuge"])
    LE, VE = sysE.equilibria[0].phases
    bE = bubble_P(LE, VE, 350.0, [0.3, 0.7])
    ps = np.exp(LE.ln_f0(350.0, bE.P))
    # Poynting factor is included, so compare with it
    report("Raoult's law", abs(bE.P / (0.3 * ps[0] + 0.7 * ps[1]) - 1) < 1e-9,
           f"{bE.P / 1e5:.5f} bar")

    print("\n5. Sugar-water at 1 atm (hand values from stage 2)")
    sysS = build_system(DEMOS["3. Sugar / water at 1 atm (instructor's example)"])
    resS = solve(sysS)
    LS, VS = sysS.equilibria[0].phases
    bS = bubble_T(LS, VS, 101325.0, [0.95, 0.05])
    report("boiling point", abs(bS.T - 374.617) < 1e-3, f"{bS.T:.3f} K")
    ice = sysS.equilibria[1].phases[1]
    Tf = sle_temperature(LS, ice, 101325.0, [0.95, 0.05])
    report("freezing point", abs(Tf - 267.956) < 1e-3, f"{Tf:.3f} K")
    eu = [r for r in resS if r.equilibrium == "Eutectic"]
    report("eutectic found", len(eu) == 1, f"T = {eu[0].data['T']:.2f} K" if eu else "")

    print("\n6. Gravity: hydrostatic pressure follows from mu_i + M_i g z = const")

    def grav(height):
        return (PR_TEXT.replace("sweep = x", "P = 25 bar\nz = 0.4\nfield = gravity\n"
                                f"height = 0, {height} m\nreference = 0 m")
                .replace("T = 300 K", "T = 320 K")
                + "propane.M = 44.10 g/mol\nn-butane.M = 58.12 g/mol\n")
    sysG = build_system(grav(100))                # short column: liquid only
    LG, VG = sysG.equilibria[0].phases
    prof = solve(sysG)[0].profiles[0]
    dz = np.diff(prof.position)
    P_int = prof.P[0] - np.concatenate(
        [[0], np.cumsum(0.5 * (prof.density[1:] + prof.density[:-1]) * 9.80665 * dz)])
    err = np.max(np.abs(P_int - prof.P)) / 1e5
    report("100 m column stays liquid", set(prof.phase) == {"liquid"}, "")
    report("P(z) vs integral of rho g", err < 1e-3,
           f"max difference {err:.1e} bar over a {(prof.P[0] - prof.P[-1]) / 1e5:.3f} bar drop")
    inv = np.array([LG.ln_f(320.0, Pp, xx) + sysG.M * 9.80665 * zz / (R * 320.0)
                    for Pp, xx, zz in zip(prof.P, prof.x, prof.position)])
    report("ln f_i + M_i g z / RT constant", np.max(np.ptp(inv, axis=0)) < 1e-8,
           f"spread {np.max(np.ptp(inv, axis=0)):.1e}")
    report("heavier n-butane enriched at the bottom", prof.x[0][1] > prof.x[-1][1],
           f"x_n-butane {prof.x[0][1]:.6f} (0 m) -> {prof.x[-1][1]:.6f} (100 m)")

    grav_text = grav(500)                          # tall column: gas cap expected
    sysG5 = build_system(grav_text)
    LG5, VG5 = sysG5.equilibria[0].phases
    prof5 = solve(sysG5)[0].profiles[0]
    report("500 m column forms a gas cap", prof5.phase[0] == "liquid" and prof5.phase[-1] == "vapor",
           next((s for s in prof5.notes if "contact" in s), "no contact found"))
    inv5 = np.array([(LG5 if ph == "liquid" else VG5).ln_f(320.0, Pp, xx)
                     + sysG5.M * 9.80665 * zz / (R * 320.0)
                     for Pp, xx, zz, ph in zip(prof5.P, prof5.x, prof5.position, prof5.phase)])
    report("invariant holds across the contact", np.max(np.ptp(inv5, axis=0)) < 1e-7,
           f"spread {np.max(np.ptp(inv5, axis=0)):.1e}")

    print("\n7. Centrifuge profile does not depend on the number of grid points")
    LC, VC = sysE.equilibria[0].phases
    walls = []
    for npts in (3, 11, 51):
        pc = saturation(LC, VC, [0.5, 0.5], 0.0, T=350.0)
        pc.phase = "two-phase"
        pr = column_profile(sysE, LC, VC, pc, npts=npts)
        walls.append(pr.P[-1])
        spurious = [s for s in pr.notes if "contact" in s]
    report("same wall pressure for 3, 11, 51 points",
           np.ptp(walls) / walls[0] < 1e-9 and not spurious,
           f"{walls[0] / 1e5:.4f} bar at r = 0.10 m")

    print("\n" + ("ALL TESTS PASSED" if ok_all else "SOME TESTS FAILED"))

    print("\n" + "=" * 78 + "\nDEMO OUTPUT\n" + "=" * 78)
    for title, text in [("PR propane / n-butane", PR_TEXT),
                        ("Sugar-water", DEMOS["3. Sugar / water at 1 atm (instructor's example)"]),
                        ("Gravity column, propane / n-butane", grav_text),
                        ("Centrifuge, ethanol-water", DEMOS["5. Ethanol / water, user Psat relations, centrifuge"])]:
        print(f"\n### {title}")
        s = build_system(text)
        print_results(s, solve(s))
