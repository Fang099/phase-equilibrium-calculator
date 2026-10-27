"""
thermo_models.py  --  Stage 2 of the Phase Equilibrium Calculator

Takes the input that input_checker.py has accepted and builds, for every
phase, a model that returns the fugacity of each component, ln f_i(T, P, x).
Since mu_i = mu_i(ref) + RT ln f_i, equal ln f_i in all phases is the
equality of chemical potentials; the solver (stage 3) only uses this.

CUBIC EOS IN HELMHOLTZ FORM
All three cubics share  P = RT/(v - b) - a / ((v + d1 b)(v + d2 b)):
    vdW: d1 = d2 = 0      SRK: d1 = 1, d2 = 0      PR: d1 = 1 + sqrt2, d2 = 1 - sqrt2
Integrating (P - nRT/V) over volume gives the residual Helmholtz energy
(the fundamental relation A(T, V, n) minus the ideal-gas part):
    A_res / RT = -n ln(1 - B/V) - D / (RT B (d1 - d2)) ln[(V + d1 B)/(V + d2 B)]
    (vdW limit d1 = d2 = 0:  - D / (RT V))
with B = sum n_i b_i, D = n^2 a_mix = sum_i sum_j n_i n_j a_ij,
a_ij = sqrt(a_i a_j)(1 - k_ij)   (van der Waals one-fluid mixing rules).
Everything else is a derivative of this one function:
    P        = nRT/V - (dA_res/dV)_T,n
    ln phi_i = (d(A_res/RT)/dn_i)_T,V,nj - ln Z
    G_res/RT = A_res/(nRT) + Z - 1 - ln Z  = sum_i x_i ln phi_i

OTHER ROUTES (when an EOS is not used)
    Ideal gas               ln f_i = ln(y_i P)
    Liquid solution         ln f_i = ln x_i + ln gamma_i + ln f_i0,
                            f_i0 = Psat_i exp[V_L,i (P - Psat_i)/RT]   (Poynting)
    Activity coefficients come from an excess Gibbs energy G^E/RT
    (gamma_i = 1 for an ideal solution):
        Margules  G^E/RT = x1 x2 (A21 x1 + A12 x2)
        van Laar  G^E/RT = A12 A21 x1 x2 / (A12 x1 + A21 x2)
        NRTL      G^E/RT = sum_i x_i (sum_j tau_ji G_ji x_j) / (sum_k G_ki x_k),
                  tau_ij = dg_ij / RT,  G_ij = exp(-alpha_ij tau_ij)
    with A12 = ln gamma_1 at infinite dilution in 2 (and A21 the reverse);
    ln gamma_i = d(n G^E/RT)/dn_i, checked numerically in the tests.
    Pure solid i            ln f_i(S) = ln f_i0(L) + ln(f_S/f_L)    (Gibbs-Helmholtz):
        ln(f_S/f_L) = -(dHfus/R)(1/T - 1/Tm) + (dCp/R)(Tm/T - 1) - (dCp/R) ln(Tm/T)
                      + (V_S - V_L)(P - P0)/RT        (Tm, dHfus taken at P0 = 1 atm)
    A nonvolatile component has no vapor pressure; its pure-liquid reference
    fugacity is unknown, so it is set to 1 Pa at 1 atm. This is consistent
    because it is only ever compared with its own pure solid, where the
    reference cancels.

EXTERNAL FIELD
    Phi = g (z - z_ref)            gravity
    Phi = -omega^2 (r^2 - r_ref^2)/2   centrifugal
    Equilibrium: ln f_i + M_i Phi / RT is the same everywhere.

Run `python thermo_models.py` for the verification tests.
"""

import math
from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from input_checker import R, G_STANDARD, check_input, psat, psat_source

P0 = 101325.0   # Pa, pressure at which Tm, dHfus, Tb are taken


# ==========================================================================
# 1. Cubic equations of state, written as a Helmholtz energy
# ==========================================================================
def _m_pr(w):
    # Peng-Robinson (1976); 1978 form for heavy components (omega > 0.49)
    if w <= 0.49:
        return 0.37464 + 1.54226 * w - 0.26992 * w ** 2
    return 0.379642 + 1.48503 * w - 0.164423 * w ** 2 + 0.016666 * w ** 3


EOS_CONSTANTS = {
    "vdW": dict(d1=0.0, d2=0.0, Oa=27 / 64, Ob=1 / 8, m=None),
    "SRK": dict(d1=1.0, d2=0.0, Oa=0.42748023, Ob=0.08664035,
                m=lambda w: 0.480 + 1.574 * w - 0.176 * w ** 2),
    "PR": dict(d1=1 + math.sqrt(2), d2=1 - math.sqrt(2), Oa=0.45723553,
               Ob=0.07779607, m=_m_pr),
}


class CubicEOS:
    """vdW, SRK or PR for a mixture, built on A_res(T, V, n)."""

    def __init__(self, model, Tc=None, Pc=None, omega=None, kij=None,
                 a_vdw=None, b_vdw=None):
        if model not in EOS_CONSTANTS:
            raise ValueError(f"unknown EOS '{model}'")
        self.model = model
        c = EOS_CONSTANTS[model]
        self.d1, self.d2 = c["d1"], c["d2"]
        n = len(Tc) if Tc is not None else len(a_vdw)
        self.n = n
        self.Tc = np.array(Tc, float) if Tc is not None else None
        self.Pc = np.array(Pc, float) if Pc is not None else None
        self.kij = np.zeros((n, n)) if kij is None else np.array(kij, float)

        # b_i, and a_i at the critical point (temperature dependence via alpha)
        self.ac = np.zeros(n)
        self.b = np.zeros(n)
        self.mk = np.zeros(n)
        for i in range(n):
            if a_vdw is not None and a_vdw[i] is not None:      # vdW a, b given
                self.ac[i], self.b[i] = a_vdw[i], b_vdw[i]
            else:
                self.ac[i] = c["Oa"] * (R * self.Tc[i]) ** 2 / self.Pc[i]
                self.b[i] = c["Ob"] * R * self.Tc[i] / self.Pc[i]
            if c["m"] is not None:
                self.mk[i] = c["m"](omega[i])
        # critical constants and omega are kept for initial estimates (Wilson K)
        self.omega = np.array(omega if omega is not None else [0.0] * n, float)
        if self.Tc is None:                      # from vdW a, b
            self.Tc = 8 * self.ac / (27 * R * self.b)
            self.Pc = self.ac / (27 * self.b ** 2)

    # ---- pure and mixture parameters ----
    def a_pure(self, T):
        if self.model == "vdW":
            return self.ac.copy()
        alpha = (1 + self.mk * (1 - np.sqrt(T / self.Tc))) ** 2
        return self.ac * alpha

    def a_matrix(self, T):
        a = self.a_pure(T)
        return np.sqrt(np.outer(a, a)) * (1 - self.kij)

    def mix(self, T, n):
        """Return D = n^2 a_mix, its gradient dD/dn_i, and B = n b_mix."""
        aij = self.a_matrix(T)
        Di = 2 * aij @ n
        D = 0.5 * n @ Di
        B = self.b @ n
        return D, Di, B

    # ---- the fundamental relation and its derivatives (T, V, n) ----
    def _log_term(self, V, B):
        """g = ln[(V + d1 B)/(V + d2 B)] / (d1 - d2); vdW limit: B / V."""
        if self.model == "vdW":
            return B / V
        return math.log((V + self.d1 * B) / (V + self.d2 * B)) / (self.d1 - self.d2)

    def helmholtz_res(self, T, V, n):
        """Residual Helmholtz energy A_res(T, V, n) in J (V in m3, n in mol)."""
        n = np.asarray(n, float)
        D, _, B = self.mix(T, n)
        ntot = n.sum()
        F = -ntot * math.log(1 - B / V) - D / (R * T * B) * self._log_term(V, B)
        return R * T * F

    def pressure(self, T, V, n):
        """P = nRT/V - dA_res/dV, written out in closed form."""
        n = np.asarray(n, float)
        D, _, B = self.mix(T, n)
        return (n.sum() * R * T / (V - B)
                - D / ((V + self.d1 * B) * (V + self.d2 * B)))

    def ln_phi_TV(self, T, V, n):
        """ln phi_i = d(A_res/RT)/dn_i - ln Z, analytic derivatives."""
        n = np.asarray(n, float)
        D, Di, B = self.mix(T, n)
        ntot = n.sum()
        b = self.b
        P = self.pressure(T, V, n)
        Z = P * V / (ntot * R * T)
        dF1 = -math.log(1 - B / V) + ntot * b / (V - B)
        if self.model == "vdW":
            dF2 = -Di / (R * T * V)
        else:
            d1, d2 = self.d1, self.d2
            g = self._log_term(V, B)
            # d/dn_i of g = b_i (d1/(V+d1B) - d2/(V+d2B)) / (d1 - d2)
            dln = b * (d1 / (V + d1 * B) - d2 / (V + d2 * B)) / (d1 - d2)
            dF2 = -(Di / (R * T * B) - D * b / (R * T * B ** 2)) * g \
                  - D / (R * T * B) * dln
        return dF1 + dF2 - math.log(Z)

    # ---- (T, P, x) interface: solve the cubic for Z ----
    def z_roots(self, T, P, x):
        """Real compressibility-factor roots with v > b, ascending."""
        x = np.asarray(x, float)
        D, _, B = self.mix(T, x)            # per mole of mixture
        A_ = D * P / (R * T) ** 2
        B_ = B * P / (R * T)
        d1, d2 = self.d1, self.d2
        s, p = d1 + d2, d1 * d2
        coeffs = [1.0,
                  (s - 1) * B_ - 1,
                  A_ + p * B_ ** 2 - s * B_ * (B_ + 1),
                  -(A_ * B_ + p * B_ ** 2 * (B_ + 1))]
        roots = np.roots(coeffs)
        real = sorted(r.real for r in roots
                      if abs(r.imag) < 1e-9 * max(1, abs(r.real)) and r.real > B_)
        return real

    def molar_volume(self, T, P, x, phase="stable"):
        """Molar volume (m3/mol). phase: 'liquid' (smallest root),
        'vapor' (largest root) or 'stable' (root with the lowest G_res)."""
        Z = self.z_roots(T, P, x)
        if not Z:
            raise ValueError("no physical root of the cubic")
        if phase == "liquid":
            z = Z[0]
        elif phase == "vapor":
            z = Z[-1]
        else:
            x = np.asarray(x, float)
            z = min(Z, key=lambda zz: x @ self.ln_phi_TV(T, zz * R * T / P, x))
        return z * R * T / P

    def ln_phi(self, T, P, x, phase="stable"):
        x = np.asarray(x, float) / np.sum(x)
        v = self.molar_volume(T, P, x, phase)
        return self.ln_phi_TV(T, v, x)


# ==========================================================================
# 2. Phase models: common interface  ln_f(T, P, x) -> array, molar_volume
# ==========================================================================
def ln_psat(p, T):
    """ln(Psat / Pa), computed in log form so it stays finite at low T."""
    from input_checker import from_si, to_si
    if "antoine" in p:
        a = p["antoine"]
        t = from_si(T, a["T_unit"], "temperature")
        arg = a["A"] - a["B"] / (t + a["C"])
        ln_val = arg * math.log(10) if a["base"] == "log10" else arg
        return ln_val + math.log(to_si(1.0, a["P_unit"], "pressure"))
    if "Psat" in p:
        return math.log(psat(p, T))
    return math.log(P0) - p["dHvap"] / R * (1 / T - 1 / p["Tb"])


def _safe_log(x):
    x = np.asarray(x, float)
    with np.errstate(divide="ignore"):
        return np.where(x > 0, np.log(np.where(x > 0, x, 1.0)), -np.inf)


class EOSPhase:
    """Liquid or vapor described by a cubic EOS."""

    def __init__(self, eos, kind, names):
        self.eos, self.kind, self.names = eos, kind, names
        self.present = np.ones(len(names), bool)

    def ln_phi(self, T, P, x):
        return self.eos.ln_phi(T, P, x, "liquid" if self.kind == "L" else "vapor")

    def ln_f(self, T, P, x):
        return _safe_log(x) + self.ln_phi(T, P, x) + math.log(P)

    def molar_volume(self, T, P, x):
        return self.eos.molar_volume(T, P, x, "liquid" if self.kind == "L" else "vapor")


class IdealGasPhase:
    """Ideal-gas vapor; nonvolatile components are absent."""

    def __init__(self, names, present):
        self.kind, self.names = "V", names
        self.present = np.asarray(present, bool)

    def ln_phi(self, T, P, x):
        return np.zeros(len(self.names))

    def ln_f(self, T, P, x):
        lf = _safe_log(x) + math.log(P)
        return np.where(self.present, lf, -np.inf)

    def molar_volume(self, T, P, x):
        return R * T / P


# --------------------------------------------------------------------------
# Activity-coefficient models (excess Gibbs energy)
# --------------------------------------------------------------------------
class Margules:
    """Two-parameter Margules model for a binary; A12, A21 dimensionless."""
    name = "Margules"

    def __init__(self, A12, A21):
        self.A12, self.A21 = A12, A21

    def ge_rt(self, T, x):
        x1, x2 = x
        return x1 * x2 * (self.A21 * x1 + self.A12 * x2)

    def ln_gamma(self, T, x):
        x1, x2 = x
        A12, A21 = self.A12, self.A21
        return np.array([x2 ** 2 * (A12 + 2 * (A21 - A12) * x1),
                         x1 ** 2 * (A21 + 2 * (A12 - A21) * x2)])


class VanLaar:
    """van Laar model for a binary; A12, A21 dimensionless and positive
    (or both negative)."""
    name = "van Laar"

    def __init__(self, A12, A21):
        if A12 * A21 <= 0:
            raise ValueError("van Laar needs A12 and A21 of the same sign")
        self.A12, self.A21 = A12, A21

    def ge_rt(self, T, x):
        x1, x2 = x
        d = self.A12 * x1 + self.A21 * x2
        return self.A12 * self.A21 * x1 * x2 / d

    def ln_gamma(self, T, x):
        x1, x2 = x
        A12, A21 = self.A12, self.A21
        d = A12 * x1 + A21 * x2
        return np.array([A12 * (A21 * x2 / d) ** 2, A21 * (A12 * x1 / d) ** 2])


class NRTL:
    """NRTL for any number of components. dg[i, j] = g_ij - g_jj (J/mol),
    alpha symmetric; tau_ij = dg_ij / RT."""
    name = "NRTL"

    def __init__(self, dg, alpha):
        self.dg = np.asarray(dg, float)
        self.alpha = np.asarray(alpha, float)

    def _tg(self, T):
        tau = self.dg / (R * T)
        return tau, np.exp(-self.alpha * tau)

    def ge_rt(self, T, x):
        x = np.asarray(x, float)
        tau, G = self._tg(T)
        num = (tau * G).T @ x             # sum_j tau_ji G_ji x_j
        den = G.T @ x                     # sum_k G_ki x_k
        return float(x @ (num / den))

    def ln_gamma(self, T, x):
        x = np.asarray(x, float)
        tau, G = self._tg(T)
        den = G.T @ x                     # den_j = sum_k x_k G_kj
        num = (tau * G).T @ x             # num_j = sum_m x_m tau_mj G_mj
        term1 = num / den
        term2 = (G * (tau - term1[None, :])) @ (x / den)
        return term1 + term2


class SolutionLiquid:
    """Liquid solution: pure-liquid reference fugacities f_i0 from vapor
    pressures, corrected by activity coefficients (ideal if activity=None)."""

    def __init__(self, names, props, activity=None):
        self.kind, self.names, self.props = "L", names, props
        self.activity = activity
        self.present = np.ones(len(names), bool)
        self.absolute = [psat_source(p) is not None for p in props]

    def ln_gamma(self, T, x):
        if self.activity is None:
            return np.zeros(len(self.names))
        x = np.clip(np.asarray(x, float), 0.0, None)
        return self.activity.ln_gamma(T, x / x.sum())

    def ge_rt(self, T, x):
        return 0.0 if self.activity is None else self.activity.ge_rt(T, x)

    def ln_f0(self, T, P):
        """ln of the pure-liquid fugacity of each component at T, P."""
        out = np.zeros(len(self.names))
        for i, p in enumerate(self.props):
            VL = p.get("V_L", 0.0)
            if self.absolute[i]:
                lps = ln_psat(p, T)          # log form: no underflow at low T
                out[i] = lps + VL * (P - math.exp(lps)) / (R * T)    # Poynting
            else:
                out[i] = VL * (P - P0) / (R * T)    # arbitrary reference, 1 Pa at P0
        return out

    def ln_f(self, T, P, x):
        return _safe_log(x) + self.ln_gamma(T, x) + self.ln_f0(T, P)

    def ln_phi(self, T, P, x):
        """Effective fugacity coefficient, f_i / (x_i P) = gamma_i f_i0 / P."""
        return self.ln_gamma(T, x) + self.ln_f0(T, P) - math.log(P)

    def molar_volume(self, T, P, x):
        VL = [p.get("V_L") for p in self.props]
        if any(v is None for v in VL):
            return None
        return float(np.dot(x, VL))       # ideal mixing of volumes


IdealSolutionLiquid = SolutionLiquid      # name used by earlier code


class PureSolidPhase:
    """Pure solid of component i, referenced to the same pure liquid that
    the liquid-phase model uses, so the reference cancels at equilibrium."""

    def __init__(self, index, names, props, liquid_ref):
        self.kind, self.names, self.index = f"S({names[index]})", names, index
        self.p, self.liquid_ref = props[index], liquid_ref
        self.present = np.zeros(len(names), bool)
        self.present[index] = True

    def ln_ratio(self, T, P):
        """ln(f_solid / f_pure_liquid) at T, P."""
        p = self.p
        Tm, dH, dCp = p["Tm"], p["dHfus"], p.get("dCp_fus", 0.0)
        val = (-dH / R * (1 / T - 1 / Tm)
               + dCp / R * (Tm / T - 1) - dCp / R * math.log(Tm / T))
        if "V_S" in p and "V_L" in p:
            val += (p["V_S"] - p["V_L"]) * (P - P0) / (R * T)
        return val

    def ln_f(self, T, P, x=None):
        out = np.full(len(self.names), -np.inf)
        out[self.index] = self.liquid_ref(T, P)[self.index] + self.ln_ratio(T, P)
        return out

    def molar_volume(self, T, P, x=None):
        return self.p.get("V_S")


# ==========================================================================
# 3. External field
# ==========================================================================
@dataclass
class Field:
    kind: str                         # 'gravity' | 'centrifugal'
    strength: float                   # g (m/s2) or omega (rad/s)
    extent: tuple                     # (start, end) in m
    reference: float                  # position where the given state applies

    def potential(self, pos):
        """Phi(pos) - Phi(reference), J/kg."""
        if self.kind == "gravity":
            return self.strength * (pos - self.reference)
        return -0.5 * self.strength ** 2 * (pos ** 2 - self.reference ** 2)

    def field_term(self, M, T, pos):
        """M_i Phi / RT for each component (dimensionless)."""
        return np.asarray(M) * self.potential(pos) / (R * T)


# ==========================================================================
# 4. Build everything from the checker's report
# ==========================================================================
@dataclass
class Equilibrium:
    label: str
    phase_names: tuple                 # e.g. ('L', 'V') or ('L', 'S(water)')
    phases: tuple                      # the two phase-model objects
    route: str
    calculation: Optional[str]


@dataclass
class ThermoSystem:
    names: list
    M: Optional[np.ndarray]
    equilibria: list
    field: Optional[Field]
    state: dict                        # T, P, x(L), x(V), z, beta, sweep (SI)
    eos: Optional[CubicEOS] = None
    notes: list = field(default_factory=list)

    def describe(self):
        out = [f"Components: {', '.join(self.names)}"]
        if self.eos:
            out.append(f"EOS: {self.eos.model}")
        if self.field:
            f = self.field
            unit = "m/s2" if f.kind == "gravity" else "rad/s"
            out.append(f"Field: {f.kind}, strength {f.strength:.4g} {unit}, "
                       f"extent {f.extent} m, reference {f.reference} m")
        for e in self.equilibria:
            out.append(f"- {e.label}: phases {e.phase_names}, route: {e.route}")
            if e.calculation:
                out.append(f"    -> {e.calculation}")
        given = {k: v for k, v in self.state.items() if v is not None}
        out.append(f"State: {given}")
        return "\n".join(out)


def build_activity(inp):
    """Activity model chosen in the input, or None for an ideal solution."""
    names, n = inp.components, len(inp.components)
    if inp.activity == "ideal":
        return None
    get = lambda i, j, q: inp.opairs[(names[i], names[j])][q]
    if inp.activity in ("margules", "vanlaar"):
        cls = Margules if inp.activity == "margules" else VanLaar
        return cls(get(0, 1, "A"), get(1, 0, "A"))
    dg = np.zeros((n, n))
    alpha = np.full((n, n), 0.3)
    for i in range(n):
        for j in range(n):
            if i != j:
                dg[i, j] = get(i, j, "dg")
                a = inp.pairs.get(frozenset((names[i], names[j])), {}).get("alpha")
                if a is not None:
                    alpha[i, j] = a
    return NRTL(dg, alpha)


def build_system(text):
    """Check the input; if sufficient, build the phase models.
    Raises ValueError with the checker's report if not."""
    rep = check_input(text)
    if rep.status != "SUFFICIENT":
        raise ValueError("input is not sufficient:\n" + str(rep))
    inp = rep.inputs
    names = inp.components
    props = [inp.props[c] for c in names]
    n = len(names)

    M = None
    if all("M" in p for p in props):
        M = np.array([p["M"] for p in props])

    # EOS (built once if any equilibrium uses it)
    eos = None
    if any(e.route_kind == "eos" for e in rep.equilibria if e.status == "SUFFICIENT"):
        kij = np.zeros((n, n))
        for i in range(n):
            for j in range(i + 1, n):
                k = inp.pairs.get(frozenset((names[i], names[j])), {}).get("kij", 0.0)
                kij[i, j] = kij[j, i] = k
        use_ab = [("a" in p and "b" in p and not ("Tc" in p and "Pc" in p))
                  for p in props]
        if inp.model == "vdW" and any(use_ab):
            # mix of (a, b) and (Tc, Pc) inputs: convert Tc, Pc to a, b
            a = [p["a"] if u else 27 / 64 * (R * p["Tc"]) ** 2 / p["Pc"]
                 for p, u in zip(props, use_ab)]
            b = [p["b"] if u else R * p["Tc"] / (8 * p["Pc"])
                 for p, u in zip(props, use_ab)]
            eos = CubicEOS("vdW", a_vdw=a, b_vdw=b, kij=kij)
        else:
            eos = CubicEOS(inp.model, Tc=[p["Tc"] for p in props],
                           Pc=[p["Pc"] for p in props],
                           omega=[p.get("omega", 0.0) for p in props], kij=kij)

    liquid_ideal = SolutionLiquid(names, props, build_activity(inp))
    equilibria = []
    for e in rep.equilibria:
        if e.status != "SUFFICIENT":
            continue
        if e.route_kind == "eos":
            phases = (EOSPhase(eos, "L", names), EOSPhase(eos, "V", names))
        elif e.route_kind == "vp":
            present = [not p.get("nonvolatile") for p in props]
            phases = (liquid_ideal, IdealGasPhase(names, present))
        else:                                    # pure solid
            idx = names.index(e.phases[1][2:-1])
            phases = (liquid_ideal,
                      PureSolidPhase(idx, names, props, liquid_ideal.ln_f0))
        equilibria.append(Equilibrium(e.label, e.phases, phases, e.route, e.calc))

    fld = None
    if inp.field_kind:
        if inp.field_kind == "gravity":
            strength = inp.g if inp.g is not None else G_STANDARD
            ext = tuple(inp.height)
        else:
            strength, ext = inp.omega, tuple(inp.radius)
        ref = inp.reference if inp.reference is not None else ext[0]
        fld = Field(inp.field_kind, strength, ext, ref)

    state = dict(T=inp.T, P=inp.P, xL=inp.comp.get("L"), xV=inp.comp.get("V"),
                 z=inp.comp.get("z"), beta=inp.beta, sweep=inp.sweep)
    return ThermoSystem(names, M, equilibria, fld, state, eos,
                        notes=[w for w in rep.warnings])


# ==========================================================================
# 5. Verification tests
# ==========================================================================
def _num_checks(eos, T, V, n, h=1e-6):
    """Compare analytic P and ln phi with numerical derivatives of A_res."""
    n = np.asarray(n, float)
    A = lambda VV, nn: eos.helmholtz_res(T, VV, nn)
    dAdV = (A(V * (1 + h), n) - A(V * (1 - h), n)) / (2 * V * h)
    P_num = n.sum() * R * T / V - dAdV
    P_an = eos.pressure(T, V, n)
    Z = P_an * V / (n.sum() * R * T)
    lnphi_num = np.zeros(len(n))
    for i in range(len(n)):
        dn = np.zeros(len(n))
        dn[i] = h * n.sum()
        lnphi_num[i] = (A(V, n + dn) - A(V, n - dn)) / (2 * dn[i] * R * T) - math.log(Z)
    lnphi_an = eos.ln_phi_TV(T, V, n)
    return abs(P_num / P_an - 1), np.max(np.abs(lnphi_num - lnphi_an))


def _bisect(f, a, b, tol=1e-10):
    fa = f(a)
    for _ in range(200):
        m = 0.5 * (a + b)
        fm = f(m)
        if fa * fm <= 0:
            b = m
        else:
            a, fa = m, fm
        if b - a < tol:
            break
    return 0.5 * (a + b)


if __name__ == "__main__":
    from input_checker import DEMOS
    ok_all = True

    def report(name, ok, detail):
        global ok_all
        ok_all &= ok
        print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}")

    # propane / n-butane constants (typical literature values)
    Tc, Pc, w = [369.8, 425.1], [42.48e5, 37.96e5], [0.152, 0.200]
    n_mix = np.array([0.4, 0.6])

    print("1. Helmholtz-energy derivatives: analytic vs numerical")
    for model in ("vdW", "SRK", "PR"):
        eos = CubicEOS(model, Tc, Pc, w, kij=[[0, 0.003], [0.003, 0]])
        # compressed liquid at 20 bar, vapor at 3 bar (both single-root states)
        for label, P_, ph in (("liquid", 20e5, "liquid"), ("vapor", 3e5, "vapor")):
            V = eos.molar_volume(300.0, P_, n_mix, ph)
            eP, ephi = _num_checks(eos, 300.0, V, n_mix)
            report(f"{model:3s} {label:11s}", eP < 1e-7 and ephi < 1e-7,
                   f"rel. error in P {eP:.1e}, abs. error in ln phi {ephi:.1e}")

    print("\n2. Identity  sum x_i ln phi_i = A_res/nRT + Z - 1 - ln Z")
    for model in ("vdW", "SRK", "PR"):
        eos = CubicEOS(model, Tc, Pc, w)
        T, P = 300.0, 5e5
        v = eos.molar_volume(T, P, n_mix, "vapor")
        Z = P * v / (R * T)
        lhs = n_mix @ eos.ln_phi_TV(T, v, n_mix)
        rhs = eos.helmholtz_res(T, v, n_mix) / (R * T) + Z - 1 - math.log(Z)
        report(f"{model:3s}", abs(lhs - rhs) < 1e-10, f"difference {abs(lhs - rhs):.1e}")

    print("\n3. Pure propane at 300 K: saturation pressure from equal fugacities")
    for model in ("SRK", "PR"):
        eos = CubicEOS(model, [369.8], [42.48e5], [0.152])
        P = 9e5
        for _ in range(100):                     # P_new = P * phi_L / phi_V
            lL = eos.ln_phi(300.0, P, [1.0], "liquid")[0]
            lV = eos.ln_phi(300.0, P, [1.0], "vapor")[0]
            P *= math.exp(lL - lV)
            if abs(lL - lV) < 1e-12:
                break
        report(f"{model:3s}", abs(P / 9.98e5 - 1) < 0.02,
               f"Psat = {P / 1e5:.3f} bar (experimental about 9.98 bar)")

    print("\n4. Latent-heat route: pure-solid fugacity (benzene, x = 0.95, 1 atm)")
    BN = DEMOS["3. Benzene / naphthalene at 1 atm from latent heats only"]
    sysw = build_system(BN)
    x = np.array([0.95, 0.05])
    P = 101325.0
    sle = next(e for e in sysw.equilibria if e.phase_names[1] == "S(benzene)")
    L, S = sle.phases
    Tf = _bisect(lambda T: L.ln_f(T, P, x)[0] - S.ln_f(T, P)[0], 240, 278.68)
    Tf_hand = 1 / (1 / 278.68 - R * math.log(0.95) / 9870)
    report("freezing point", abs(Tf - Tf_hand) < 1e-6,
           f"{Tf:.3f} K (hand calculation {Tf_hand:.3f} K)")
    sysn = build_system(BN + "naphthalene.nonvolatile = yes\n")
    L, V = next(e for e in sysn.equilibria if e.phase_names[1] == "V").phases
    y = np.array([1.0, 0.0])
    Tb = _bisect(lambda T: L.ln_f(T, P, x)[0] - V.ln_f(T, P, y)[0], 340, 380)
    Tb_hand = 1 / (1 / 353.24 + R * math.log(0.95) / 30720)
    report("boiling point with a nonvolatile solute", abs(Tb - Tb_hand) < 1e-6,
           f"{Tb:.3f} K (hand calculation {Tb_hand:.3f} K)")
    report("nonvolatile component absent from vapor", V.ln_f(Tb, P, y)[1] == -np.inf,
           "ln f = -inf")

    print("\n4b. Activity models: ln gamma_i = d(n G^E/RT)/dn_i, Gibbs-Duhem, limits")
    dgm = np.array([[0, 2000, -500], [800, 0, 1500], [300, 1200, 0]], float)
    alm = np.array([[0, .3, .4], [.3, 0, .2], [.4, .2, 0]])
    cases = [(Margules(1.2, 0.7), [0.3, 0.7]), (VanLaar(1.68, 0.92), [0.3, 0.7]),
             (NRTL(dgm, alm), [0.2, 0.5, 0.3])]
    for m, xm in cases:
        xm = np.array(xm)
        num = np.zeros(len(xm))
        for i in range(len(xm)):
            dn = np.zeros(len(xm))
            dn[i] = 1e-6
            F = lambda nn: nn.sum() * m.ge_rt(340.0, nn / nn.sum())
            num[i] = (F(xm + dn) - F(xm - dn)) / 2e-6
        err = np.max(np.abs(num - m.ln_gamma(340.0, xm)))
        # Gibbs-Duhem: sum_i x_i d(ln gamma_i) = 0 along any composition change
        dx = np.zeros(len(xm))
        dx[0], dx[1] = 1e-6, -1e-6
        gd = xm @ (m.ln_gamma(340.0, xm + dx) - m.ln_gamma(340.0, xm - dx)) / 2e-6
        report(f"{m.name:9s}", err < 1e-8 and abs(gd) < 1e-7,
               f"derivative error {err:.1e}, Gibbs-Duhem residual {abs(gd):.1e}")
    vl = VanLaar(1.68, 0.92)
    report("infinite dilution: ln gamma_1(x1 -> 0) = A12",
           abs(vl.ln_gamma(340.0, [0.0, 1.0])[0] - 1.68) < 1e-12
           and abs(Margules(1.2, 0.7).ln_gamma(340.0, [0.0, 1.0])[0] - 1.2) < 1e-12, "")

    print("\n5. Field potential (centrifuge, 20000 rpm, r from 0.05 to 0.10 m)")
    sysc = build_system(DEMOS["5. Ethanol / water, user Psat relations, centrifuge"])
    f = sysc.field
    phi = f.potential(0.10)
    hand = -0.5 * (20000 * 2 * math.pi / 60) ** 2 * (0.10 ** 2 - 0.05 ** 2)
    report("Phi(0.10 m)", abs(phi / hand - 1) < 1e-12, f"{phi:.4g} J/kg")
    print("   M_i Phi / RT at the wall, 350 K:",
          np.round(f.field_term(sysc.M, 350.0, 0.10), 4), "for", sysc.names)

    print("\n6. build_system on the PR demo")
    print(build_system(DEMOS["1. Propane / n-butane, PR, Pxy diagram at 300 K"]).describe())

    print("\nALL TESTS PASSED" if ok_all else "\nSOME TESTS FAILED")
