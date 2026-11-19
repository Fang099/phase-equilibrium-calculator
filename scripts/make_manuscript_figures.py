"""Figures and table values for the manuscript, computed with the project's code.

    python make_figures.py /path/to/phase-equilibrium-calculator

writes figures/*.pdf and numbers.json (every number quoted in the text).
"""
import json
import math
import os
import sys
import warnings

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CODE = sys.argv[1] if len(sys.argv) > 1 else "../Codes"
sys.path.insert(0, os.path.abspath(CODE))
warnings.filterwarnings("ignore")
np.seterr(all="ignore")

import compound_db                                    # noqa: E402
import equilibrium_solver as es                       # noqa: E402
import pure_fluid as pf                               # noqa: E402
import web_api                                        # noqa: E402
from input_checker import check_input                 # noqa: E402
from thermo_models import build_system                # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
os.makedirs(FIG, exist_ok=True)
NUM = {}

plt.rcParams.update({"font.size": 8, "axes.labelsize": 8, "legend.fontsize": 6.5, "xtick.labelsize": 7,
                     "ytick.labelsize": 7, "axes.linewidth": 0.6, "lines.linewidth": 1.2,
                     "font.family": "DejaVu Sans", "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
                     "axes.spines.top": False, "axes.spines.right": False})
C_L, C_V, C_3, C_K, C_G = "#2a6fdb", "#e8590c", "#1baf7a", "#222831", "#8a94a0"
ONE = 3.33     # single column width, inches
TWO = 7.0


VAL = os.path.join(CODE, "validation")


def read_csv(name):
    import csv as _csv
    with open(os.path.join(VAL, name)) as fh:
        return list(_csv.DictReader([r for r in fh if not r.startswith("#")]))


REF = read_csv("reference_saturation.csv")


def ref_of(cid):
    rows = [r for r in REF if r["compound"] == cid]
    return (np.array([float(r["T_K"]) for r in rows]), np.array([float(r["Psat_Pa"]) for r in rows]),
            np.array([float(r["rhoL_kg_m3"]) for r in rows]), np.array([float(r["rhoV_kg_m3"]) for r in rows]))


def save(fig, name):
    fig.savefig(os.path.join(FIG, name + ".pdf"))
    fig.savefig(os.path.join(FIG, name + ".png"), dpi=200)
    plt.close(fig)


def label(ax, s, x=-0.02, y=1.02):
    ax.text(x, y, s, transform=ax.transAxes, va="bottom", ha="right", fontweight="bold", fontsize=8)


# ---------------------------------------------------------------------------
# 1. pure fluids: vapor pressure against Antoine, isotherms, densities
# ---------------------------------------------------------------------------
fig, axs = plt.subplots(1, 3, figsize=(TWO, 2.2))
ax = axs[0]
for cid, col in (("water", C_L), ("propane", C_V), ("benzene", C_3)):
    rec = compound_db.find(cid)
    for m, ls in (("PR", "-"), ("vdW", ":")):
        e, M, _ = pf.make_eos(cid, m)
        cur = pf.saturation_curve(e, npts=50)
        ax.plot(1000 / cur.T[:-1], cur.P[:-1] / 1e5, ls, color=col, lw=1.1 if m == "PR" else 0.9)
    Tr_, Pr_, _, _ = ref_of(cid)
    ax.plot(1000 / Tr_, Pr_ / 1e5, "o", ms=2.6, mfc="none", mec=col, label=rec["name"])
ax.set_yscale("log")
ax.set_xlabel("1000 / T (K$^{-1}$)")
ax.set_ylabel("$P^{sat}$ (bar)")
ax.plot([], [], "k-", label="Peng–Robinson")
ax.plot([], [], "k:", label="van der Waals")
ax.plot([], [], "ko", ms=2.6, mfc="none", label="reference EOS")
ax.legend(frameon=False, loc="upper right", ncol=1)
label(ax, "(a)")

ax = axs[1]
for m, col in (("vdW", C_G), ("PR", C_L)):
    e, M, _ = pf.make_eos("carbon-dioxide", m)
    Tc, Pc, vc = pf.eos_critical_point(e)
    for tr in (0.85, 0.9, 0.95, 1.0):
        it = pf.isotherm(e, tr * Tc, npts=300, V_max_factor=6)
        ax.plot(it["V"] / vc, it["P"] / Pc, color=col, lw=0.9 if tr < 1 else 1.3, alpha=0.9)
        if it["Psat"]:
            ax.plot([it["vL"] / vc, it["vV"] / vc], [it["Psat"] / Pc] * 2, "--", color=col, lw=0.7)
    cur = pf.saturation_curve(e, npts=50)
    ax.plot(np.r_[cur.vL / vc, (cur.vV / vc)[::-1]], np.r_[cur.P / Pc, (cur.P / Pc)[::-1]], color=C_3, lw=1.0)
ax.set_xlim(0.3, 5); ax.set_ylim(-0.5, 1.6)
ax.axhline(0, color="k", lw=0.4)
ax.set_xlabel("$V/V_c$"); ax.set_ylabel("$P/P_c$")
ax.plot([], [], color=C_L, label="PR"); ax.plot([], [], color=C_G, label="vdW")
ax.plot([], [], "--", color="k", lw=0.7, label="equal-area line"); ax.plot([], [], color=C_3, label="two-phase dome")
ax.legend(frameon=False, loc="upper right")
label(ax, "(b)")
it = pf.isotherm(pf.make_eos("carbon-dioxide", "vdW")[0], 0.9 * 304.13, npts=300)
NUM["vdW_area_mismatch"] = it["maxwell_area"]
e9 = pf.make_eos(None, "vdW", Tc=300.0, Pc=50e5, omega=0)[0]
NUM["vdW_Pr_at_Tr09"] = pf.saturation_pressure(e9, 270.0)[0] / 50e5

ax = axs[2]
for cid, col in (("propane", C_V), ("water", C_L)):
    e, M, _ = pf.make_eos(cid, "PR")
    cur = pf.saturation_curve(e, npts=60)
    ax.plot(M / cur.vL, cur.T, color=col, label=compound_db.find(cid)["name"] + ", PR")
    ax.plot(M / cur.vV, cur.T, color=col)
    ax.plot(M / cur.vc, cur.Tc, "o", color=col, ms=3)
for cid, col in (("propane", C_V), ("water", C_L)):
    Tr_, _, rl_, rv_ = ref_of(cid)
    ax.plot(rl_, Tr_, "o", ms=2.6, mfc="none", mec=col); ax.plot(rv_, Tr_, "o", ms=2.6, mfc="none", mec=col)
ax.plot([], [], "ko", ms=2.6, mfc="none", label="reference EOS")
ax.set_xlabel("Density (kg m$^{-3}$)"); ax.set_ylabel("T (K)")
ax.legend(frameon=False, loc="lower center")
label(ax, "(c)")
fig.tight_layout(w_pad=1.2)
save(fig, "fig_pure")
e, M, _ = pf.make_eos("water", "PR")
NUM["water_rhoL_PR_373"] = M / pf.saturation_pressure(e, 373.15)[1]
NUM["water_rhoL_PR_300"] = M / pf.saturation_pressure(e, 300.0)[1]
cur = pf.saturation_curve(e, npts=80)
i = int(np.argmin(abs(cur.T - 373.15)))
NUM["water_dHvap_PR_373"] = float(np.interp(373.15, cur.T, cur.dHvap)) / 1000

# vapor-pressure table
rows = []
for cid in ("water", "methanol", "ethanol", "acetone", "benzene", "n-hexane", "propane", "methane", "nitrogen",
            "carbon-dioxide"):
    rec = compound_db.find(cid)
    row = {"name": rec["name"], "Tb": rec["Tb"]["value"] if rec["Tb"] else None}
    a = rec["antoine"]
    T7 = 0.7 * rec["Tc"]["value"]
    Pa = 10 ** (a["A"] - a["B"] / (T7 + a["C"])) if a and a["Tmin"] - 5 <= T7 <= a["Tmax"] + 5 else None
    for m in ("vdW", "SRK", "PR"):
        e, M, _ = pf.make_eos(cid, m)
        try:
            row["Tb_" + m] = pf.saturation_temperature(e, 101325.0)
        except ValueError:
            row["Tb_" + m] = None
        row["dP_" + m] = 100 * (pf.saturation_pressure(e, T7)[0] / Pa - 1) if Pa else None
    rows.append(row)
NUM["pure_table"] = rows
dev = [abs(r["Tb_" + m] - r["Tb"]) for r in rows for m in ("SRK", "PR") if r["Tb"] and r["Tb_" + m]]
NUM["Tb_maxdev_SRK_PR"] = max(dev)
NUM["Tb_vdW_range"] = [min(r["Tb"] - r["Tb_vdW"] for r in rows if r["Tb"]), max(r["Tb"] - r["Tb_vdW"] for r in rows if r["Tb"])]

# ---------------------------------------------------------------------------
# 2. pure fluid in a field
# ---------------------------------------------------------------------------
e, M, _ = pf.make_eos("propane", "PR")
fld = pf.PureField("gravity", 9.81, 0.0, 20.0, 0.0)
col = pf.column(e, M, 300.0, 10.0e5, fld, npts=201)
NUM["propane_col"] = {"Psat": col.Psat / 1e5, "interface": col.interface, "rhoL": float(col.rho[0]),
                      "rhoV": float(col.rho[-1]), "P_top": float(col.P[-1] / 1e5),
                      "h_est": (10.0e5 - col.Psat) / (col.rho[0] * 9.81)}
fig, axs = plt.subplots(1, 3, figsize=(TWO, 2.2))
ax = axs[0]
liq = np.array([p_ == "liquid" for p_ in col.phase])
ax.plot(col.P[liq] / 1e5, col.position[liq], color=C_L, label="liquid")
ax.plot(col.P[~liq] / 1e5, col.position[~liq], color=C_V, label="vapor")
ax.axhline(col.interface, color=C_G, lw=0.6, ls="--")
ax.axvline(col.Psat / 1e5, color=C_G, lw=0.6, ls=":")
ax.text(col.Psat / 1e5, 19, " $P^{sat}$", color=C_G, fontsize=6.5, va="top")
ax.text(9.93, col.interface + 0.4, "interface", color=C_G, fontsize=6.5)
ax.set_xlabel("P (bar)"); ax.set_ylabel("Height (m)")
ax.text(9.999, 2.3, r"$\rho_L$=%.0f kg m$^{-3}$" % col.rho[0], fontsize=6, ha="right", color=C_L)
ax.text(9.912, 10, r"$\rho_V$=%.1f kg m$^{-3}$" % col.rho[-1], fontsize=6, ha="left", color=C_V)
ax.legend(frameon=False, loc="center right")
label(ax, "(a)")

ax = axs[1]
Ts = np.linspace(200, 365, 40)
for gfac, ls in ((1, ":"), (5, "--"), (10, "-")):
    f_ = pf.PureField("gravity", 9.81 * gfac, 0.0, 100.0, 0.0)
    T_, Ps, P1, P2 = pf.band(e, M, f_, Ts)
    ax.plot(T_, (P2 - Ps) / 1e5, ls, color=C_V, label=f"last vapor, {gfac} g")
    if gfac == 10:
        ax.fill_between(T_, (P1 - Ps) / 1e5, (P2 - Ps) / 1e5, color=C_3, alpha=0.18, lw=0)
        NUM["band_10g_300K"] = float(np.interp(300, T_, (P2 - Ps) / 1e5))
ax.axhline(0, color=C_L, lw=1.2, label="first liquid (all g)")
ax.set_xlabel("T (K)"); ax.set_ylabel("$P_{bottom} - P^{sat}$ (bar)")
ax.legend(frameon=False, loc="upper right")
label(ax, "(b)")

ax = axs[2]
fc = pf.PureField("centrifugal", 1.0, 0.02, 0.15, 0.02)
e2, M2, _ = pf.make_eos("carbon-dioxide", "PR")
T = 280.0
Ps2 = pf.saturation_pressure(e2, T)[0]
rpms = np.linspace(2000, 30000, 60)
for Pin, ls in ((0.90 * Ps2, "-"), (0.97 * Ps2, "--")):
    pos = []
    ln_fs = pf._ln_f(e2, T, Ps2, "vapor")
    ln_fr = pf._ln_f(e2, T, Pin, "vapor")
    for rpm in rpms:
        f_ = pf.PureField("centrifugal", rpm * 2 * math.pi / 60, 0.02, 0.15, 0.02)
        phi_star = f_.phi(0.02) + 8.314462618 * T * (ln_fr - ln_fs) / M2
        pos.append(f_.position_at(phi_star))
    pos = np.array([np.nan if p_ is None else p_ for p_ in pos])
    ax.plot(rpms / 1000, pos * 100, ls, color=C_3, label=f"$P_{{in}}$ = {Pin / 1e5:.1f} bar")
ax.set_xlabel("Speed (1000 rpm)"); ax.set_ylabel("Interface radius (cm)")
ax.set_ylim(2, 15)
ax.legend(frameon=False, loc="upper right", title=f"CO$_2$, {T:.0f} K, PR", title_fontsize=6.5)
label(ax, "(c)")
fig.tight_layout(w_pad=1.2)
save(fig, "fig_field_pure")
NUM["co2_Psat_280"] = Ps2 / 1e5

# hydrostatic and barometric checks
eN, MN, _ = pf.make_eos("nitrogen", "PR")
cN = pf.column(eN, MN, 300.0, 1e5, pf.PureField("gravity", 9.81, 0.0, 1000.0, 0.0), npts=11)
NUM["N2_top"] = cN.P[-1] / 1e5
NUM["N2_barometric"] = math.exp(-MN * 9.81 * 1000 / (8.314462618 * 300))
eW, MW, _ = pf.make_eos("water", "PR")
cW = pf.column(eW, MW, 300.0, 10e5, pf.PureField("gravity", 9.81, 0.0, 50.0, 0.0), npts=11)
NUM["W_top"] = cW.P[-1] / 1e5
NUM["W_hydro"] = (10e5 - cW.rho[0] * 9.81 * 50.0) / 1e5

# ---------------------------------------------------------------------------
# 3. binary VLE
# ---------------------------------------------------------------------------
def diag(comps, model, hold, val, params=None, npts=81):
    req = {"components": comps, "model": model, "hold": hold, hold: val, "npts": npts}
    if params:
        req["params"] = params
    r = web_api.compute(req)
    assert r["ok"], r.get("error")
    return r


fig, axs = plt.subplots(2, 2, figsize=(TWO, 4.3))


def exp_panel(ax, fname, comps, xlabel, legend=False):
    data = read_csv(fname)
    r = diag(comps, "nrtl", "P", float(data[0]["P_kPa"]) / 100)
    ri = diag(comps, "ideal", "P", float(data[0]["P_kPa"]) / 100)
    c = r["columns"]; ci = ri["columns"]
    ax.plot(c["x1"], c["T"], color=C_L, label="bubble, NRTL"); ax.plot(c["y1"], c["T"], color=C_V, label="dew, NRTL")
    ax.plot(ci["x1"], ci["T"], ":", color=C_L, label="ideal solution"); ax.plot(ci["y1"], ci["T"], ":", color=C_V)
    x = [float(d["x1"]) for d in data]; y = [float(d["y1"]) for d in data]; T = [float(d["T_K"]) for d in data]
    ax.plot(x, T, "o", ms=2.8, mfc="none", mec=C_L, label="$x$, measured")
    ax.plot(y, T, "^", ms=2.8, mfc="none", mec=C_V, label="$y$, measured")
    ax.set_xlabel(xlabel); ax.set_ylabel("T (K)")
    if legend:
        ax.legend(frameon=False, loc="upper right")
    return r


r = exp_panel(axs[0, 0], "vle_ethanol_water_101kPa_Kamihama2012.csv", ["ethanol", "water"], "$x_1, y_1$ ethanol", True)
label(axs[0, 0], "(a)")
NUM["etoh_az"] = r["azeotrope"]
r = exp_panel(axs[0, 1], "vle_chloroform_methanol_101kPa_Li2014.csv", ["chloroform", "methanol"], "$x_1, y_1$ chloroform")
label(axs[0, 1], "(b)")
NUM["chl_az"] = r["azeotrope"]
r = exp_panel(axs[1, 0], "vle_acetone_ethanol_101kPa_KuTu2005.csv", ["acetone", "ethanol"], "$x_1, y_1$ acetone")
label(axs[1, 0], "(c)")
# CO2 + ethane: E-PPR78 prediction (kept for the text; no figure panel)
r = diag(["carbon-dioxide", "ethane"], "PR", "T", 263.15)
NUM["co2_az"] = r["azeotrope"]; NUM["co2_k"] = compound_db.PPR78Kij("carbon-dioxide", "ethane")(263.15)
r0 = diag(["carbon-dioxide", "ethane"], "PR", "T", 263.15, {"kij": 0.0})
NUM["co2_az_k0"] = r0["azeotrope"]

ax = axs[1, 1]
txt = open(os.path.join(CODE, "examples", "benzene_naphthalene.txt")).read()
system = build_system(txt)
res = es.solve(system)
eu = None
for rr in res:
    d = rr.data
    if isinstance(d, es.Diagram):
        ax.plot(d.x1, d.var, color=C_L); ax.plot(d.y1, d.var, color=C_V)
    elif isinstance(d, es.SLECurve):
        xs = d.x_solid if d.solid == "benzene" else 1 - d.x_solid
        keep = d.T >= 0
        ax.plot(xs, d.T, color=C_3 if d.solid == "benzene" else C_K, lw=1.0)
    elif isinstance(d, dict) and rr.equilibrium == "Eutectic":
        eu = d
xe = eu["x"] if eu["component"] == "benzene" else 1 - eu["x"]
ax.axhline(eu["T"], color=C_G, lw=0.6, ls="--")
# measured: VLE (Chemeo, ~99 kPa) and naphthalene solubility in benzene (molality -> mole fraction)
xv = [0.6213, 0.5224, 0.4129, 0.2909, 0.1542]; Tv = [368.55, 375.75, 385.15, 397.75, 421.95]
ax.plot(xv, Tv, "o", ms=2.8, mfc="none", mec=C_L, label="bubble T, measured")
mol = np.array([3.3133, 3.9031, 4.5226, 5.4333, 6.3668, 7.4257, 8.8264, 11.0421, 13.1887, 16.3235])
Tm = np.array([285.70, 290.30, 294.50, 299.50, 303.70, 307.90, 312.10, 317.40, 321.50, 325.80])
xn = mol / (mol + 1000 / 78.11)
ax.plot(1 - xn, Tm, "s", ms=2.6, mfc="none", mec=C_K, label="naphthalene solubility, measured")
ax.set_xlabel("$x$ benzene"); ax.set_ylabel("T (K)"); ax.set_ylim(240, 500)
ax.legend(frameon=False, loc="upper right")
ax.text(0.03, eu["T"] + 4, "eutectic", fontsize=6.5, color=C_G)
label(ax, "(d)")
fig.tight_layout(h_pad=0.8, w_pad=1.2)
save(fig, "fig_binary")
NUM["eutectic"] = {"x_benzene": xe, "T": eu["T"]}
# compare the solubility curve with the measured points
sle = [rr.data for rr in res if isinstance(rr.data, es.SLECurve) and rr.data.solid == "naphthalene"][0]
xn_calc = np.interp(Tm, sle.T[::-1] if sle.T[0] > sle.T[-1] else sle.T,
                    (sle.x_solid[::-1] if sle.T[0] > sle.T[-1] else sle.x_solid))
NUM["sle_dev"] = {"max_abs_dx": float(np.max(np.abs(xn_calc - xn))), "mean_abs_dx": float(np.mean(np.abs(xn_calc - xn)))}
vle = [rr.data for rr in res if isinstance(rr.data, es.Diagram)][0]
Tcalc = np.interp(xv, vle.x1, vle.var)
NUM["bn_vle_dev"] = {"max_abs_dT": float(np.max(np.abs(Tcalc - np.array(Tv))))}

# azeotropes of every NRTL pair (summary) and the three others in the table
az = {}
for a, b in (("methanol", "acetone"), ("ethanol", "benzene"), ("methanol", "benzene"), ("methanol", "toluene"),
             ("ethanol", "ethyl-acetate"), ("methanol", "chloroform")):
    r = diag([a, b], "nrtl", "P", 1.01325, npts=61)
    az[f"{a}|{b}"] = r["azeotrope"]
NUM["az_more"] = az

# ---------------------------------------------------------------------------
# 4. field on a binary, and the ternary diagram
# ---------------------------------------------------------------------------
fig, axs = plt.subplots(1, 2, figsize=(TWO, 2.6), gridspec_kw={"width_ratios": [1, 1.15]})
ax = axs[0]
r = web_api.compute({"components": ["propane", "n-butane"], "model": "PR", "hold": "field", "T": 320, "P": 25.0,
                     "z1": 0.4, "field": {"kind": "gravity", "g": 9.80665, "height": 500}})
c = r["columns"]
pos = np.array(c["position"]); P = np.array(c["P"], float); x1 = np.array(c["x1"], float)
liq = np.array([p_ == "liquid" for p_ in c["phase"]])
ax.plot(P[liq], pos[liq], color=C_L, label="liquid")
ax.plot(P[~liq], pos[~liq], color=C_V, label="vapor (gas cap)")
ax.set_xlabel("P (bar)"); ax.set_ylabel("Height (m)")
ax3 = ax.twiny(); ax3.plot(x1, pos, color=C_K, lw=0.7, ls="-.")
ax3.set_xlabel("$x_{propane}$ of the local phase, dash-dot", fontsize=6.5, labelpad=2); ax3.tick_params(labelsize=6)
ax.legend(frameon=False, loc="upper right")
label(ax, "(a)", x=-0.12, y=1.0)
NUM["gascap"] = r["contacts"]

ax = axs[1]
ex = [x for x in web_api.MIXTURE_EXAMPLES if x["id"] == "solvents"][0]
m = web_api.mixture(dict(ex, map=True, divisions=14))
pts = m["map"]["points"]
A, B, Cc = np.array([0.5, math.sqrt(3) / 2]), np.array([0.0, 0.0]), np.array([1.0, 0.0])
xy = lambda v: v[0] * A + v[1] * B + v[2] * Cc
Xs = np.array([xy(p_["x"]) for p_ in pts if p_["T"]]); Ts = np.array([p_["T"] for p_ in pts if p_["T"]])
tc = ax.tricontourf(Xs[:, 0], Xs[:, 1], Ts, levels=14, cmap="Purples")
for p_ in pts:
    if p_["y"] and round(p_["x"][0] * 14) % 2 == 0 and round(p_["x"][1] * 14) % 2 == 0:
        p1, p2 = xy(p_["x"]), xy(p_["y"])
        ax.annotate("", xy=p2, xytext=p1, arrowprops=dict(arrowstyle="-|>", lw=0.5, color=C_K, mutation_scale=5))
pz, px, py = xy(m["z"]), xy(m["x"]), xy(m["y"])
ax.plot([px[0], py[0]], [px[1], py[1]], "-", color=C_K, lw=1.3)
ax.plot(*px, "o", color=C_L, ms=5, mec="w"); ax.plot(*py, "o", color=C_V, ms=5, mec="w"); ax.plot(*pz, "s", color=C_K, ms=4, mec="w")
ax.plot(*np.array([A, B, Cc, A]).T, color=C_K, lw=0.8)
ax.text(*A + [0, 0.03], "acetone", ha="center", fontsize=7); ax.text(*B - [0.02, 0.07], "methanol", ha="left", fontsize=7)
ax.text(*Cc - [-0.02, 0.07], "water", ha="right", fontsize=7)
ax.set_aspect("equal"); ax.axis("off")
cb = fig.colorbar(tc, ax=ax, shrink=0.75, pad=0.02); cb.set_label("bubble T at 1 atm (K)", fontsize=7); cb.ax.tick_params(labelsize=6)
label(ax, "(b)")
fig.tight_layout(w_pad=1.0)
save(fig, "fig_field_ternary")
NUM["mix"] = {}
for ex in web_api.MIXTURE_EXAMPLES:
    mm = web_api.mixture(ex)
    zz, xx, yy, bb = map(np.array, (mm["z"], mm["x"], mm["y"], [mm["beta"]]))
    NUM["mix"][ex["id"]] = {"label": ex["label"], "names": mm["names"], "T": ex["T"], "P": ex["P"], "z": mm["z"],
                            "x": mm["x"], "y": mm["y"], "K": mm["K"], "beta": mm["beta"], "bounds": mm["bounds"],
                            "balance": float(np.max(np.abs((1 - bb) * xx + bb * yy - zz)))}

# k_ij table
ktab = []
for a, b in (("carbon-dioxide", "ethane"), ("carbon-dioxide", "methane"), ("methane", "ethane"),
             ("propane", "n-butane"), ("carbon-dioxide", "n-hexane"), ("methane", "nitrogen"),
             ("benzene", "cyclohexane"), ("carbon-dioxide", "hydrogen-sulfide")):
    k = compound_db.PPR78Kij(a, b)
    ktab.append({"a": compound_db.find(a)["name"], "b": compound_db.find(b)["name"],
                 "k250": k(250.0), "k298": k(298.15), "k350": k(350.0),
                 "dechema": compound_db.pair(a, b).get("kij_pr")})
NUM["ktab"] = ktab

# checker verdicts
NUM["verdicts"] = []
for fn in sorted(os.listdir(os.path.join(CODE, "examples"))):
    if fn.endswith(".txt"):
        NUM["verdicts"].append([fn, check_input(open(os.path.join(CODE, "examples", fn)).read()).status])

# ---------------------------------------------------------------------------
# 5. validation against measured and reference data (scripts/validate_real_data.py)
# ---------------------------------------------------------------------------
import importlib.util                                  # noqa: E402
spec = importlib.util.spec_from_file_location("vrd", os.path.join(CODE, "scripts", "validate_real_data.py"))
vrd = importlib.util.module_from_spec(spec); spec.loader.exec_module(vrd)
NUM["val_pure"] = vrd.pure_saturation()
NUM["val_act"] = vrd.binary_activity()
NUM["val_eos"] = vrd.binary_eos()
NUM["val_field"] = vrd.field()

fig, axs = plt.subplots(1, 3, figsize=(TWO, 2.3))
ax = axs[0]
fr = read_csv("reference_field.csv")
for gf, mk in ((1, "o"), (5, "s"), (10, "^")):
    rr = [r_ for r_ in fr if r_["case"] == "band_propane_100m" and float(r_["factor"]) == gf]
    ax.plot([float(r_["T_K"]) for r_ in rr], [float(r_["value"]) / 1e5 for r_ in rr], mk, ms=3, mfc="none",
            mec=C_K, label=f"reference EOS, {gf} g")
    for model, col in (("PR", C_L), ("SRK", C_V)):
        e_, M_, _ = pf.make_eos("propane", model)
        T_, Ps_, P1_, P2_ = pf.band(e_, M_, pf.PureField("gravity", 9.81 * gf, 0.0, 100.0, 0.0), np.linspace(200, 366, 30))
        ax.plot(T_, (P2_ - Ps_) / 1e5, color=col, lw=0.9, label=model if gf == 10 else None)
ax.set_xlabel("T (K)"); ax.set_ylabel("band width, $P_{bottom}-P^{sat}$ (bar)")
ax.set_ylim(0, 88)
ax.legend(frameon=False, loc="upper center", fontsize=5.3, ncol=2)
label(ax, "(a)")

ax = axs[1]
me = read_csv("vle_methane_ethane_Han2012.csv")
for T0, col in ((126.01, C_L), (133.98, C_3), (140.01, C_V)):
    dd = [d for d in me if abs(float(d["T_K"]) - T0) < 0.005]
    ax.plot([float(d["x1"]) for d in dd], [float(d["P_kPa"]) / 100 for d in dd], "o", ms=2.8, mfc="none", mec=col)
    ax.plot([float(d["y1"]) for d in dd], [float(d["P_kPa"]) / 100 for d in dd], "^", ms=2.8, mfc="none", mec=col)
    r = diag(["methane", "ethane"], "PR", "T", T0)
    c = r["columns"]
    ax.plot(c["x1"], c["P"], color=col, lw=0.9, label=f"{T0:.0f} K")
    ax.plot(c["y1"], c["P"], color=col, lw=0.9)
ax.set_xlim(0.2, 1.0); ax.set_ylim(0, 7)
ax.set_xlabel("$x_1, y_1$ methane"); ax.set_ylabel("P (bar)")
ax.legend(frameon=False, loc="upper left", title="PR, E-PPR78; symbols: data", title_fontsize=5.5, fontsize=5.5)
label(ax, "(b)")

ax = axs[2]
pts = read_csv("results/pure_saturation_points.csv")
for cid, col, mk in (("water", C_L, "o"), ("methanol", C_3, "s"), ("acetone", C_V, "^"),
                     ("propane", C_K, "D"), ("carbon-dioxide", C_G, "v")):
    rr = [p_ for p_ in pts if p_["compound"] == cid and p_["model"] == "PR"]
    Tc_ = compound_db.find(cid)["Tc"]["value"]
    ax.plot([float(p_["T_K"]) / Tc_ for p_ in rr], [float(p_["dev_rhoL_pct"]) for p_ in rr], mk + "-", ms=2.6,
            lw=0.7, color=col, mfc="none", label=compound_db.find(cid)["name"])
ax.axhline(0, color="k", lw=0.4)
ax.set_xlabel("$T/T_c$"); ax.set_ylabel(r"PR $\rho^L$ deviation (%)")
ax.set_ylim(-30, 14)
ax.legend(frameon=False, loc="center", bbox_to_anchor=(0.42, 0.62), fontsize=5.0, ncol=2)
label(ax, "(c)")
fig.tight_layout(w_pad=1.0)
save(fig, "fig_validation")

# table-of-contents graphic (3.25 x 1.75 in): the coexistence curve widening into a band
fig, ax = plt.subplots(figsize=(3.25, 1.75))
e_, M_, _ = pf.make_eos("propane", "PR")
Tg = np.linspace(200, 368, 60)
for gf, al in ((10, 0.18), (5, 0.28), (1, 0.6)):
    T_, Ps_, P1_, P2_ = pf.band(e_, M_, pf.PureField("gravity", 9.81 * gf, 0.0, 100.0, 0.0), Tg)
    ax.fill_between(T_, Ps_ / 1e5, P2_ / 1e5, color=C_L, alpha=al, lw=0, label=f"{gf} g")
ax.plot(T_, Ps_ / 1e5, color=C_K, lw=1.3, label="no field: $P^{sat}(T)$")
ax.set_xlabel("T (K)"); ax.set_ylabel("P at bottom (bar)")
ax.set_title("Propane, 100 m column: two-phase band", fontsize=7)
ax.legend(frameon=False, fontsize=6, loc="upper left")
save(fig, "toc")

with open(os.path.join(HERE, "numbers.json"), "w") as fh:
    json.dump(NUM, fh, indent=1, default=float)
print(json.dumps({k: v for k, v in NUM.items() if k not in ("pure_table", "mix", "ktab", "verdicts")},
                 indent=1, default=float))
