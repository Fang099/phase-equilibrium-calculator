"""
plotter.py  --  Stage 4 of the Phase Equilibrium Calculator

Draws the results of equilibrium_solver.py with matplotlib:

    Pxy / Txy diagrams       bubble and dew curves, two-phase region, azeotrope,
                             and a marked state point or tie line
    xy diagram               vapor against liquid composition
    complete T-x diagram     vapor-liquid curves, solid-liquid curves and the
                             eutectic together, at one pressure
    composition chart        liquid and vapor compositions for any number of
                             components (flash, bubble or dew point)
    field profiles           pressure and composition against height or radius,
                             with gas-liquid contacts
    3D P-T-x surface         bubble and dew surfaces of a binary
    3D Gibbs-energy surface  Gibbs energy of mixing of the liquid against x and T

Everything is chosen from the input automatically:

    python plotter.py examples/benzene_naphthalene.txt            # saves PNGs to plots/
    python plotter.py examples/ethanol_water_vanlaar.txt --show   # also opens windows
    python plotter.py examples/propane_butane_pr_pxy.txt --no-3d  # skip the 3D surfaces

From Python:

    from thermo_models import build_system
    from equilibrium_solver import solve
    from plotter import plot_results
    system = build_system(open("examples/benzene_naphthalene.txt").read())
    figures = plot_results(system, solve(system))      # {name: Figure}
"""

import argparse
import math
import os
import re
import sys

import numpy as np
import matplotlib

if "--show" not in sys.argv:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt                        # noqa: E402
from mpl_toolkits.mplot3d import Axes3D                # noqa: E402,F401

from thermo_models import R, SolutionLiquid, build_system           # noqa: E402
from equilibrium_solver import (Diagram, SLECurve, VLEPoint, Profile,  # noqa: E402
                                SolverError, saturation, solve, vle_diagram)

# --------------------------------------------------------------------------
# Style: one palette for every figure. The three series colors were checked
# for color-vision-deficiency separation; text never uses them.
# --------------------------------------------------------------------------
COLORS = {"liquid": "#2a78d6", "vapor": "#eb6834", "third": "#1baf7a",
          "fourth": "#6250d6", "ink": "#13233a", "ink2": "#46566a",
          "muted": "#76828f", "grid": "#dbe2e8", "axis": "#b9c3cc",
          "surface": "#fcfcfb"}
SERIES = [COLORS["liquid"], COLORS["vapor"], COLORS["third"], COLORS["fourth"]]

plt.rcParams.update({
    "figure.facecolor": COLORS["surface"], "axes.facecolor": COLORS["surface"],
    "savefig.facecolor": COLORS["surface"], "figure.dpi": 110, "savefig.dpi": 160,
    "font.size": 10.5, "axes.titlesize": 12, "axes.titleweight": "bold",
    "axes.titlelocation": "left", "axes.titlepad": 10, "axes.labelsize": 10.5,
    "axes.labelcolor": COLORS["ink2"], "axes.edgecolor": COLORS["axis"],
    "axes.linewidth": 0.8, "axes.grid": True, "grid.color": COLORS["grid"],
    "grid.linewidth": 0.8, "axes.axisbelow": True, "axes.spines.top": False,
    "axes.spines.right": False, "xtick.color": COLORS["muted"],
    "ytick.color": COLORS["muted"], "xtick.labelcolor": COLORS["ink2"],
    "ytick.labelcolor": COLORS["ink2"], "lines.linewidth": 2.0,
    "lines.solid_capstyle": "round", "legend.frameon": False,
    "legend.fontsize": 9.5, "text.color": COLORS["ink"],
})


def _label(ax, x, y, text, **kw):
    kw.setdefault("ha", "center")
    kw.setdefault("va", "center")
    ax.text(x, y, text, color=COLORS["ink2"], fontsize=9.5, **kw)


def _interp_inverse(xs, vs, x):
    """Value of a curve given as (xs, vs) at composition x (xs need not be sorted)."""
    order = np.argsort(xs)
    return float(np.interp(x, np.asarray(xs)[order], np.asarray(vs)[order]))


def _scale(diagram):
    """Display units: bar for pressure, K for temperature."""
    if diagram.kind == "Pxy":
        return diagram.var / 1e5, "Pressure (bar)", "bar"
    return diagram.var, "Temperature (K)", "K"


# --------------------------------------------------------------------------
# Vapor-liquid diagrams
# --------------------------------------------------------------------------
def plot_vle_diagram(d, ax=None, point=None, title=None):
    """Pxy or Txy diagram with the two-phase region shaded.
    point: optional VLEPoint to mark (bubble/dew point or flash tie line)."""
    if ax is None:
        _, ax = plt.subplots(figsize=(7, 4.8))
    v, ylabel, unit = _scale(d)
    c1 = d.names[0]
    # False when the vapor is a single pure component (a nonvolatile partner):
    # the dew "curve" is then the vertical line y1 = 0 or 1
    two_curves = not np.allclose(d.y1, d.y1[0])
    # two-phase region: polygon bubble curve -> dew curve reversed
    poly_x = np.concatenate([d.x1, d.y1[::-1]])
    poly_v = np.concatenate([v, v[::-1]])
    ax.fill(poly_x, poly_v, color=COLORS["liquid"], alpha=0.10, lw=0, label="Two-phase region")
    ax.plot(d.x1, v, color=COLORS["liquid"], label="Bubble curve (liquid)")
    ax.plot(d.y1, v, color=COLORS["vapor"], label="Dew curve (vapor)")

    lo, hi = float(np.nanmin(v)), float(np.nanmax(v))
    if point is not None:                     # keep the marked state in view
        pv = point.P / 1e5 if d.kind == "Pxy" else point.T
        lo, hi = min(lo, pv), max(hi, pv)
    span = hi - lo if hi > lo else 1.0
    ax.set_ylim(lo - 0.08 * span, hi + 0.12 * span)
    ax.set_xlim(0, 1)
    up_is_liquid = d.kind == "Pxy"           # liquid at high P, or at low T
    ylo, yhi = lo - 0.08 * span, hi + 0.12 * span
    xs = np.linspace(0.1, 0.9, 81)
    vb = np.array([_interp_inverse(d.x1, v, x) for x in xs])
    vd = np.array([_interp_inverse(d.y1, v, x) for x in xs]) if two_curves else vb
    away = np.ones_like(xs, bool)
    if d.azeotrope:
        away = np.abs(xs - d.azeotrope[0]) > 0.2
    # liquid side of the bubble curve, vapor side of the dew curve
    liq_room = (yhi - vb) if up_is_liquid else (vb - ylo)
    vap_room = (vd - ylo) if up_is_liquid else (yhi - vd)
    for room, curve, text, sign in ((liq_room, vb, "Liquid", 1 if up_is_liquid else -1),
                                    (vap_room, vd, "Vapor", -1 if up_is_liquid else 1)):
        r_ = np.where(away, room, -np.inf)
        k = int(np.argmax(r_))
        if r_[k] > 0.1 * span:
            _label(ax, xs[k], curve[k] + sign * min(0.5 * r_[k], 0.2 * span), text)
    if two_curves:
        gap = np.where(away, np.abs(vb - vd), -np.inf)
        k = int(np.argmax(gap))
        if gap[k] > 0.12 * span:
            _label(ax, xs[k], 0.5 * (vb[k] + vd[k]), "Liquid + vapor")

    if d.azeotrope:
        xa, va = d.azeotrope
        va = va / 1e5 if d.kind == "Pxy" else va
        ax.plot([xa], [va], "o", ms=7, color=COLORS["ink2"], mec=COLORS["surface"], mew=2, zorder=5)
        toward_middle = 1 if va < 0.5 * (lo + hi) else -1
        ax.annotate(f"Azeotrope\nx = {xa:.3f}, {va:.4g} {unit}", (xa, va),
                    xytext=(-14, 30 * toward_middle), textcoords="offset points",
                    ha="right", va="bottom" if toward_middle > 0 else "top",
                    fontsize=9, color=COLORS["ink2"],
                    arrowprops=dict(arrowstyle="-", color=COLORS["muted"], lw=0.8))
    if point is not None:
        _mark_point(ax, d, point)

    ax.set_xlabel(f"Mole fraction of {c1}")
    ax.set_ylabel(ylabel)
    fixed = (f"T = {d.fixed:.2f} K" if d.kind == "Pxy" else f"P = {d.fixed / 1e5:.4g} bar")
    ax.set_title(title or f"{d.kind} diagram, {d.names[0]} / {d.names[1]}, {fixed}")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=3)
    return ax


def _mark_point(ax, d, pt):
    """Mark a solved state: a tie line from liquid to vapor at its T or P."""
    val = pt.P / 1e5 if d.kind == "Pxy" else pt.T
    if pt.phase == "two-phase":
        x, y = pt.x[0], pt.y[0]
        if np.isfinite(x) and np.isfinite(y):
            ax.plot([x, y], [val, val], color=COLORS["ink2"], lw=1.2, zorder=4)
            ax.plot([x], [val], "o", ms=7, color=COLORS["liquid"], mec=COLORS["surface"], mew=2, zorder=5)
            ax.plot([y], [val], "o", ms=7, color=COLORS["vapor"], mec=COLORS["surface"], mew=2, zorder=5)
            ax.annotate(f"x = {x:.3f}", (x, val), xytext=(0, 10), textcoords="offset points",
                        ha="center", fontsize=9, color=COLORS["ink2"])
            ax.annotate(f"y = {y:.3f}", (y, val), xytext=(0, -16), textcoords="offset points",
                        ha="center", fontsize=9, color=COLORS["ink2"])
    else:
        z = pt.x[0] if np.isfinite(pt.x[0]) else pt.y[0]
        ax.plot([z], [val], "o", ms=7, color=COLORS["ink2"], mec=COLORS["surface"], mew=2, zorder=5)
        ax.annotate(f"{pt.phase.capitalize()} only at this state", (z, val), xytext=(8, 8),
                    textcoords="offset points",
                    fontsize=9, color=COLORS["ink2"])


def plot_xy(d, ax=None):
    """Vapor mole fraction against liquid mole fraction, with the y = x line."""
    if ax is None:
        _, ax = plt.subplots(figsize=(4.8, 4.8))
    ax.plot([0, 1], [0, 1], color=COLORS["axis"], lw=1)
    ax.plot(d.x1, d.y1, color=COLORS["liquid"])
    if d.azeotrope:
        ax.plot([d.azeotrope[0]], [d.azeotrope[0]], "o", ms=7, color=COLORS["ink2"],
                mec=COLORS["surface"], mew=2)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.set_xlabel(f"x, liquid mole fraction of {d.names[0]}")
    ax.set_ylabel(f"y, vapor mole fraction of {d.names[0]}")
    fixed = (f"T = {d.fixed:.2f} K" if d.kind == "Pxy" else f"P = {d.fixed / 1e5:.4g} bar")
    ax.set_title(f"xy diagram, {fixed}")
    return ax


# --------------------------------------------------------------------------
# Complete T-x diagram (vapor-liquid + solid-liquid) at one pressure
# --------------------------------------------------------------------------
def _sle_x1(curve, names):
    """Solid-liquid curve in terms of the mole fraction of component 1."""
    xs = np.asarray(curve.x_solid)
    return (xs if names.index(curve.solid) == 0 else 1 - xs), np.asarray(curve.T)


def plot_tx_complete(names, P, txy=None, sle=(), eut=None, ax=None):
    """Phase diagram at fixed P: boiling curves, the liquidus of each solid,
    and the eutectic line below which both solids coexist."""
    if ax is None:
        _, ax = plt.subplots(figsize=(7.4, 5.2))
    Tall = []
    if txy is not None:
        poly_x = np.concatenate([txy.x1, txy.y1[::-1]])
        poly_T = np.concatenate([txy.var, txy.var[::-1]])
        ax.fill(poly_x, poly_T, color=COLORS["vapor"], alpha=0.10, lw=0)
        ax.plot(txy.x1, txy.var, color=COLORS["vapor"], label="Boiling (bubble) curve")
        if not np.allclose(txy.y1, txy.y1[0]):
            ax.plot(txy.y1, txy.var, color=COLORS["vapor"], lw=1.3, alpha=0.8,
                    label="Dew curve")
        Tall += list(txy.var)
    ex1 = None
    if eut is not None:
        ex1 = eut["x"] if names.index(eut["component"]) == 0 else 1 - eut["x"]
    sle_colors = [COLORS["liquid"], COLORS["third"]]
    for k, c in enumerate(sle):
        x1, T = _sle_x1(c, names)
        if ex1 is not None:                       # liquidus ends at the eutectic
            keep = (x1 >= ex1 - 1e-12) if names.index(c.solid) == 0 else (x1 <= ex1 + 1e-12)
            x1, T = x1[keep], T[keep]
            x1 = np.append(x1, ex1)
            T = np.append(T, eut["T"])
        ax.plot(x1, T, color=sle_colors[k % 2], label=f"{c.solid[0].upper() + c.solid[1:]} crystallizes")
        Tall += list(T)
    if eut is not None:
        ax.axhline(eut["T"], color=COLORS["ink2"], lw=1.1)
        ax.plot([ex1], [eut["T"]], "o", ms=7, color=COLORS["ink2"], mec=COLORS["surface"],
                mew=2, zorder=5)
        right = ex1 > 0.5                   # label on the side with more room
        ax.annotate(f"Eutectic: x = {ex1:.3f}\n{eut['T']:.1f} K", (ex1, eut["T"]),
                    xytext=(0, -8), textcoords="offset points",
                    ha="right" if right else "left", va="top", fontsize=9,
                    color=COLORS["ink2"])
        Tall.append(eut["T"])
    lo, hi = min(Tall), max(Tall)
    span = hi - lo
    ax.set_ylim(lo - (0.2 if eut is not None else 0.08) * span, hi + 0.08 * span)
    ax.set_xlim(0, 1)
    if eut is not None:
        _label(ax, ex1 / 2 if ex1 > 0.5 else (1 + ex1) / 2, eut["T"] - 0.08 * span,
               f"Solid {names[0]} + solid {names[1]}")
        # liquid + solid regions: between each liquidus and the eutectic line,
        # labelled where that region is tallest
        for c in sle:
            x1, T = sorted_xy(_sle_x1(c, names))
            own = (x1 >= ex1) if names.index(c.solid) == 0 else (x1 <= ex1)
            if own.sum() < 2:
                continue
            h = T[own] - eut["T"]
            j = int(np.argmax(h))
            if h[j] > 0.12 * span:
                xs_ = min(max(x1[own][j], 0.2), 0.8)
                Tl = float(np.interp(xs_, x1, T))
                if Tl - eut["T"] > 0.12 * span:
                    _label(ax, xs_, 0.5 * (Tl + eut["T"]), f"Liquid + solid {c.solid}")
    # liquid label between the highest liquidus and the boiling curve at mid composition
    xm = 0.5 if ex1 is None else 0.5 * (ex1 + (1 if ex1 < 0.5 else 0))
    top = _interp_inverse(txy.x1, txy.var, xm) if txy is not None else hi
    bottom = max([float(np.interp(xm, *sorted_xy(_sle_x1(c, names)))) for c in sle] or [lo])
    if top > bottom:
        _label(ax, xm, 0.5 * (top + bottom), "Liquid")
    if txy is not None:
        _label(ax, 0.5, hi + 0.03 * span, "Vapor")
    ax.set_xlabel(f"Mole fraction of {names[0]}")
    ax.set_ylabel("Temperature (K)")
    ax.set_title(f"Phase diagram of {names[0]} / {names[1]} at P = {P / 1e5:.4g} bar")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=4)
    return ax


def sorted_xy(xy):
    x, y = xy
    o = np.argsort(x)
    return x[o], y[o]


# --------------------------------------------------------------------------
# Compositions of a single state (any number of components)
# --------------------------------------------------------------------------
def plot_compositions(pt, names, ax=None, title=None):
    if ax is None:
        _, ax = plt.subplots(figsize=(max(5, 1.2 * len(names) + 2.5), 4))
    idx = np.arange(len(names))
    sets = [(pt.x, COLORS["liquid"], "Liquid x"), (pt.y, COLORS["vapor"], "Vapor y")]
    sets = [s_ for s_ in sets if np.all(np.isfinite(s_[0]))]
    w = 0.3
    bars = []
    for k, (vals, color, lab) in enumerate(sets):
        off = (k - (len(sets) - 1) / 2) * w
        bars.append(ax.bar(idx + off, vals, w * 0.9, color=color, label=lab))
    for b in bars:
        ax.bar_label(b, fmt="%.3f", padding=3, fontsize=8.5, color=COLORS["ink2"])
    ax.set_xticks(idx, names)
    ax.set_ylim(0, 1.08)
    ax.set_ylabel("Mole fraction")
    ax.grid(axis="x", visible=False)
    extra = f", vapor fraction {pt.beta:.3f}" if pt.phase == "two-phase" else f" ({pt.phase})"
    ax.set_title(title or f"T = {pt.T:.2f} K, P = {pt.P / 1e5:.4g} bar{extra}")
    ax.legend(loc="upper right")
    return ax


# --------------------------------------------------------------------------
# Field profiles
# --------------------------------------------------------------------------
def plot_profile(prof, comp=0, title=None):
    """Pressure and composition against position. Height is drawn vertically
    for gravity, radius horizontally for a centrifuge."""
    vertical = prof.kind == "gravity"
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.6) if not vertical else (9.4, 5.4))
    pos = prof.position
    phases = np.array(prof.phase)
    contacts = [float(m.group(1)) for s in prof.notes
                for m in [re.search(r"contact at ([-\d.e+]+) m", s)] if m]
    pos_label = "Height (m)" if vertical else "Radius (m)"
    quantities = [(prof.P / 1e5, "Pressure (bar)"),
                  (prof.x[:, comp], f"Mole fraction of {prof.names[comp]}")]
    cinfo = prof.extra.get("contacts", [])
    for q, (ax, (val, vlabel)) in enumerate(zip(axes, quantities)):
        for ph, color in (("liquid", COLORS["liquid"]), ("vapor", COLORS["vapor"])):
            m = phases == ph
            if not m.any():
                continue
            pp, vv = list(pos[m]), list(val[m])
            for c in cinfo:                    # both phases meet at the contact
                pp.append(c["position"])
                vv.append(c["P"] / 1e5 if q == 0 else c["x_" + ph][comp])
            o = np.argsort(pp)
            pp, vv = np.asarray(pp)[o], np.asarray(vv)[o]
            if vertical:
                ax.plot(vv, pp, color=color, label=ph.capitalize())
            else:
                ax.plot(pp, vv, color=color, label=ph.capitalize())
        for zc in contacts:
            if vertical:
                ax.axhline(zc, color=COLORS["ink2"], lw=1)
                ax.annotate(f"gas-liquid contact, {zc:.4g} m", (1, zc), xycoords=("axes fraction", "data"),
                            xytext=(-4, -5), textcoords="offset points", ha="right", va="top",
                            fontsize=8.5, color=COLORS["ink2"])
            else:
                ax.axvline(zc, color=COLORS["ink2"], lw=1)
        if vertical:
            ax.set_xlabel(vlabel)
            ax.set_ylabel(pos_label)
        else:
            ax.set_xlabel(pos_label)
            ax.set_ylabel(vlabel)
    axes[0].legend(loc="best")
    T = prof.extra.get("T")
    lab = prof.extra.get("label", "")
    fig.suptitle(title or f"{prof.kind.capitalize()} field profile"
                 + (f", T = {T:.2f} K" if T else "") + (f", {lab}" if lab else ""),
                 x=0.01, ha="left", fontsize=12, fontweight="bold", color=COLORS["ink"])
    fig.tight_layout()
    return fig


# --------------------------------------------------------------------------
# 3D surfaces
# --------------------------------------------------------------------------
def ptx_surface(L, V, T_values, nx=21):
    """Bubble pressure P(T, x) and the matching vapor composition y(T, x).
    Points that cannot be solved (e.g. above a critical point) are NaN."""
    xs = np.linspace(0, 1, nx)
    Pb = np.full((len(T_values), nx), np.nan)
    Y = np.full_like(Pb, np.nan)
    for i, T in enumerate(T_values):
        guess = None
        for j, x1 in enumerate(xs[::-1]):
            jj = nx - 1 - j
            try:
                pt = saturation(L, V, [x1, 1 - x1], 0.0, T=T, guess=guess)
            except (SolverError, ValueError, ZeroDivisionError, OverflowError):
                continue
            guess = pt.P
            Pb[i, jj], Y[i, jj] = pt.P, pt.y[0]
    X, TT = np.meshgrid(xs, T_values)
    return X, Y, TT, Pb


def _auto_T_range(d, L, V):
    """Temperatures for the P-T-x surface: around the fixed T of a Pxy
    diagram, or spanning the boiling points of a Txy diagram."""
    if d.kind == "Pxy":
        T0 = d.fixed
        return np.linspace(T0 - 0.12 * T0, T0 + 0.08 * T0, 11)
    Tlo, Thi = float(np.nanmin(d.var)), float(np.nanmax(d.var))
    pad = 0.08 * (Thi - Tlo) + 5
    return np.linspace(Tlo - pad, Thi + pad, 11)


def plot_ptx_surface(L, V, T_values, names, nx=21, ax=None):
    X, Y, TT, Pb = ptx_surface(L, V, T_values, nx)
    if ax is None:
        fig = plt.figure(figsize=(8, 6.4))
        ax = fig.add_subplot(projection="3d")
    ax.plot_surface(X, TT, Pb / 1e5, color=COLORS["liquid"], alpha=0.45, lw=0,
                    antialiased=True, shade=True)
    ax.plot_wireframe(Y, TT, Pb / 1e5, color=COLORS["vapor"], lw=0.6, alpha=0.7,
                      rstride=0, cstride=2)
    for i in range(0, len(T_values), max(1, len(T_values) // 5)):   # isotherms = Pxy curves
        ax.plot(X[i], TT[i], Pb[i] / 1e5, color=COLORS["liquid"], lw=1.6)
        ax.plot(Y[i], TT[i], Pb[i] / 1e5, color=COLORS["vapor"], lw=1.6)
    ax.set_xlabel(f"Mole fraction of {names[0]}", labelpad=8)
    ax.set_ylabel("Temperature (K)", labelpad=8)
    ax.set_zlabel("Pressure (bar)", labelpad=6)
    ax.set_title(f"Bubble surface (blue) and dew surface (orange), {names[0]} / {names[1]}\n"
                 f"bold lines are isotherms: each pair is one Pxy diagram")
    ax.view_init(elev=16, azim=-72)
    _style_3d(ax)
    return ax


def gmix_surface(L, T_values, P, nx=41):
    """Gibbs energy of mixing of the liquid, in kJ/mol, on an (x, T) grid:
    dG_mix = RT sum_i x_i [ln f_i(x) - ln f_i(pure i)]."""
    xs = np.linspace(1e-4, 1 - 1e-4, nx)
    G = np.zeros((len(T_values), nx))
    for i, T in enumerate(T_values):
        pure = [L.ln_f(T, P, np.eye(2)[k])[k] for k in range(2)]
        for j, x1 in enumerate(xs):
            x = np.array([x1, 1 - x1])
            G[i, j] = R * T * float(np.dot(x, L.ln_f(T, P, x) - pure)) / 1e3
    X, TT = np.meshgrid(xs, T_values)
    return X, TT, G


def plot_gmix_surface(L, T_values, P, names, ax=None):
    X, TT, G = gmix_surface(L, T_values, P)
    if ax is None:
        fig = plt.figure(figsize=(8, 6.4))
        ax = fig.add_subplot(projection="3d")
    ax.plot_surface(X, TT, G, cmap="Blues_r", alpha=0.9, lw=0, antialiased=True)
    for i in range(0, len(T_values), max(1, len(T_values) // 5)):
        ax.plot(X[i], TT[i], G[i], color=COLORS["ink2"], lw=0.8)
    ax.set_xlabel(f"Mole fraction of {names[0]}", labelpad=8)
    ax.set_ylabel("Temperature (K)", labelpad=8)
    ax.set_zlabel("Gibbs energy of mixing (kJ/mol)", labelpad=8)
    model = getattr(getattr(L, "activity", None), "name", "ideal solution")
    ax.set_title(f"Liquid Gibbs energy of mixing ({model}), P = {P / 1e5:.4g} bar")
    ax.view_init(elev=24, azim=-52)
    _style_3d(ax)
    return ax


def _style_3d(ax):
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.set_facecolor(COLORS["surface"])
        axis.pane.set_edgecolor(COLORS["grid"])
        axis._axinfo["grid"]["color"] = COLORS["grid"]
        axis._axinfo["grid"]["linewidth"] = 0.6


# --------------------------------------------------------------------------
# Choose and draw every relevant figure for a solved system
# --------------------------------------------------------------------------
def plot_results(system, results, surfaces=True):
    """Return {name: Figure} for everything the results contain."""
    figs = {}
    names = system.names
    st = system.state
    binary = len(names) == 2
    txy, sle, eut = None, [], None

    for r in results:
        d = r.data
        if isinstance(d, Diagram):
            L, V = r.phases
            fig, ax = plt.subplots(figsize=(7, 4.8))
            plot_vle_diagram(d, ax)
            fig.tight_layout()
            figs[f"{d.kind.lower()}"] = fig
            fig2, ax2 = plt.subplots(figsize=(4.8, 4.8))
            plot_xy(d, ax2)
            fig2.tight_layout()
            figs["xy"] = fig2
            if d.kind == "Txy":
                txy = d
            if surfaces:
                Ts = _auto_T_range(d, L, V)
                fig3 = plt.figure(figsize=(8, 6.4))
                plot_ptx_surface(L, V, Ts, names, ax=fig3.add_subplot(projection="3d"))
                fig3.tight_layout()
                figs["ptx_surface"] = fig3
                if isinstance(L, SolutionLiquid):
                    Pg = d.fixed if d.kind == "Txy" else float(np.nanmedian(d.var))
                    fig4 = plt.figure(figsize=(8, 6.4))
                    plot_gmix_surface(L, Ts, Pg, names, ax=fig4.add_subplot(projection="3d"))
                    fig4.tight_layout()
                    figs["gmix_surface"] = fig4
        elif isinstance(d, VLEPoint):
            L, V = r.phases
            if binary:
                try:
                    if st["T"] is not None:        # Pxy at the given temperature
                        dd = vle_diagram(L, V, T=d.T)
                    else:
                        dd = vle_diagram(L, V, P=d.P)
                    fig, ax = plt.subplots(figsize=(7, 4.8))
                    plot_vle_diagram(dd, ax, point=d)
                    fig.tight_layout()
                    figs["state_on_diagram"] = fig
                except SolverError:
                    pass
            if d.phase == "two-phase" or not binary:
                fig, ax = plt.subplots(figsize=(max(5, 1.2 * len(names) + 2.5), 4))
                plot_compositions(d, names, ax)
                fig.tight_layout()
                figs["compositions"] = fig
        elif isinstance(d, SLECurve):
            sle.append(d)
        elif isinstance(d, dict) and "component" in d:
            eut = d
        for k, prof in enumerate(r.profiles):
            figs[f"profile_{k + 1}"] = plot_profile(prof)

    if binary and sle and st["P"] is not None:
        fig, ax = plt.subplots(figsize=(7.4, 5.2))
        plot_tx_complete(names, st["P"], txy, sle, eut, ax)
        fig.tight_layout()
        figs["tx_complete"] = fig
    return figs


def save_figures(figs, outdir, stem):
    os.makedirs(outdir, exist_ok=True)
    paths = []
    for name, fig in figs.items():
        path = os.path.join(outdir, f"{stem}_{name}.png")
        fig.savefig(path, bbox_inches="tight")
        paths.append(path)
    return paths


def main(argv=None):
    ap = argparse.ArgumentParser(description="Solve an input file and plot the results.")
    ap.add_argument("input", help="input file, e.g. examples/benzene_naphthalene.txt")
    ap.add_argument("--out", default="plots", help="folder for the PNG files (default: plots)")
    ap.add_argument("--show", action="store_true", help="open the figures in windows")
    ap.add_argument("--no-3d", action="store_true", help="skip the 3D surfaces")
    a = ap.parse_args(argv)
    from input_checker import check_input
    text = open(a.input).read()
    rep = check_input(text)
    if rep.status != "SUFFICIENT":
        print(rep)
        return 1
    system = build_system(text)
    results = solve(system)
    figs = plot_results(system, results, surfaces=not a.no_3d)
    stem = os.path.splitext(os.path.basename(a.input))[0]
    for p in save_figures(figs, a.out, stem):
        print("saved", p)
    if a.show:
        plt.show()
    return 0


if __name__ == "__main__":
    sys.exit(main())
