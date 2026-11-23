"""Figures from the QC simulations. Draws figures only: every performance measure is computed in R by
R/compute_metrics.R (run it first).

Inputs : tables/scenario_metrics.csv, Table_S1_effort.csv, thresholds.csv (written by R/compute_metrics.R)
Outputs: figures/Figure_1, _2, _3, _4, _5, _6 and Figure_S2 (.png, 300 dpi; see style.py for submission)
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon, Arc, FancyArrowPatch
import style as S

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB = os.path.join(ROOT, "tables")          # tables written by R/compute_metrics.R
FIG = os.path.join(ROOT, "figures")
S.apply()

m = pd.read_csv(os.path.join(TAB, "scenario_metrics.csv"), dtype={"design": str})
eff = pd.read_csv(os.path.join(TAB, "Table_S1_effort.csv"))
eff["design"] = eff["design"].map({v: k for k, v in {
    "65": "65% overlap", "10": "10% overlap", "0": "0% overlap", "200": "200 m", "gimbal": "Gimbal"}.items()})
thr = pd.read_csv(os.path.join(TAB, "thresholds.csv"), dtype={"design": str})

ALT, NADIR_HFOV, NADIR_W = 120.0, 32.2, 34.7
G_HFOV, G_TILT, G_W, G_SIGMA = 25.0, 55.0, 289.7, 160.0
DENS = [50, 100, 200, 300, 400, 500, 600, 700, 800, 900]


def design_handles(designs, open_note=False, flag_note=False):
    h = [Line2D([], [], color=S.DESIGN_COLORS[d], marker=S.DESIGN_MARKERS[d], ms=5, lw=1.5,
                label=S.DESIGN_LABELS[d]) for d in designs]
    if open_note:
        h.append(Line2D([], [], color=S.MUTED, marker="o", mfc="white", ms=5, lw=0,
                        label="Open symbol: >5% failed replicates"))
    if flag_note:
        h.append(Line2D([], [], marker=S.DESIGN_MARKERS["gimbal"], mfc="white", mec=FLAG_EDGE, mew=0.9, ms=5,
                        lw=0, label="Gimbal, < 20 detections (not reliably estimable)"))
    return h


# Gimbal scenarios with < 20 mean detections (compute_metrics.py: reliably_estimable == False)
FLAG_EDGE = "#8C8C8C"
def flagged(s):
    return ~s["reliably_estimable"].astype(bool)
ANY_FLAG = bool(flagged(m).any())


def cap_label(v):
    return f"{v:.0f}" if v >= 10 else f"{v:.1f}"


# ---- Figure 1: nadir camera geometry -------------------------------------------------------
def figure_1():
    fig, ax = plt.subplots(figsize=(4.2, 3.6))
    ax.grid(False)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_xticks([]); ax.set_yticks([])
    ax.add_patch(Polygon([[0, ALT], [-NADIR_W, 0], [NADIR_W, 0]], closed=True, fc="#009E73", alpha=0.15, ec="none"))
    ax.plot([0, -NADIR_W], [ALT, 0], color=S.INK, lw=1)
    ax.plot([0, NADIR_W], [ALT, 0], color=S.INK, lw=1)
    ax.plot([0, 0], [ALT, 0], color=S.MUTED, lw=0.8, ls="--")
    ax.plot([-75, 75], [0, 0], color=S.INK, lw=1.5)
    ax.add_patch(Polygon([[-6, ALT + 2], [6, ALT + 2], [4, ALT - 2], [-4, ALT - 2]], fc=S.INK))
    ax.plot([-10, 10], [ALT + 3, ALT + 3], color=S.INK, lw=1.2)
    r = 22
    ax.add_patch(Arc((0, ALT), 2 * r, 2 * r, theta1=270 - NADIR_HFOV / 2, theta2=270, color=S.INK, lw=0.8))
    ax.add_patch(Arc((0, ALT), 2 * r, 2 * r, theta1=270, theta2=270 + NADIR_HFOV / 2, color=S.INK, lw=0.8))
    ax.text(-9, ALT - r + 2, "16.1°", ha="right", va="center", fontsize=7)
    ax.text(9, ALT - r + 2, "16.1°", ha="left", va="center", fontsize=7)
    ax.text(0, ALT + 9, f"Horizontal FOV {NADIR_HFOV}° (640 × 512 px thermal)", ha="center", fontsize=7.5)
    ax.annotate("", xy=(-62, 0), xytext=(-62, ALT), arrowprops=dict(arrowstyle="<->", color=S.INK, lw=0.8))
    ax.text(-64, ALT / 2, "120 m AGL", rotation=90, ha="right", va="center", fontsize=7.5)
    ax.annotate("", xy=(0, -9), xytext=(NADIR_W, -9), arrowprops=dict(arrowstyle="<->", color=S.INK, lw=0.8))
    ax.text(NADIR_W / 2, -15, "w = 34.7 m", ha="center", va="top", fontsize=7.5)
    ax.annotate("", xy=(-NADIR_W, -27), xytext=(NADIR_W, -27), arrowprops=dict(arrowstyle="<->", color=S.INK, lw=0.8))
    ax.text(0, -33, "Swath 69.4 m (≈10.8 cm per pixel)", ha="center", va="top", fontsize=7.5)
    ax.text(2, 40, "nadir", ha="left", fontsize=6.5, color=S.MUTED, rotation=90, va="bottom")
    ax.set_xlim(-80, 80); ax.set_ylim(-45, ALT + 18); ax.set_aspect("equal")
    S.save(fig, "Figure_1", FIG); plt.close(fig)


# ---- Figure 2: detection functions and gimbal geometry ----------------------------------------
def figure_2():
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.0, 2.9), gridspec_kw=dict(width_ratios=[1.25, 1]))
    x = np.linspace(0, 320, 801)
    a.plot([0, NADIR_W, NADIR_W], [1, 1, 0], color=S.DESIGN_COLORS["0"], lw=2,
           label="Nadir designs: uniform, g(x) = 1 to 34.7 m")
    xs = x[x <= G_W]
    a.plot(xs, np.exp(-xs ** 2 / (2 * G_SIGMA ** 2)), color=S.DESIGN_COLORS["gimbal"], lw=2,
           label="Gimbal: half-normal, σ = 160 m, w = 289.7 m")
    a.plot([G_W, G_W], [0, np.exp(-G_W ** 2 / (2 * G_SIGMA ** 2))], color=S.DESIGN_COLORS["gimbal"], lw=1, ls=":")
    a.text(G_W - 4, np.exp(-G_W ** 2 / (2 * G_SIGMA ** 2)) + 0.05, f"g(w) = {np.exp(-G_W**2/(2*G_SIGMA**2)):.2f}",
           ha="center", fontsize=7)
    a.set_xlim(0, 320); a.set_ylim(0, 1.08)
    a.set_xlabel("Perpendicular distance from transect (m)")
    a.set_ylabel("Detection probability g(x)")
    a.legend(loc="upper right", fontsize=6.8)
    a.set_title("a  Detection functions", loc="left")

    b.grid(False); b.set_xticks([0, 110, 289.7]); b.set_xticklabels(["0", "110", "289.7"])
    b.set_yticks([0, 120]); b.set_yticklabels(["0", "120"])
    inner, outer = G_TILT - G_HFOV / 2, G_TILT + G_HFOV / 2
    xi, xo = ALT * np.tan(np.radians(inner)), ALT * np.tan(np.radians(outer))
    b.add_patch(Polygon([[0, ALT], [0, 0], [xo, 0]], fc=S.DESIGN_COLORS["gimbal"], alpha=0.08, ec="none"))
    b.add_patch(Polygon([[0, ALT], [xi, 0], [xo, 0]], fc=S.DESIGN_COLORS["gimbal"], alpha=0.30, ec="none"))
    b.plot([0, xi], [ALT, 0], color=S.INK, lw=0.8); b.plot([0, xo], [ALT, 0], color=S.INK, lw=0.8)
    b.plot([0, 0], [ALT, 0], color=S.MUTED, lw=0.8, ls="--")
    xa = ALT * np.tan(np.radians(G_TILT))
    b.plot([0, xa], [ALT, 0], color=S.MUTED, lw=0.8, ls=":")
    b.add_patch(Arc((0, ALT), 90, 90, theta1=270, theta2=270 + G_TILT, color=S.INK, lw=0.8))
    b.text(5, 22, "55° max tilt", fontsize=6.8)
    b.annotate("25° FOV at max tilt", xy=(205, 20), xytext=(215, 80), fontsize=6.8, ha="center",
               arrowprops=dict(arrowstyle="-", color=S.MUTED, lw=0.6))
    b.plot([0, 300], [0, 0], color=S.INK, lw=1.2)
    b.add_patch(Polygon([[-7, ALT + 3], [7, ALT + 3], [5, ALT - 2], [-5, ALT - 2]], fc=S.INK, clip_on=False))
    b.set_xlim(-15, 310); b.set_ylim(-5, 145)
    b.set_xlabel("Distance from transect (m)"); b.set_ylabel("Height (m)")
    b.set_title("b  Gimbal geometry (one side)", loc="left")
    b.text(150, 128, "Swept: 0–289.7 m by tilting 0–55°", ha="center", fontsize=6.8, color=S.MUTED)
    fig.tight_layout()
    S.save(fig, "Figure_2", FIG); plt.close(fig)


# ---- Figure 3: survey effort -------------------------------------------------------------------
def figure_3():
    fig, ax = plt.subplots(figsize=(4.6, 3.6))
    for d in S.DESIGN_ORDER:
        e = eff[eff["design"] == d].sort_values("blocks")
        ax.plot(e["blocks"], e["transect_length_km"], color=S.DESIGN_COLORS[d], marker=S.DESIGN_MARKERS[d], ms=5)
    ax.set_yscale("log"); ax.set_xticks(S.BLOCKS)
    ax.set_xlabel("Number of 2 × 2 km blocks (one block per survey day)")
    ax.set_ylabel("Total transect length (km), log scale")
    fig.legend(handles=design_handles(S.DESIGN_ORDER), loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=(0, 0.13, 1, 1))
    S.save(fig, "Figure_3", FIG); plt.close(fig)


def panel_grid(height=5.2):
    fig, axes = plt.subplots(2, 2, figsize=(7.0, height), sharex=True, sharey=True)
    return fig, axes.ravel()


# ---- Figure 4: relative estimated density ------------------------------------------------------
def figure_4(ymax=4.0):
    fig, axes = panel_grid(height=5.6)
    offs = {"0": -0.24, "200": 0.0, "gimbal": 0.24}
    pos = {d: i for i, d in enumerate(DENS)}
    for ax, b in zip(axes, S.BLOCKS):
        ax.axhline(1.0, color=S.INK, lw=0.8, ls="--", zorder=1)
        for d in S.SIMULATED:
            s = m[(m["design"] == d) & (m["blocks"] == b)].sort_values("true_density")
            x = s["true_density"].map(pos) + offs[d]
            lo, hi = s["rel_q025"], s["rel_q975"]
            fl = flagged(s).to_numpy()
            col = np.where(fl, FLAG_EDGE, S.DESIGN_COLORS[d])
            ax.vlines(x, lo, np.minimum(hi, ymax), colors=col, lw=1.2)
            ax.plot(x[~fl], s["rel_median"][~fl], ls="none", marker=S.DESIGN_MARKERS[d], ms=4.5,
                    color=S.DESIGN_COLORS[d], mec="white", mew=0.5, zorder=3)
            ax.plot(x[fl], s["rel_median"][fl], ls="none", marker=S.DESIGN_MARKERS[d], ms=4.5,
                    mfc="white", mec=FLAG_EDGE, mew=0.9, zorder=3)
            for xi, h, c in zip(x, hi, col):
                if h > ymax:
                    ax.plot(xi, ymax, marker="^", ms=4, color=c, clip_on=False, zorder=4)
                    ax.text(xi, ymax + 0.12, cap_label(h), ha="center", va="bottom", fontsize=5.5, color=S.MUTED)
        ax.set_title(f"{b} blocks", loc="left", pad=12)
        ax.set_ylim(0, ymax); ax.set_xlim(-0.6, len(DENS) - 0.4)
        ax.set_xticks(range(len(DENS))); ax.set_xticklabels(DENS)
        ax.grid(axis="x", visible=False)
    for ax in axes[2:]:
        ax.set_xlabel(S.DENS_LABEL)
    for ax in axes[::2]:
        ax.set_ylabel("Estimated / true density")
    h = design_handles(S.SIMULATED, flag_note=ANY_FLAG)
    if (m.loc[m["design"].isin(S.SIMULATED), "rel_q975"] > ymax).any():   # arrows drawn
        h.append(Line2D([], [], color=S.MUTED, marker="^", lw=0, ms=4,
                        label=f"97.5% limit above {ymax:.0f} (value given)"))
    fig.legend(handles=h, loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=(0, 0.07, 1, 0.98))
    S.save(fig, "Figure_4", FIG); plt.close(fig)


# ---- Figures 5 and S2: empirical CV and relative RMSE ------------------------------------------
def cv_like(metric, ylabel, name, ref_line=None, poisson=False, ycap=3.0):
    """Flagged gimbal scenarios: open grey-edged symbols, not joined by the line; values above ycap are
    drawn at ycap with an arrow and the value. ycap is above every non-flagged value."""
    ok = m[~flagged(m) & m["design"].isin(S.SIMULATED)]
    if ok[metric].max() >= ycap:
        raise SystemExit(f"{name}: a non-flagged value exceeds the cap {ycap}")
    ymin = 0.7 * m.loc[m["design"].isin(S.SIMULATED), metric].min()
    fig, axes = panel_grid()
    for ax, b in zip(axes, S.BLOCKS):
        if ref_line is not None:
            ax.axhline(ref_line, color=S.INK, lw=0.8, ls="--", zorder=1)
        if poisson:   # approximations computed in R (compute_metrics.R), 0% overlap design
            s0 = m[(m["design"] == "0") & (m["blocks"] == b)].sort_values("true_density")
            ax.plot(s0["true_density"], s0["cv_poisson"], color="#9A9A9A", lw=0.8, zorder=1)
            ax.plot(s0["true_density"], s0["cv_binomial"], color="black", lw=0.8, zorder=2.5)  # above the 0% line
        for d in S.SIMULATED:
            s = m[(m["design"] == d) & (m["blocks"] == b)].sort_values("true_density")
            fl = flagged(s)
            f = s[fl]
            below = f[f[metric] <= ycap]
            ax.plot(below["true_density"], below[metric], ls="none", marker=S.DESIGN_MARKERS[d],
                    ms=4.5, mfc="white", mec=FLAG_EDGE, mew=0.9, zorder=3)
            above = f[f[metric] > ycap]
            if len(above):   # an empty unclipped artist would sit at the figure origin and break the layout
                ax.plot(above["true_density"], [ycap] * len(above), ls="none", marker="^", ms=4, color=FLAG_EDGE,
                        clip_on=False, zorder=4)
            for i, (xv, v) in enumerate(zip(above["true_density"], above[metric])):
                # stagger the labels of neighbouring capped points so they do not overlap
                t = ax.annotate(cap_label(v), (xv, ycap), xytext=(0, 4 + 6 * (i % 2)), textcoords="offset points",
                                ha="center", va="bottom", fontsize=5.5, color=S.MUTED, annotation_clip=False)
                t.set_in_layout(False)
            s = s[~fl]
            ax.plot(s["true_density"], s[metric], color=S.DESIGN_COLORS[d], lw=1.5, zorder=2)
            for fail, sub in s.groupby("failed_gt5pct"):
                ax.plot(sub["true_density"], sub[metric], ls="none", marker=S.DESIGN_MARKERS[d], ms=4.5,
                        color=S.DESIGN_COLORS[d], mfc="white" if fail else S.DESIGN_COLORS[d],
                        mec=S.DESIGN_COLORS[d] if fail else "white", mew=0.8 if fail else 0.5, zorder=3)
        ax.set_yscale("log"); ax.set_ylim(ymin, ycap); ax.set_title(f"{b} blocks", loc="left", pad=22)
        ax.set_xticks([50, 200, 400, 600, 800])
    for ax in axes[2:]:
        ax.set_xlabel(S.DENS_LABEL)
    for ax in axes[::2]:
        ax.set_ylabel(ylabel)
    # the open-symbol entry only when some scenario has >5% failed replicates
    h = design_handles(S.SIMULATED, open_note=bool(m["failed_gt5pct"].any()), flag_note=ANY_FLAG)
    if ANY_FLAG and (m.loc[flagged(m), metric] > ycap).any():
        h.append(Line2D([], [], color=FLAG_EDGE, marker="^", lw=0, ms=4, label=f"Value above {ycap:g} (value given)"))
    if ref_line is not None:
        h.append(Line2D([], [], color=S.INK, lw=0.8, ls="--", label=f"{ylabel.split(' (')[0]} = {ref_line}"))
    if poisson:
        h.append(Line2D([], [], color="#9A9A9A", lw=0.8, label="Poisson approximation (detections only)"))
        h.append(Line2D([], [], color="black", lw=0.8, label="Binomial approximation (detections and proportion covered)"))
    fig.legend(handles=h, loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.04))
    fig.tight_layout(rect=(0, 0.08, 1, 0.98))
    S.save(fig, name, FIG); plt.close(fig)


# ---- Figure 6: density threshold against effort (total transect length) -------------------------
def figure_6():
    fig, ax = plt.subplots(figsize=(5.8, 3.9))
    NOT_MET = {"65": 1130, "10": 1095, "0": 1060, "200": 1025, "gimbal": 990}
    LAB_OFF = {"65": (0, 5), "10": (6, 3), "0": (-7, -9), "200": (0, -10), "gimbal": (0, 5)}
    ax.axhspan(960, 1160, color="#F2F2F2", zorder=0)
    for d in S.DESIGN_ORDER:
        t = thr[thr["design"] == d].sort_values("blocks")
        met = t[t["threshold_density"].notna()]
        ax.plot(met["transect_length_km"], met["threshold_density"], color=S.DESIGN_COLORS[d],
                marker=S.DESIGN_MARKERS[d], ms=5.5, mec="white", mew=0.5, zorder=3)
        nm = t[t["threshold_density"].isna()]
        ax.plot(nm["transect_length_km"], [NOT_MET[d]] * len(nm), ls="none", marker=S.DESIGN_MARKERS[d], ms=5,
                mfc="white", mec=S.DESIGN_COLORS[d], mew=1, zorder=3)
        for _, r in t.iterrows():
            met_r = pd.notna(r["threshold_density"])
            y = r["threshold_density"] if met_r else NOT_MET[d]
            off = LAB_OFF[d] if met_r else (7, -2.5)
            ax.annotate(f"{int(r['blocks'])}", (r["transect_length_km"], y), xytext=off, textcoords="offset points",
                        ha="center", fontsize=5.5, color=S.MUTED)
    ax.set_xscale("log"); ax.set_xlim(15, 7000); ax.set_ylim(0, 1160)
    ax.set_yticks([0, 100, 300, 500, 700, 900, 1060]); ax.set_yticklabels(["0", "100", "300", "500", "700", "900", "not met"])
    ax.set_xlabel("Total transect length (km), log scale")
    ax.set_ylabel("Lowest density meeting CV ≤ 0.15 and |bias| ≤ 10%\n(individuals per 1,000 km$^2$)")
    h = design_handles(S.DESIGN_ORDER) + [Line2D([], [], color=S.MUTED, marker="o", mfc="white", lw=0, ms=5,
                                                 label="Not met at any simulated density (≤ 900)")]
    ax.legend(handles=h, loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize=6.8)
    ax.text(0, -0.25, "Numbers next to symbols: number of blocks. 65% and 10% overlap use the 0% overlap estimates.\n"
            "Gimbal scenarios with < 20 detections (not reliably estimable) never count as meeting the criteria.",
            transform=ax.transAxes, ha="left", fontsize=6, color=S.MUTED)
    fig.tight_layout()
    S.save(fig, "Figure_6", FIG); plt.close(fig)


if __name__ == "__main__":
    figure_1(); figure_2(); figure_3(); figure_4()
    cv_like("emp_cv", "Empirical CV (SD / mean)", "Figure_5", ref_line=0.15, poisson=True)
    cv_like("rrmse", "Relative RMSE (RMSE / true)", "Figure_S2")
    figure_6()
    print("Figures written to", FIG)
