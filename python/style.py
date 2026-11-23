"""Shared figure style: design order, colours (Okabe-Ito, colour-blind safe), markers, labels."""
import matplotlib as mpl
DESIGN_ORDER = ["65", "10", "0", "200", "gimbal"]
DESIGN_LABELS = {"65": "65% overlap", "10": "10% overlap", "0": "0% overlap",
                 "200": "200 m spacing", "gimbal": "Gimbal"}
DESIGN_SHORT = {"65": "65%", "10": "10%", "0": "0%", "200": "200 m", "gimbal": "Gimbal"}
DESIGN_COLORS = {"65": "#E69F00", "10": "#56B4E9", "0": "#009E73", "200": "#0072B2", "gimbal": "#D55E00"}
DESIGN_MARKERS = {"65": "^", "10": "v", "0": "o", "200": "s", "gimbal": "D"}
SIMULATED = ["0", "200", "gimbal"]
BLOCKS = [4, 9, 16, 25]
DENS_LABEL = "True density (individuals per 1,000 km$^2$)"
INK, MUTED, GRID = "#222222", "#666666", "#DDDDDD"
def apply():
    mpl.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 8, "axes.titlesize": 8.5, "axes.labelsize": 8,
        "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7, "axes.edgecolor": MUTED,
        "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK,
        "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
        "grid.linewidth": 0.5, "axes.axisbelow": True, "lines.linewidth": 1.5, "savefig.dpi": 300,
        "legend.frameon": False})
# Export: PNG at 300 dpi. For journal submission set DPI = 600 and TIFF = True.
DPI = 300
TIFF = False
def save(fig, name, outdir):
    import os
    os.makedirs(outdir, exist_ok=True)
    fig.savefig(os.path.join(outdir, f"{name}.png"), dpi=DPI, bbox_inches="tight")
    if TIFF:
        fig.savefig(os.path.join(outdir, f"{name}.tiff"), dpi=DPI, bbox_inches="tight",
                    pil_kwargs={"compression": "tiff_lzw"})
