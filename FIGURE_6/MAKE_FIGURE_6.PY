"""Recreate manuscript Figure 6 from the repository plot-data pickle."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, LogFormatterMathtext
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "PLOT_DATA.PKL"
CHARGE = "Charge (h+/surface atom)"
RAW = "raw values"
MODELS = [
    ("Neutral linear 2D", "neutral", "linear", 2),
    ("Full linear 2D", "full", "linear", 2),
    ("Neutral nonlinear 2D", "neutral", "nonlinear", 2),
    ("Full nonlinear 2D", "full", "nonlinear", 2),
    ("Full linear 4D", "full", "linear", 4),
    ("Original Julia 4D", "full", "nonlinear", 4),
]
COLORS = {2: "#0029FF", 4: "#FF0500"}
LINE_STYLES = {"linear": "-", "nonlinear": (0, (5, 2.5))}
MARKERS = {"neutral": "o", "full": "s"}
RAW_COLOR = "#252525"

df = pd.read_pickle(SOURCE).sort_values(CHARGE)
assert len(df) == 23
assert np.allclose(df[CHARGE].to_numpy(), np.round(np.arange(-0.11, 0.111, 0.01), 2))
values = df[[name for name, _, _, _ in MODELS]].to_numpy(dtype=float)
assert np.isfinite(values).all() and (values > 0).all()
raw = df.loc[df[RAW].notna()]
assert len(raw) == 5
neutral = df.loc[np.isclose(df[CHARGE], 0)].iloc[0]

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "mathtext.fontset": "dejavusans",
        "font.size": 20.5,
        "axes.labelsize": 21.5,
        "axes.linewidth": 0.9,
        "xtick.labelsize": 18.8,
        "ytick.labelsize": 18.8,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.major.size": 4,
        "ytick.major.size": 4,
        "savefig.dpi": 1200,
    }
)

fig = plt.figure(figsize=(10.6, 10.3), facecolor="white")
axes = [fig.add_axes([0.135, 0.595, 0.37, 0.235]),
        fig.add_axes([0.595, 0.595, 0.37, 0.235])]


def style_box(ax):
    ax.grid(False)
    ax.tick_params(axis="both", pad=4)
    for spine in ax.spines.values():
        spine.set_visible(True)


for ax, subset, tag in [
    (axes[0], df.loc[df[CHARGE] <= 0].iloc[::-1], "(a)"),
    (axes[1], df.loc[df[CHARGE] >= 0.01], "(b)"),
]:
    for name, training, form, n_descriptors in MODELS:
        color = COLORS[n_descriptors]
        marker = MARKERS[training]
        ax.plot(subset[CHARGE], subset[name], lw=2.05, color=color,
                ls=LINE_STYLES[form], marker=marker, markersize=4.0,
                markerfacecolor="white", markeredgecolor=color,
                markeredgewidth=0.9, solid_capstyle="round", zorder=2)
    raw_subset = raw.loc[raw[CHARGE].isin(subset[CHARGE])]
    ax.scatter(raw_subset[CHARGE], raw_subset[RAW], marker="o", s=53,
               facecolor=RAW_COLOR, edgecolor="white", lw=0.75, zorder=5)
    ax.set_yscale("log")
    ax.yaxis.set_major_formatter(LogFormatterMathtext())
    tag_x, tag_align = (0.98, "right") if tag == "(a)" else (0.02, "left")
    ax.text(tag_x, 0.04, tag, transform=ax.transAxes, va="bottom", ha=tag_align,
            fontsize=22, fontweight="bold")
    style_box(ax)

axes[0].set_xlim(-0.115, 0.004)
axes[0].set_xticks([-0.11, -0.05, 0])
axes[0].set_xticklabels([r"$-\frac{2}{18}$", r"$-\frac{1}{18}$", "0"])
axes[0].set_ylim(1e-8, 5e2)
axes[0].yaxis.set_major_locator(FixedLocator([1e-7, 1e-4, 1e-1, 1e2]))
axes[0].set_ylabel(r"H$_2$ TOF (min$^{-1}$)", labelpad=6)
axes[0].yaxis.set_label_coords(-0.233, 0.5)

axes[1].set_xlim(0.006, 0.114)
axes[1].set_xticks([0.05, 0.11])
axes[1].set_xticklabels([r"$+\frac{1}{18}$", r"$+\frac{2}{18}$"])
axes[1].set_ylim(1e-36, 5e2)
axes[1].yaxis.set_major_locator(FixedLocator(10.0 ** np.array([-35., -28., -21., -14., -7., 2.])))
fig.text(0.55, 0.545, "Charge",
         ha="center", va="center", fontsize=19)

heat_ax = fig.add_axes([0.135, 0.30, 0.83, 0.18])
family_names = [
    "Neutral\nLinear\n$\Delta E^{\mathrm{ad}}_{\mathrm{CH}^*}$,\n$\Delta E^{\mathrm{ad}}_{\mathrm{HCO}^*}$",
    "Full\nLinear\n$\Delta E^{\mathrm{ad}}_{\mathrm{CO}_2^*}$,\n$\Delta E^{\mathrm{ad}}_{\mathrm{COOH}^*}$",
    "Neutral\nNonlinear\n$\Delta E^{\mathrm{ad}}_{\mathrm{OH}^*}$,\n$\Delta E^{\mathrm{ad}}_{\mathrm{COOH}^*}$",
    "Full\nNonlinear\n$\Delta E^{\mathrm{ad}}_{\mathrm{CO}^*_{\mathrm{FCC}}}$,\n$\Delta E^{\mathrm{ad}}_{\mathrm{H}_2\mathrm{O}^*}$",
    "Full\nLinear\n$\Delta E^{\mathrm{ad}}_{\mathrm{CH}^*}$,\n$\Delta E^{\mathrm{ad}}_{\mathrm{HCO}^*}$,\n$\Delta E^{\mathrm{ad}}_{\mathrm{H}_2\mathrm{O}^*}$,\n$\Delta E^{\mathrm{ad}}_{\mathrm{O}^*}$",
    "Full\nNonlinear\n$\Delta E^{\mathrm{ad}}_{\mathrm{CH}^*}$,\n$\Delta E^{\mathrm{ad}}_{\mathrm{HCO}^*}$,\n$\Delta E^{\mathrm{ad}}_{\mathrm{H}_2\mathrm{O}^*}$,\n$\Delta E^{\mathrm{ad}}_{\mathrm{O}^*}$",
]
ratio = np.log10(raw[[name for name, *_ in MODELS]].to_numpy(dtype=float)
                 / raw[RAW].to_numpy(dtype=float)[:, None])
cmap=plt.get_cmap("RdBu_r").copy()
cmap.set_bad("#dddddd")
heat = heat_ax.imshow(np.ma.masked_invalid(ratio), cmap=cmap,
                       norm=TwoSlopeNorm(vmin=-3, vcenter=0, vmax=3),
                       aspect="auto", interpolation="nearest")
heat_ax.set_xticks(range(len(family_names)), family_names,
                   fontsize=14.0, linespacing=1.25)
charge_values = raw[CHARGE].to_numpy(dtype=float)
heat_ax.set_yticks(range(len(raw)),
                   [r"$-\frac{2}{18}$", r"$-\frac{1}{18}$", "0",
                    r"$+\frac{1}{18}$", r"$+\frac{2}{18}$"])
heat_ax.set_ylabel("Charge", labelpad=8)
heat_ax.yaxis.set_label_coords(-0.104, 0.5)
heat_ax.set_xticks(np.arange(-0.5, len(family_names), 1), minor=True)
heat_ax.set_yticks(np.arange(-0.5, len(raw), 1), minor=True)
heat_ax.grid(which="minor", color="white", linewidth=1.4)
heat_ax.tick_params(axis="both", length=0, pad=5)
heat_ax.tick_params(which="minor", bottom=False, left=False)
for i in range(ratio.shape[0]):
    for j in range(ratio.shape[1]):
        value = ratio[i, j]
        shown = 0.0 if abs(value) < 0.05 else value
        label = f"{shown:+.1f}"
        heat_ax.text(j, i, label, ha="center", va="center",
                     fontsize=18.0, fontweight="semibold",
                     color="white" if abs(value) > 1.7 else "#202020")
style_box(heat_ax)
heat_ax.text(0.0, 1.03, "(c)", transform=heat_ax.transAxes,
             va="bottom", ha="left", fontsize=22, fontweight="bold")
color_ax = fig.add_axes([0.28, 0.070, 0.54, 0.022])
colorbar = fig.colorbar(heat, cax=color_ax, orientation="horizontal", extend="both")
colorbar.set_ticks([-3, -2, -1, 0, 1, 2, 3])
colorbar.ax.tick_params(labelsize=16, length=3, pad=2)
colorbar.ax.xaxis.set_label_position("top")
colorbar.set_label(r"$\log_{10}(\mathrm{TOF}_{\mathrm{SR}}/\mathrm{TOF}_{\mathrm{DFT}})$",
                   fontsize=18, labelpad=8)

legend = [
    Line2D([0], [0], color=COLORS[2], lw=2.2, label="2 descriptors"),
    Line2D([0], [0], color=COLORS[4], lw=2.2, label="4 descriptors"),
    Line2D([0], [0], color=RAW_COLOR, lw=2.2, ls=LINE_STYLES["linear"], label="Linear"),
    Line2D([0], [0], color=RAW_COLOR, lw=2.2, ls=LINE_STYLES["nonlinear"], label="Nonlinear"),
    Line2D([0], [0], marker=MARKERS["neutral"], ls="none", markersize=6,
           markerfacecolor="white", markeredgecolor=RAW_COLOR, label="Neutral-only fit"),
    Line2D([0], [0], marker=MARKERS["full"], ls="none", markersize=6,
           markerfacecolor="white", markeredgecolor=RAW_COLOR, label="Full-dataset fit"),
    Line2D([0], [0], marker="o", ls="none", markersize=6,
           markerfacecolor=RAW_COLOR, markeredgecolor=RAW_COLOR, label="DFT energetics"),
]
# Matplotlib fills legend columns first, placing each paired attribute vertically.
legend_grid = [legend[0], legend[1],
               legend[2], legend[3],
               legend[4], legend[5]]
fig.legend(handles=legend_grid, loc="upper center", bbox_to_anchor=(0.55, 0.987),
           ncol=3, frameon=False, fontsize=18.5, handlelength=1.8,
           handletextpad=0.35, columnspacing=0.8, labelspacing=0.3)
fig.legend(handles=[legend[6]], loc="upper center", bbox_to_anchor=(0.55, 0.925),
           ncol=1, frameon=False, fontsize=18.5, handlelength=1.8,
           handletextpad=0.35, columnspacing=0.8)

output = HERE / "FIGURE_6.PNG"
fig.savefig(output, dpi=1200, facecolor="white", bbox_inches="tight", pad_inches=0.08)
plt.close(fig)
print(output)
