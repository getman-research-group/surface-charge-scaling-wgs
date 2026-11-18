"""Create Supporting Information Figure S14 from the energy pickle.

Run from anywhere with: python FIGURE_S14/MAKE_FIGURE_S14.PY
The output is FIGURE_S14.PNG at 1200 DPI.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FormatStrFormatter
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "ADSORPTION_TS_ENERGIES.PKL"
CHARGES = (-0.11, -0.05, 0.0, 0.05, 0.11)
COLORS = ("#748eba", "#b4c7e7", "#FFFFFF", "#f8cbad", "#d69163")
CHARGE_COLOR = dict(zip(CHARGES, COLORS))
CHARGE_LABELS = (
    r"$-\frac{2}{18}$", r"$-\frac{1}{18}$", r"$0$",
    r"$+\frac{1}{18}$", r"$+\frac{2}{18}$",
)
METAL_MARKERS = {"Pt": "o", "Cu": "s", "Ru": "^", "Ni": "D", "Au": "v"}
ADSORBATES = (
    ("H*", r"$\mathrm{H}^{*}$"),
    ("O*", r"$\mathrm{O}^{*}$"),
    ("CH*", r"$\mathrm{CH}^{*}$"),
    ("CO* FCC", r"$\mathrm{CO}^{*}_{\mathrm{FCC}}$"),
    ("CO2*", r"$\mathrm{CO}_{2}^{*}$"),
    ("H2O*", r"$\mathrm{H}_{2}\mathrm{O}^{*}$"),
    ("HCO*", r"$\mathrm{HCO}^{*}$"),
    ("COOH*", r"$\mathrm{COOH}^{*}$"),
)


def plot_series(ax, rows, species, marker):
    rows = rows.sort_values("Charge")
    x = rows["Charge"].to_numpy(dtype=float)
    y = rows[species].to_numpy(dtype=float)
    finite = np.isfinite(y)
    ax.plot(x, np.ma.masked_where(~finite, y), color="#707070", lw=1.15, zorder=1)
    for charge, energy in zip(x[finite], y[finite]):
        ax.scatter(
            charge, energy, marker=marker, s=69,
            facecolor=CHARGE_COLOR[round(float(charge), 2)],
            edgecolor="#343434", linewidth=1.0, zorder=3,
        )


def style_axes(ax):
    ax.set_xlim(-0.125, 0.125)
    ax.set_xticks(CHARGES)
    ax.set_xticklabels(CHARGE_LABELS)
    ax.tick_params(axis="both", direction="in", labelsize=12.5, length=4)
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_visible(True)
    ax.set_box_aspect(0.86)


def add_legends(fig):
    marker_handles = [
        Line2D(
            [], [], ls="none", marker=shape, markersize=8.5,
            markerfacecolor="white", markeredgecolor="#343434",
            markeredgewidth=1.0, label=name,
        )
        for name, shape in METAL_MARKERS.items()
    ]
    fig.legend(
        handles=marker_handles, loc="lower center", ncol=len(marker_handles),
        bbox_to_anchor=(0.5, 0.072), frameon=False, fontsize=13.5,
        handletextpad=0.35, columnspacing=1.35,
    )
    color_handles = [
        Line2D(
            [], [], ls="none", marker="o", markersize=8.5,
            markerfacecolor=CHARGE_COLOR[charge], markeredgecolor="#343434",
            markeredgewidth=1.0, label=label,
        )
        for charge, label in zip(CHARGES, CHARGE_LABELS)
    ]
    fig.legend(
        handles=color_handles, loc="lower center", ncol=5,
        bbox_to_anchor=(0.5, 0.019), frameon=False, fontsize=13.5,
        handletextpad=0.35, columnspacing=1.25,
    )
    fig.text(0.105, 0.052, "Charge", ha="left", va="center", fontsize=13.5)


def main():
    data = pd.read_pickle(DATA)
    if len(data) != 25 or set(np.round(data.Charge.unique(), 2)) != set(CHARGES):
        raise ValueError("Expected five metals at five charge states")
    if not data.groupby("Metal").size().eq(5).all():
        raise ValueError("Each metal must have five charge states")
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 15,
        "mathtext.fontset": "dejavusans", "axes.linewidth": 0.95,
    })
    fig, axes = plt.subplots(2, 4, figsize=(13.2, 8.2))
    fig.subplots_adjust(
        left=0.075, right=0.985, top=0.91, bottom=0.235,
        wspace=0.26, hspace=0.28,
    )
    for index, (ax, (species, label)) in enumerate(zip(axes.flat, ADSORBATES)):
        for metal, marker in METAL_MARKERS.items():
            plot_series(ax, data.loc[data.Metal == metal], species, marker)
        style_axes(ax)
        ax.set_title(label, fontsize=17, fontweight="bold", pad=5)
        ax.text(
            0.0, 1.025, f"({chr(97 + index)})", transform=ax.transAxes,
            va="bottom", fontsize=15.5, fontweight="bold",
        )
        if index < 4:
            ax.tick_params(labelbottom=False)
        minimum, maximum = data[species].min(), data[species].max()
        pad = max(0.07, (maximum - minimum) * 0.11)
        ax.set_ylim(minimum - pad, maximum + pad)
        ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    fig.supylabel(
        r"$\Delta E^{\mathrm{ad}}_{s*}$ (eV)",
        x=0.011, y=(0.91 + 0.235) / 2, fontsize=17,
    )
    fig.supxlabel(
        r"Surface charge ($h^+$/surface metal atom)", y=0.18, fontsize=17,
    )
    add_legends(fig)
    fig.savefig(HERE / "FIGURE_S14.PNG", dpi=1200, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
