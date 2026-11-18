"""Create Supporting Information Figure S15 from the energy pickle.

Run from anywhere with: python FIGURE_S15/MAKE_FIGURE_S15.PY
The output is FIGURE_S15.PNG at 1200 DPI.
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
TS_MARKERS = {"HO--H": "o", "OC--OH": "s", "O--CO": "D", "H--CO2": "v"}
TS_LABELS = {
    "HO--H": r"HO$-$H",
    "OC--OH": r"OC$-$OH",
    "O--CO": r"O$-$CO",
    "H--CO2": r"H$-$CO$_2$",
}


def plot_series(ax, rows, species, marker):
    rows = rows.sort_values("Charge")
    x = rows["Charge"].to_numpy(dtype=float)
    y = rows[species].to_numpy(dtype=float)
    finite = np.isfinite(y)
    # Mask missing saddle points so the line does not span them.
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
    ax.set_box_aspect(1.0)
    ax.set_ylim(0.10, 4.12)
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))


def add_legends(fig):
    marker_handles = [
        Line2D(
            [], [], ls="none", marker=shape, markersize=8.5,
            markerfacecolor="white", markeredgecolor="#343434",
            markeredgewidth=1.0, label=TS_LABELS[name],
        )
        for name, shape in TS_MARKERS.items()
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
    fig, axes = plt.subplots(1, 3, figsize=(12.8, 6.2))
    fig.subplots_adjust(
        left=0.08, right=0.985, top=0.91, bottom=0.31, wspace=0.24,
    )
    for index, (ax, metal) in enumerate(zip(axes, ("Cu", "Ni", "Au"))):
        rows = data.loc[data.Metal == metal]
        for species, marker in TS_MARKERS.items():
            plot_series(ax, rows, species, marker)
        style_axes(ax)
        ax.set_title(f"{metal}(111)", fontsize=17, fontweight="bold", pad=5)
        ax.text(
            0.0, 1.025, f"({chr(97 + index)})", transform=ax.transAxes,
            va="bottom", fontsize=15.5, fontweight="bold",
        )
    fig.supylabel("Activation energy (eV)", x=0.012, fontsize=17)
    fig.supxlabel(
        r"Surface charge ($h^+$/surface metal atom)", y=0.205, fontsize=17,
    )
    add_legends(fig)
    fig.savefig(HERE / "FIGURE_S15.PNG", dpi=1200, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
