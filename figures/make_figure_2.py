"""Generate manuscript Figure 2 from the original full-precision pickle.

Run ``python figures/make_figure_2.py`` from the repository root. The PNG and
PDF are written beside this script. This script does not edit the manuscript.

The PNG was verified byte-for-byte against manuscript Figure 2 with Python
3.10.13, NumPy 2.2.6, pandas 2.3.1, and Matplotlib 3.10.3. Its SHA-256 is
0a5839fb9fdd957caeb3441be3b5f7d49cfdcbf7dda8dd85a93041cbbb9161d9.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "adsorption_TS_energies.pkl"
CHARGES = (-0.11, -0.05, 0.0, 0.05, 0.11)
COLORS = ("#748eba", "#b4c7e7", "#FFFFFF", "#f8cbad", "#d69163")
COLOR = dict(zip(CHARGES, COLORS))
METAL_MARKERS = {"Pt": "o", "Cu": "s", "Ru": "^", "Ni": "D", "Au": "v"}
TS_MARKERS = {"HO--H": "o", "OC--OH": "s", "O--CO": "D", "H--CO2": "v"}
TS_LABELS = {
    "HO--H": r"HO$-$H", "OC--OH": r"OC$-$OH",
    "O--CO": r"O$-$CO", "H--CO2": r"H$-$CO$_2$",
}


def draw_series(ax, rows, column, marker):
    rows = rows.sort_values("Charge")
    x = rows["Charge"].to_numpy(float)
    y = rows[column].to_numpy(float)
    finite = np.isfinite(y)
    # A missing saddle point should not be bridged by a line.
    ax.plot(x, np.ma.masked_where(~finite, y), color="#646464", lw=1.25, zorder=1)
    for charge, energy in zip(x[finite], y[finite]):
        ax.scatter(
            charge, energy, marker=marker, s=78,
            facecolor=COLOR[float(charge)], edgecolor="#343434",
            linewidth=1.15, zorder=3,
        )


def marker_legend(ax, labels, markers, title, location, columns=1):
    handles = [
        Line2D([], [], ls="none", marker=markers[name], markersize=8.4,
               markerfacecolor="white", markeredgecolor="#343434",
               markeredgewidth=1.0, label=labels[name])
        for name in markers
    ]
    ax.legend(
        handles=handles, title=title, ncol=columns, frameon=True,
        loc=location, borderaxespad=0.5,
        handlelength=0.85 if columns == 2 else 2.0,
        handletextpad=0.3 if columns == 2 else 0.5,
        columnspacing=0.35 if columns == 2 else 2.0,
        labelspacing=0.24, borderpad=0.2 if columns == 2 else 0.32,
        fontsize=14.5, title_fontsize=14.5,
        facecolor="white", edgecolor="none", framealpha=0.9,
    )


def main():
    data = pd.read_pickle(DATA)
    assert len(data) == 25
    assert set(data["Charge"].unique()) == set(CHARGES)
    pt = data.loc[data.Metal == "Pt"]
    assert len(pt) == 5

    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 15,
        "mathtext.fontset": "dejavusans", "axes.linewidth": 0.95,
        "xtick.major.size": 4, "ytick.major.size": 4,
        "xtick.direction": "in", "ytick.direction": "in",
    })
    fig, axes = plt.subplots(1, 2, figsize=(8.0, 4.0), sharex=True)
    # Match the square panel boxes used in the manuscript's parity figure.
    for ax, left in zip(axes, (0.11, 0.61)):
        ax.set_position((left, 0.25, 0.35, 0.70))

    for metal, marker in METAL_MARKERS.items():
        draw_series(axes[0], data.loc[data.Metal == metal], "OH*", marker)
    for ts, marker in TS_MARKERS.items():
        draw_series(axes[1], pt, ts, marker)

    for ax in axes:
        ax.set_xlim(-0.125, 0.125)
        ax.set_xticks(CHARGES)
        ax.set_xticklabels((r"$-\frac{2}{18}$", r"$-\frac{1}{18}$", "0",
                            r"$+\frac{1}{18}$", r"$+\frac{2}{18}$"))
        ax.tick_params(labelsize=15, pad=4)
        ax.grid(False)
        for spine in ax.spines.values():
            spine.set_visible(True)
        ax.set_box_aspect(1.0)

    axes[0].set_ylim(-4.72, -1.78)
    axes[1].set_ylim(0.48, 2.04)
    axes[1].set_yticks((0.5, 1.0, 1.5, 2.0))
    axes[1].set_yticklabels(("0.5", "1.0", "1.5", "2.0"))
    axes[0].set_ylabel(r"$\Delta E^{\mathrm{ad}}_{\mathrm{OH}^{*}}$ (eV)", fontsize=17)
    axes[1].set_ylabel(r"$\Delta E^{\mathrm{act}}$ (eV)", fontsize=17)
    for index, ax in enumerate(axes):
        ax.text(0.975 if index == 0 else 0.80, 0.965,
                f"({chr(ord('a') + index)})",
                transform=ax.transAxes, va="top", ha="right",
                fontsize=18, fontweight="bold",
                bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.4})

    marker_legend(axes[0], {name: name for name in METAL_MARKERS},
                  METAL_MARKERS, "Metals", "lower left", columns=2)
    marker_legend(axes[1], TS_LABELS, TS_MARKERS, "Pt TSs", "upper left")
    fig.supxlabel(r"Surface charge ($h^+$/surface metal atom)",
                 y=0.095, fontsize=17)
    fig.text(0.035, 0.025, "Charge", ha="left", va="center",
             fontsize=14.5)
    for x, charge in zip((0.25, 0.39, 0.53, 0.67, 0.81), CHARGES):
        fig.add_artist(Line2D([x], [0.025], ls="none", marker="o",
                              markersize=8.4, markerfacecolor=COLOR[charge],
                              markeredgecolor="#343434", markeredgewidth=1.0,
                              transform=fig.transFigure))
        label = {-.11: r"$-\frac{2}{18}$", -.05: r"$-\frac{1}{18}$",
                 0.0: "0", .05: r"$+\frac{1}{18}$",
                 .11: r"$+\frac{2}{18}$"}[charge]
        fig.text(x + 0.025, 0.025, label,
                 ha="left", va="center", fontsize=14.5)

    fig.savefig(HERE / "fig_02_charge_energies.pdf", bbox_inches="tight")
    fig.savefig(HERE / "fig_02_charge_energies.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
