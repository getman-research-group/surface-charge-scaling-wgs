"""Make SI Figure S16 from the saved interpolation analysis.

Run the analysis first:
    python PT_DESCRIPTOR_INTERPOLATION/ANALYZE_PT_DESCRIPTOR_INTERPOLATION.PY
Then run this script:
    python FIGURE_S16/MAKE_FIGURE_S16.PY
The output is FIGURE_S16.PNG at 1200 DPI.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter, MaxNLocator
import pandas as pd


HERE = Path(__file__).resolve().parent
ANALYSIS = HERE.parent / "PT_DESCRIPTOR_INTERPOLATION"
CURVES = ANALYSIS / "INTERPOLATION_CURVES.PKL"
NODES = ANALYSIS / "INTERPOLATION_NODES.PKL"
SERIES = (
    ("CH*", r"CH$^*$"),
    ("HCO*", r"HCO$^*$"),
    ("O*", r"O$^*$"),
    ("H2O*", r"H$_2$O$^*$"),
    ("CO2*", r"CO$_2^*$"),
    ("COOH*", r"COOH$^*$"),
    ("H*", r"H$^*$"),
    ("CO* FCC", r"CO$^*_{\mathrm{FCC}}$"),
)
CHARGE_TICKS = (-0.11, -0.05, 0.0, 0.05, 0.11)
CHARGE_LABELS = (
    r"$-\frac{2}{18}$", r"$-\frac{1}{18}$", r"$0$",
    r"$+\frac{1}{18}$", r"$+\frac{2}{18}$",
)


def main():
    curves = pd.read_pickle(CURVES)
    nodes = pd.read_pickle(NODES)
    if len(curves) != 4509 or len(nodes) != 45:
        raise ValueError("Run the S16 interpolation analysis before plotting")
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 11,
        "axes.titlesize": 12,
        "axes.labelsize": 11,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 10.5,
    })
    fig, axes = plt.subplots(2, 4, figsize=(13.2, 6.8), sharex=True)
    fig.subplots_adjust(
        left=0.075, right=0.985, top=0.875, bottom=0.14,
        wspace=0.35, hspace=0.28,
    )
    for ax, (descriptor, title) in zip(axes.flat, SERIES):
        curve = curves.loc[curves.DESCRIPTOR.eq(descriptor)].sort_values("CHARGE")
        dft = nodes.loc[nodes.DESCRIPTOR.eq(descriptor)].sort_values("CHARGE")
        if len(curve) != 501 or len(dft) != 5:
            raise ValueError(f"Missing interpolation data for {descriptor}")
        ax.plot(
            curve.CHARGE, curve.LAGRANGE_EV,
            color="#5A3D91", lw=2.2,
            label="Fourth-order Lagrange interpolation",
        )
        ax.plot(
            curve.CHARGE, curve.QUADRATIC_EV,
            color="#D4874E", lw=1.9, ls="--", label="Quadratic fit",
        )
        ax.scatter(
            dft.CHARGE, dft.DFT_EV,
            s=40, facecolor="white", edgecolor="black",
            lw=1.15, zorder=3, label="DFT-calculated energy on Pt(111)",
        )
        ax.set_title(title, fontweight="bold", pad=6)
        ax.set_xlim(-0.115, 0.115)
        ax.set_xticks(CHARGE_TICKS, labels=CHARGE_LABELS)
        ax.yaxis.set_major_locator(MaxNLocator(nbins=5, steps=[1, 2, 5, 10]))
        ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
        ax.tick_params(direction="in", top=True, right=True)
        for spine in ax.spines.values():
            spine.set_visible(True)
    fig.supxlabel("Charge", y=0.045, fontsize=12)
    fig.supylabel(r"$\Delta E^{\mathrm{ad}}_{s*}$ (eV)", x=0.015, fontsize=12)
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(
        handles, labels, loc="upper center",
        bbox_to_anchor=(0.5, 0.99), ncol=3, frameon=False,
    )
    fig.savefig(HERE / "FIGURE_S16.PNG", dpi=1200, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
