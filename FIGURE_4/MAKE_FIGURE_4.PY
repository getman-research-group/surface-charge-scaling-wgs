"""Make manuscript Figure 4 from the archived two-descriptor result PKLs.

Run from anywhere with: python FIGURE_4/MAKE_FIGURE_4.PY
The output is FIGURE_4.PNG at 1200 DPI.
"""

from itertools import combinations
from pathlib import Path
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, NullFormatter
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
ENERGIES = ROOT / "ADSORPTION_TS_ENERGIES.PKL"
ADSORBATES = (
    "H*", "O*", "OH*", "CH*", "CO* FCC", "CO2*", "H2O*", "HCO*", "COOH*"
)
TS = ("HO--H", "O--CO", "H--CO2", "OC--OH")
MODELS = {
    "Neutral only": ("LINEAR_NEUTRAL", "NONLINEAR_NEUTRAL"),
    "All charges": ("LINEAR_ALL_CHARGES", "NONLINEAR_ALL_CHARGES"),
}


def result_name(pair):
    def slug(species):
        return re.sub(r"[^A-Z0-9]+", "_", species.upper()).strip("_")

    return f"RESULTS_{slug(pair[0])}_{slug(pair[1])}.PKL"


def pooled_rmse(frame, counts, column):
    weights = frame["TARGET"].map(counts)
    if weights.isna().any():
        raise ValueError("A target is absent from the energy data")
    errors = frame[column].to_numpy(dtype=float)
    return float(np.sqrt(np.average(errors ** 2, weights=weights.to_numpy(dtype=float))))


def load_scores():
    energies = pd.read_pickle(ENERGIES)
    pairs = list(combinations(ADSORBATES, 2))
    rows = []
    for dataset, folders in MODELS.items():
        observations = energies.loc[energies["Charge"].eq(0)] if dataset == "Neutral only" else energies
        counts = observations.count()
        for pair in pairs:
            row = {"dataset": dataset, "pair": pair}
            expected_targets = [species for species in ADSORBATES if species not in pair] + list(TS)
            for model, folder in zip(("linear", "nonlinear"), folders):
                path = ROOT / folder / result_name(pair)
                results = pd.read_pickle(path)
                if results["TARGET"].tolist() != expected_targets:
                    raise ValueError(f"Unexpected targets in {path}")
                for metric, column in (("fit", "FIT_RMSE_EV"), ("lomo", "LOMO_RMSE_EV")):
                    row[f"{model}_{metric}"] = pooled_rmse(results, counts, column)
            rows.append(row)
    scores = pd.DataFrame(rows)
    if len(scores) != 72:
        raise ValueError("Expected 36 pairs for each dataset")
    return scores


def make_figure(scores):
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 8.0,
        "axes.linewidth": 0.75, "axes.grid": False,
        "savefig.facecolor": "white", "figure.facecolor": "white",
    })
    fig, axes = plt.subplots(1, 2, figsize=(3.8, 2.30), constrained_layout=False)
    fig.subplots_adjust(left=0.13, right=0.97, bottom=0.17, top=0.78, wspace=0.35)
    styles = {
        "Neutral only": dict(facecolors="white", edgecolors="black"),
        "All charges": dict(facecolors="#555555", edgecolors="black"),
    }
    for index, kind in enumerate(("fit", "lomo")):
        ax = axes[index]
        if kind == "fit":
            ax.set_xlim(0, 0.70)
            ax.set_ylim(0, 0.70)
        else:
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlim(0.5, 20)
            ax.set_ylim(0.5, 100)
            ax.xaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10)))
            ax.yaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10)))
            ax.xaxis.set_minor_formatter(NullFormatter())
            ax.yaxis.set_minor_formatter(NullFormatter())
        lo = max(ax.get_xlim()[0], ax.get_ylim()[0])
        hi = min(ax.get_xlim()[1], ax.get_ylim()[1])
        ax.plot([lo, hi], [lo, hi], ls="--", lw=0.8, color="#777777", zorder=1)
        for dataset in ("Neutral only", "All charges"):
            subset = scores.loc[scores.dataset.eq(dataset)]
            ax.scatter(
                subset[f"linear_{kind}"], subset[f"nonlinear_{kind}"],
                marker="o" if kind == "fit" else "s",
                s=23, linewidths=0.75, alpha=0.9, zorder=2,
                label=f"{dataset}, {'fit' if kind == 'fit' else 'LOMO'}",
                **styles[dataset],
            )
        ax.text(
            0.97, 0.04, f"({chr(ord('a') + index)})",
            transform=ax.transAxes, ha="right", va="bottom",
            fontsize=9.0, fontweight="bold",
        )
        ax.tick_params(axis="both", labelsize=8.0, direction="out", length=2.4)
        ax.set_xlabel("Linear RMSE (eV)", fontsize=8.5, labelpad=3)
        ax.set_title(
            "Fit RMSE" if kind == "fit" else "Nested LOMO RMSE",
            fontsize=8.5, fontweight="bold", pad=4,
        )
        for spine in ax.spines.values():
            spine.set_visible(True)
    fig.text(0.015, 0.52, "Nonlinear RMSE (eV)", va="center", rotation=90, fontsize=8.5)
    fig.legend(
        list(axes[0].collections), ["Neutral only", "All charges"],
        loc="upper center", bbox_to_anchor=(0.55, 1.0), ncol=2,
        frameon=False, fontsize=7.8, handletextpad=0.35, columnspacing=1.0,
    )
    fig.savefig(HERE / "FIGURE_4.PNG", dpi=1200)
    plt.close(fig)


if __name__ == "__main__":
    make_figure(load_scores())
