"""Rebuild SI Figure S22 from the energy and BEP-result PKLs at 1200 DPI.

Run from any directory with:
    python FIGURE_S22/MAKE_FIGURE_S22.PY

The neutral BEP equations are evaluated on neutral data and then on all
charges unchanged. The last row uses the equations fitted to all charges.
"""

from pathlib import Path
from runpy import run_path
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FormatStrFormatter
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent.parent
OUT = Path(__file__).resolve().parent / "FIGURE_S22.PNG"
REACTION_ENERGIES = run_path(str(ROOT / "BEP_COMMON.PY"))["reaction_energies"]
CHARGES = (-0.11, -0.05, 0.0, 0.05, 0.11)
COLORS = ("#748eba", "#b4c7e7", "#FFFFFF", "#f8cbad", "#d69163")
CHARGE_LABELS = (
    r"$-\frac{2}{18}$", r"$-\frac{1}{18}$", "0",
    r"$+\frac{1}{18}$", r"$+\frac{2}{18}$",
)
MARKERS = {"Pt": "o", "Cu": "s", "Ni": "D", "Au": "v"}
TARGETS = (
    ("HO--H", r"$\mathrm{HO{-}H}$", 2.6, (0.0, 1.0, 2.0)),
    ("O--CO", r"$\mathrm{O{-}CO}$", 4.2, (0.0, 2.0, 4.0)),
    ("H--CO2", r"$\mathrm{H{-}CO_2}$", 3.2, (0.0, 1.5, 3.0)),
    ("OC--OH", r"$\mathrm{OC{-}OH}$", 1.4, (0.0, 0.5, 1.0)),
)
COMPARISONS = (
    ("Neutral BEP to neutral", "Neutral fit to neutral"),
    ("Neutral BEP to full", "Neutral fit to full"),
    ("Full BEP to full", "Full fit to full"),
)


def coefficients(equation):
    """Read slope and intercept from a saved BEP relation."""
    match = re.fullmatch(
        r"\s*([+-]?\d+(?:\.\d+)?)\\Delta E_.*?([+-]\d+(?:\.\d+)?)\s*",
        equation,
    )
    if match is None:
        raise ValueError(f"Unrecognized BEP equation: {equation}")
    return float(match[1]), float(match[2])


def score(actual, predicted):
    actual = np.asarray(actual, dtype=float)
    predicted = np.asarray(predicted, dtype=float)
    sst = float(np.square(actual - actual.mean()).sum())
    if sst == 0:
        raise ValueError("R2 is undefined for a constant target")
    return 1.0 - float(np.square(actual - predicted).sum()) / sst, float(
        np.sqrt(np.square(actual - predicted).mean())
    )


def build_cases():
    energy_data = pd.read_pickle(ROOT / "ADSORPTION_TS_ENERGIES.PKL")
    neutral = pd.read_pickle(ROOT / "BEP_NEUTRAL/RESULTS_ALL.PKL").set_index("TARGET")
    full = pd.read_pickle(ROOT / "BEP_ALL_CHARGES/RESULTS_ALL.PKL").set_index("TARGET")
    cases = []
    equations = {}
    for target, _, _, _ in TARGETS:
        observed = energy_data.loc[energy_data[target].notna()].copy()
        reaction = REACTION_ENERGIES(observed, target)
        observed["REACTION_ENERGY"] = reaction
        for comparison, model, is_neutral in (
            (COMPARISONS[0][0], neutral, True),
            (COMPARISONS[1][0], neutral, False),
            (COMPARISONS[2][0], full, False),
        ):
            subset = observed.loc[np.isclose(observed.Charge, 0)] if is_neutral else observed
            slope, intercept = coefficients(model.at[target, "EQUATION"])
            predicted = slope * subset.REACTION_ENERGY.to_numpy(float) + intercept
            actual = subset[target].to_numpy(float)
            if comparison != COMPARISONS[1][0]:
                fitted_rmse = score(actual, predicted)[1]
                reported = float(model.at[target, "FIT_RMSE_EV"])
                # Result PKLs report fitted RMSE to three decimals.
                if abs(fitted_rmse - reported) > 0.0006:
                    raise ValueError(f"Saved BEP fit RMSE does not match {target}: {fitted_rmse}")
            equations[(comparison, target)] = r"$" + f"{slope:.2f}r{intercept:+.2f}" + r"$"
            for (_, row), estimate in zip(subset.iterrows(), predicted):
                cases.append({
                    "comparison": comparison, "target": target,
                    "metal": str(row.Metal).split("(")[0], "charge": float(row.Charge),
                    "actual_eV": float(row[target]), "predicted_eV": float(estimate),
                })
    table = pd.DataFrame(cases)
    counts = table.groupby("comparison").size().to_dict()
    if counts != {COMPARISONS[0][0]: 16, COMPARISONS[1][0]: 74, COMPARISONS[2][0]: 74}:
        raise ValueError(f"Unexpected numbers of BEP parity points: {counts}")
    return table, equations


def main():
    cases, equations = build_cases()
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "mathtext.fontset": "dejavusans",
        "axes.spines.top": True, "axes.spines.right": True,
    })
    fig = plt.figure(figsize=(13.0, 8.8), facecolor="white")
    x_positions = (0.110, 0.335, 0.560, 0.785)
    y_positions = (0.660, 0.365, 0.070)
    width, height = 0.160, 0.235

    for row, (comparison, row_label) in enumerate(COMPARISONS):
        for col, (target, target_label, limit, ticks) in enumerate(TARGETS):
            ax = fig.add_axes([x_positions[col], y_positions[row], width, height])
            panel = cases.loc[cases.comparison.eq(comparison) & cases.target.eq(target)]
            for point in panel.itertuples(index=False):
                charge_index = min(range(len(CHARGES)), key=lambda i: abs(point.charge - CHARGES[i]))
                if abs(point.charge - CHARGES[charge_index]) > 1e-8:
                    raise ValueError(f"Unexpected charge: {point.charge}")
                ax.scatter(
                    point.predicted_eV, point.actual_eV,
                    marker=MARKERS[point.metal], s=41,
                    facecolor=COLORS[charge_index], edgecolor="#282828",
                    linewidth=0.75, zorder=3,
                )
            ax.plot((0, limit), (0, limit), "--", color="#696969", lw=1, zorder=1)
            ax.set_xlim(0, limit)
            ax.set_ylim(0, limit)
            ax.set_aspect("equal", adjustable="box")
            ax.set_xticks(ticks)
            ax.set_yticks(ticks)
            ax.xaxis.set_major_formatter(FormatStrFormatter("%.2f"))
            ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
            ax.tick_params(labelsize=9.1, length=3.3, pad=2)
            for spine in ax.spines.values():
                spine.set_linewidth(0.8)
            r2, rmse = score(panel.actual_eV, panel.predicted_eV)
            metric_x = 0.965 if target == "HO--H" and row > 0 else 0.035
            ax.text(
                metric_x, 0.97,
                rf"$R^2={r2:.2f}$" + "\n" + rf"RMSE$={rmse:.2f}$ eV",
                transform=ax.transAxes,
                ha="right" if metric_x > 0.5 else "left", va="top",
                fontsize=9.1, linespacing=1.05,
            )
            ax.text(
                0.97, 0.035, f"({chr(ord('a') + 4 * row + col)})",
                transform=ax.transAxes, ha="right", va="bottom",
                fontsize=10.8, fontweight="bold",
            )
            fig.text(
                x_positions[col] + width / 2, y_positions[row] - 0.047,
                equations[(comparison, target)], ha="center", va="center", fontsize=10.5,
            )
            if row == 0:
                fig.text(
                    x_positions[col] + width / 2, 0.915, target_label,
                    ha="center", va="center", fontsize=13.0, fontweight="bold",
                )
        fig.text(0.025, y_positions[row] + height / 2, row_label,
                 ha="center", va="center", rotation=90, fontsize=11.3)

    metal_handles = [
        Line2D([], [], ls="none", marker=marker, markersize=7,
               markerfacecolor="white", markeredgecolor="#282828", label=metal)
        for metal, marker in MARKERS.items()
    ]
    charge_handles = [
        Line2D([], [], ls="none", marker="o", markersize=7,
               markerfacecolor=color, markeredgecolor="#282828", label=label)
        for color, label in zip(COLORS, CHARGE_LABELS)
    ]
    fig.legend(handles=metal_handles, ncol=4, loc="center", bbox_to_anchor=(0.55, 0.985),
               frameon=False, fontsize=11.2, columnspacing=2.4, handletextpad=0.4)
    fig.legend(handles=charge_handles, ncol=5, loc="center", bbox_to_anchor=(0.55, 0.950),
               frameon=False, fontsize=11.2, columnspacing=2.0, handletextpad=0.4)
    fig.text(0.110, 0.985, "Metals", ha="left", va="center", fontsize=11.2)
    fig.text(0.110, 0.950, "Charge", ha="left", va="center", fontsize=11.2)
    fig.text(0.070, 0.483, r"$\Delta E^{\mathrm{act}}_i$ (eV)",
             ha="center", va="center", rotation=90, fontsize=13.0)
    fig.savefig(OUT, dpi=1200, facecolor="white")
    plt.close(fig)
    print(OUT)


if __name__ == "__main__":
    main()
