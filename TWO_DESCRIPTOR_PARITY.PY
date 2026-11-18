"""Draw SI winner parity panels from repository energies and equation PKLs."""

from pathlib import Path
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.offsetbox import AnnotationBbox, HPacker, VPacker, TextArea
from matplotlib.ticker import FormatStrFormatter, MaxNLocator
import numpy as np
import pandas as pd
import sympy as sp
from sympy.parsing.sympy_parser import (
    implicit_multiplication_application, parse_expr, standard_transformations,
)

ROOT = Path(__file__).resolve().parent
ADS = ("H*", "O*", "OH*", "CH*", "CO* FCC", "CO2*", "H2O*", "HCO*", "COOH*")
TS = ("HO--H", "O--CO", "H--CO2", "OC--OH")
MARKERS = {"Pt": "o", "Cu": "s", "Ru": "^", "Ni": "D", "Au": "v"}
CHARGES = (-0.11, -0.05, 0.0, 0.05, 0.11)
COLORS = ("#748eba", "#b4c7e7", "#FFFFFF", "#f8cbad", "#d69163")
VAR_COLORS = {"x": "#7041a6", "y": "#157a52"}
ENERGY_SYMBOLS = {
    "CH*": r"\Delta E_{\mathrm{CH}^{*}}^{\mathrm{ad}}",
    "HCO*": r"\Delta E_{\mathrm{HCO}^{*}}^{\mathrm{ad}}",
    "CO2*": r"\Delta E_{\mathrm{CO}_{2}^{*}}^{\mathrm{ad}}",
    "COOH*": r"\Delta E_{\mathrm{COOH}^{*}}^{\mathrm{ad}}",
    "OH*": r"\Delta E_{\mathrm{OH}^{*}}^{\mathrm{ad}}",
    "CO* FCC": r"\Delta E_{\mathrm{CO}^{*}_{\mathrm{FCC}}}^{\mathrm{ad}}",
    "H2O*": r"\Delta E_{\mathrm{H}_{2}\mathrm{O}^{*}}^{\mathrm{ad}}",
}
TARGET_LABELS = {
    "H*": r"$\mathrm{H}^{*}$", "O*": r"$\mathrm{O}^{*}$",
    "OH*": r"$\mathrm{OH}^{*}$", "CH*": r"$\mathrm{CH}^{*}$",
    "CO* FCC": r"$\mathrm{CO}_{\mathrm{FCC}}^{*}$",
    "CO2*": r"$\mathrm{CO}_{2}^{*}$",
    "H2O*": r"$\mathrm{H}_{2}\mathrm{O}^{*}$",
    "HCO*": r"$\mathrm{HCO}^{*}$", "COOH*": r"$\mathrm{COOH}^{*}$",
    "HO--H": "HO–H TS", "O--CO": "O–CO TS",
    "H--CO2": "H–CO₂ TS", "OC--OH": "OC–OH TS",
}

# These are display labels from the SI figures. Predictions always use the
# equations stored in the result PKLs, not these rounded labels.
NEUTRAL_NONLINEAR_LABELS = {
    "H*": "1.35x-1.35y+4.13-4.64y/x",
    "O*": "1.03x+0.75y+2.84",
    "CH*": "1.33x+2.32-6.22y/x",
    "CO* FCC": "-3.96-5.08/y",
    "CO2*": "0.002x-0.021y-0.055",
    "H2O*": "-0.056xy-0.056x-0.080",
    "HCO*": "0.59x+1.71x/y-2.64",
    "HO--H": "-0.82xy-1.81x-1.81y-3.14",
    "O--CO": "6.55x-6.55y+23.82-22.99y/x",
    "H--CO2": "4.38x-4.38y+16.14-15.41y/x",
    "OC--OH": "0.81xy+1.96x+1.96y+5.47",
}
FULL_NONLINEAR_LABELS = {
    "H*": "0.12x-0.18/x-0.40",
    "O*": "0.51x+1.13y/x+0.64y-1.16",
    "OH*": "3.31y+1.91y/x+0.90xy-2.83",
    "CH*": "1.42x+0.19y-3.97",
    "CO2*": "-0.27xy-0.16x-0.09",
    "HCO*": "2.54x/(x+y)-0.18y+1.12x-2.35",
    "COOH*": "0.64x+0.58y-1.01",
    "HO--H": "0.51x(1-y)-y(2.17+1.27/x)+1.73",
    "O--CO": "-0.78x+(1.23y-1.41)/x-1.51",
    "H--CO2": "-0.44/y-1.56y+2.03y/x-1.60",
    "OC--OH": "-0.78/x-0.49xy-1.04(x+y)-1.70",
}
MODELS = {
    "S17": {"pair": ("CH*", "HCO*"), "file": "LINEAR_NEUTRAL/RESULTS_CH_HCO.PKL", "neutral_fit": True},
    "S18": {"pair": ("CO2*", "COOH*"), "file": "LINEAR_ALL_CHARGES/RESULTS_CO2_COOH.PKL", "neutral_fit": False},
    "S19": {"pair": ("OH*", "COOH*"), "file": "NONLINEAR_NEUTRAL/RESULTS_OH_COOH.PKL", "neutral_fit": True,
            "labels": NEUTRAL_NONLINEAR_LABELS, "metric_y": {"H--CO2": 0.77}, "panel_y": {"CO2*": 0.10}},
    "S20": {"pair": ("CO* FCC", "H2O*"), "file": "NONLINEAR_ALL_CHARGES/RESULTS_CO_FCC_H2O.PKL", "neutral_fit": False,
            "labels": FULL_NONLINEAR_LABELS},
}


def parse_equation(latex, pair):
    """Convert the restricted saved +,-,*,/ equation language to SymPy."""
    cleaned = latex
    for i, descriptor in enumerate(pair):
        cleaned = cleaned.replace(ENERGY_SYMBOLS[descriptor], f"X{i}")
    if "\\Delta" in cleaned:
        raise ValueError(f"Unrecognized descriptor in {latex}")
    cleaned = cleaned.replace(r"\left(", "(").replace(r"\right)", ")")
    fraction = re.compile(r"\\frac\{([^{}]+)\}\{([^{}]+)\}")
    while r"\frac" in cleaned:
        new = fraction.sub(lambda match: f"(({match[1]})/({match[2]}))", cleaned)
        if new == cleaned:
            raise ValueError(f"Unrecognized fraction in {latex}")
        cleaned = new
    cleaned = cleaned.replace("\\", "")
    x0, x1 = sp.symbols("X0 X1")
    return parse_expr(cleaned, local_dict={"X0": x0, "X1": x1},
                      transformations=standard_transformations + (implicit_multiplication_application,))


def linear_label(expression):
    x0, x1 = sp.symbols("X0 X1")
    polynomial = sp.Poly(expression, x0, x1)
    terms = ((float(polynomial.coeff_monomial(x0)), "x"),
             (float(polynomial.coeff_monomial(x1)), "y"),
             (float(polynomial.coeff_monomial(1)), ""))
    result = ""
    for coefficient, variable in terms:
        if abs(coefficient) < 1e-9 and variable:
            continue
        sign = "-" if coefficient < 0 else ("+" if result else "")
        result += f"{sign}{abs(coefficient):.2f}{variable}"
    return result


def calculate_cases(figure_id):
    config = MODELS[figure_id]
    pair = config["pair"]
    data = pd.read_pickle(ROOT / "ADSORPTION_TS_ENERGIES.PKL")
    fits = pd.read_pickle(ROOT / config["file"]).set_index("TARGET")
    targets = [target for target in (*ADS, *TS) if target not in pair]
    if len(targets) != 11 or set(targets) != set(fits.index):
        raise ValueError("The winner PKL does not have the expected 11 targets")
    cases, metrics, labels = {}, {}, {}
    x0, x1 = sp.symbols("X0 X1")
    for target in targets:
        equation = parse_equation(fits.at[target, "EQUATION"], pair)
        function = sp.lambdify((x0, x1), equation, modules="numpy")
        rows = data.loc[data[target].notna() & data[list(pair)].notna().all(axis=1)].copy()
        actual = rows[target].to_numpy(float)
        values = np.asarray(function(rows[pair[0]].to_numpy(float), rows[pair[1]].to_numpy(float)), dtype=float)
        predicted = np.broadcast_to(values, actual.shape)
        if not np.isfinite(predicted).all():
            raise ValueError(f"Non-finite predictions for {target}")
        fit_mask = np.isclose(rows.Charge.to_numpy(float), 0) if config["neutral_fit"] else np.ones(len(rows), bool)
        fit_rmse = float(np.sqrt(np.mean((predicted[fit_mask] - actual[fit_mask]) ** 2)))
        if abs(fit_rmse - float(fits.at[target, "FIT_RMSE_EV"])) > 0.0001:
            raise ValueError(f"Saved fit RMSE does not match predictions for {target}: {fit_rmse}")
        residual = predicted - actual
        sst = float(np.sum((actual - actual.mean()) ** 2))
        metrics[target] = (1.0 - float(np.sum(residual ** 2)) / sst,
                           float(np.sqrt(np.mean(residual ** 2))))
        rows["predicted_eV"] = predicted
        cases[target] = rows
        labels[target] = config.get("labels", {}).get(target) or linear_label(equation)
    return pair, targets, cases, metrics, labels


def colored_equation(expression):
    lines = []
    for line in expression.split("\n"):
        pieces = []
        for token in re.split(r"([xy])", line):
            if not token:
                continue
            math = rf"$\mathbf{{{token}}}$" if token in VAR_COLORS else "$" + token + "$"
            pieces.append(TextArea(math, textprops={
                "fontsize": 10.0, "color": VAR_COLORS.get(token, "#202020"),
                "fontfamily": "DejaVu Sans", "math_fontfamily": "dejavusans",
            }))
        lines.append(HPacker(children=pieces, align="baseline", pad=0, sep=0.25))
    return lines[0] if len(lines) == 1 else VPacker(children=lines, align="center", pad=0, sep=0.7)


def make_figure(figure_id):
    pair, targets, cases, metrics, labels = calculate_cases(figure_id)
    config = MODELS[figure_id]
    color_for = dict(zip(CHARGES, COLORS))
    plt.rcParams.update({"font.family": "DejaVu Sans", "mathtext.fontset": "dejavusans"})
    fig = plt.figure(figsize=(14, 10.8))
    slots = [(0, 0), (0, 1), (0, 2), (0, 3),
             (1, 0), (1, 1), (1, 2),
             (2, 0), (2, 1), (2, 2), (2, 3)]
    x0, dx, width = 0.095, 0.220, 0.205
    y0, dy, height = 0.635, 0.268, 0.200
    for index, (target, (r, c)) in enumerate(zip(targets, slots)):
        ax = fig.add_axes([x0 + c * dx, y0 - r * dy, width, height])
        rows = cases[target]
        values = np.r_[rows.predicted_eV.to_numpy(float), rows[target].to_numpy(float)]
        span = float(np.ptp(values))
        pad = max(0.025, span * 0.09)
        lo, hi = float(values.min() - pad), float(values.max() + pad)
        ax.plot([lo, hi], [lo, hi], ls="--", lw=1.0, color="#696969", zorder=1)
        for _, point in rows.iterrows():
            metal = str(point.Metal).split("(")[0]
            charge = min(CHARGES, key=lambda value: abs(value - point.Charge))
            if abs(charge - point.Charge) > 1e-8:
                raise ValueError(f"Unexpected charge {point.Charge}")
            ax.scatter(point.predicted_eV, point[target],
                       marker=MARKERS[metal], s=43, facecolor=color_for[charge],
                       edgecolor="#282828", linewidth=0.75, zorder=2)
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_aspect("equal", adjustable="box")
        ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
        ax.xaxis.set_major_formatter(FormatStrFormatter("%.2f"))
        ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
        ax.tick_params(labelsize=10.2, length=3.5, pad=2)
        for spine in ax.spines.values():
            spine.set_linewidth(0.8)
        ax.set_title(TARGET_LABELS[target], fontsize=13.2, fontweight="bold", pad=4)
        r2, rmse = metrics[target]
        ax.text(0.035, config.get("metric_y", {}).get(target, 0.97),
                rf"$R^2_{{\mathrm{{all}}}}={r2:.2f}$" + "\n" +
                rf"RMSE$={rmse:.2f}$ eV",
                transform=ax.transAxes, va="top", ha="left", fontsize=9.5)
        ax.text(0.975, config.get("panel_y", {}).get(target, 0.035),
                f"({chr(ord('a') + index)})", transform=ax.transAxes,
                ha="right", va="bottom", fontsize=11.5, fontweight="bold")
        ax.add_artist(AnnotationBbox(
            colored_equation(labels[target]), (0.5, -0.15), xycoords=ax.transAxes,
            frameon=False, box_alignment=(0.5, 0.5), annotation_clip=False))
    metal_handles = [Line2D([], [], ls="none", marker=shape, markersize=7,
                            markerfacecolor="white", markeredgecolor="#282828", label=metal)
                     for metal, shape in MARKERS.items()]
    charge_labels = (r"$-\frac{2}{18}$", r"$-\frac{1}{18}$", "0",
                     r"$+\frac{1}{18}$", r"$+\frac{2}{18}$")
    charge_handles = [Line2D([], [], ls="none", marker="o", markersize=7,
                             markerfacecolor=color_for[charge], markeredgecolor="#282828", label=label)
                      for charge, label in zip(CHARGES, charge_labels)]
    fig.legend(handles=metal_handles, ncol=5, loc="center", bbox_to_anchor=(0.51, 0.97),
               frameon=False, fontsize=11.5, columnspacing=2.4, handletextpad=0.4)
    fig.legend(handles=charge_handles, ncol=5, loc="center", bbox_to_anchor=(0.51, 0.928),
               frameon=False, fontsize=11.5, columnspacing=2.2, handletextpad=0.4)
    fig.text(0.08, 0.97, "Metals", ha="left", va="center", fontsize=11.5)
    fig.text(0.08, 0.928, "Charge", ha="left", va="center", fontsize=11.5)
    fig.text(0.36, 0.889, r"$\mathbf{x}=" + ENERGY_SYMBOLS[pair[0]] + "$",
             color=VAR_COLORS["x"], ha="center", va="center", fontsize=12.2)
    fig.text(0.67, 0.889, r"$\mathbf{y}=" + ENERGY_SYMBOLS[pair[1]] + "$",
             color=VAR_COLORS["y"], ha="center", va="center", fontsize=12.2)
    fig.text(0.040, 0.605, r"$\Delta E^{\mathrm{ad}}_{s*}$ (eV)", rotation=90,
             ha="center", va="center", fontsize=13.5)
    fig.text(0.040, 0.198, r"$\Delta E^{\mathrm{act}}_i$ (eV)", rotation=90,
             ha="center", va="center", fontsize=13.5)
    fig.text(0.52, 0.030, "Scaling-relation prediction (eV)",
             ha="center", va="bottom", fontsize=13.5)
    output = ROOT / f"FIGURE_{figure_id}" / f"FIGURE_{figure_id}.PNG"
    output.parent.mkdir(exist_ok=True)
    fig.savefig(output, dpi=1200, facecolor="white")
    plt.close(fig)
    return output
