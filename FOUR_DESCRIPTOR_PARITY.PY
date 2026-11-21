"""Rebuild the SI four-descriptor parity figures from repository PKLs."""

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
CHARGES = (-0.11, -0.05, 0.0, 0.05, 0.11)
COLORS = ("#748eba", "#b4c7e7", "#FFFFFF", "#f8cbad", "#d69163")
MARKERS = {"Pt": "o", "Cu": "s", "Ru": "^", "Ni": "D", "Au": "v"}
CHARGE_LABELS = (
    r"$-\frac{2}{18}$", r"$-\frac{1}{18}$", "0",
    r"$+\frac{1}{18}$", r"$+\frac{2}{18}$",
)
VARIABLE_COLORS = {"x": "#7041a6", "y": "#157a52", "z": "#bf6b1e", "w": "#27609b"}
ENERGY_TERMS = {
    "CH*": r"\Delta E_{\mathrm{CH}^{*}}^{\mathrm{ad}}",
    "HCO*": r"\Delta E_{\mathrm{HCO}^{*}}^{\mathrm{ad}}",
    "H2O*": r"\Delta E_{\mathrm{H}_{2}\mathrm{O}^{*}}^{\mathrm{ad}}",
    "CO2*": r"\Delta E_{\mathrm{CO}_{2}^{*}}^{\mathrm{ad}}",
    "O*": r"\Delta E_{\mathrm{O}^{*}}^{\mathrm{ad}}",
}
LABEL_SYMBOLS = {
    "CH*": r"\Delta E^{\mathrm{ad}}_{\mathrm{CH}^{*}}",
    "HCO*": r"\Delta E^{\mathrm{ad}}_{\mathrm{HCO}^{*}}",
    "H2O*": r"\Delta E^{\mathrm{ad}}_{\mathrm{H_2O}^{*}}",
    "CO2*": r"\Delta E^{\mathrm{ad}}_{\mathrm{CO_2}^{*}}",
    "O*": r"\Delta E^{\mathrm{ad}}_{\mathrm{O}^{*}}",
}
TARGET_LABELS = {
    "H*": r"$\mathrm{H}^{*}$", "O*": r"$\mathrm{O}^{*}$",
    "OH*": r"$\mathrm{OH}^{*}$",
    "CO* FCC": r"$\mathrm{CO}^{*}_{\mathrm{FCC}}$",
    "CO2*": r"$\mathrm{CO}_{2}^{*}$",
    "COOH*": r"$\mathrm{COOH}^{*}$",
    "HO--H": r"$\mathrm{HO{-}H}$ TS",
    "O--CO": r"$\mathrm{O{-}CO}$ TS",
    "H--CO2": r"$\mathrm{H{-}CO}_{2}$ TS",
    "OC--OH": r"$\mathrm{OC{-}OH}$ TS",
}
LINEAR_TARGETS = (
    "H*", "O*", "OH*", "CO* FCC", "COOH*",
    "HO--H", "O--CO", "H--CO2", "OC--OH",
)
NONLINEAR_TARGETS = (
    "OH*", "CO* FCC", "CO2*", "COOH*",
    "O--CO", "H--CO2", "OC--OH",
)
DISPLAY_NONLINEAR = {
    "OH*": "-0.29x+0.77w+0.73z-3.04",
    "CO* FCC": "0.08x/w-0.32z\n+0.88y-0.25",
    "CO2*": "-0.17zy-0.06x-0.22",
    "COOH*": "0.44x+0.50z+0.74",
    "O--CO": "1.55yw+3.45w\n+1.80x+12.84",
    "H--CO2": "-3.18w/y-1.30(z+w)\n+2.24y+5.65",
    "OC--OH": "-1.36w/y-0.57w+0.26x+2.49",
}
CONFIG = {
    "S21": {
        "descriptors": ("CH*", "HCO*", "H2O*", "CO2*"),
        "file": "LINEAR_FOUR_DESCRIPTORS/RESULTS_CH_HCO_H2O_CO2.PKL",
        "targets": LINEAR_TARGETS,
    },
    "S23": {
        "descriptors": ("CH*", "HCO*", "H2O*", "O*"),
        "file": "NONLINEAR_FOUR_DESCRIPTORS/RESULTS_CH_HCO_O_H2O.PKL",
        "targets": NONLINEAR_TARGETS,
    },
}


def parse_equation(latex, descriptors):
    """Read saved arithmetic LaTex equations with four adsorption energies."""
    text = latex
    for i, descriptor in enumerate(descriptors):
        text = text.replace(ENERGY_TERMS[descriptor], f"X{i}")
    if r"\Delta" in text:
        raise ValueError(f"Unrecognized descriptor in equation: {latex}")
    text = text.replace(r"\left(", "(").replace(r"\right)", ")")
    fraction = re.compile(r"\\frac\{([^{}]+)\}\{([^{}]+)\}")
    while r"\frac" in text:
        replaced = fraction.sub(lambda match: f"(({match[1]})/({match[2]}))", text)
        if replaced == text:
            raise ValueError(f"Unrecognized fraction in equation: {latex}")
        text = replaced
    text = re.sub(r"(X\d)(?=X\d)", r"\1*", text)
    if "\\" in text:
        raise ValueError(f"Unrecognized LaTex notation in equation: {latex}")
    symbols = sp.symbols("X0:4")
    return parse_expr(
        text, local_dict={str(symbol): symbol for symbol in symbols},
        transformations=standard_transformations + (implicit_multiplication_application,),
    )


def load_cases(figure_id):
    config = CONFIG[figure_id]
    descriptors = config["descriptors"]
    models = pd.read_pickle(ROOT / config["file"]).set_index("TARGET")
    data = pd.read_pickle(ROOT / "ADSORPTION_TS_ENERGIES.PKL")
    if len(models) != 9 or set(models.index).intersection(descriptors):
        raise ValueError("Four-descriptor PKL has an unexpected target set")
    if not set(config["targets"]).issubset(models.index):
        raise ValueError("Four-descriptor PKL is missing a plotted target")
    cases, metrics, equations = {}, {}, {}
    symbols = sp.symbols("X0:4")
    for target in config["targets"]:
        model = models.loc[target]
        equation = parse_equation(model.EQUATION, descriptors)
        observed = data.loc[data[target].notna() & data[list(descriptors)].notna().all(axis=1)].copy()
        function = sp.lambdify(symbols, equation, modules="numpy")
        raw = function(*[observed[name].to_numpy(float) for name in descriptors])
        predicted = np.broadcast_to(np.asarray(raw, dtype=float), (len(observed),))
        actual = observed[target].to_numpy(float)
        if not np.isfinite(predicted).all() or len(observed) != int(model.N):
            raise ValueError(f"Invalid predictions or observation count for {target}")
        rmse = float(np.sqrt(np.mean(np.square(predicted - actual))))
        r2 = 1.0 - float(np.square(predicted - actual).sum()) / float(np.square(actual - actual.mean()).sum())
        tolerance = 0.00006 if figure_id == "S21" else 1e-9
        if abs(rmse - float(model.FIT_RMSE_EV)) > tolerance:
            raise ValueError(f"Saved fit RMSE does not match {target}: {rmse}")
        if abs(r2 - float(model.FIT_R2)) > (0.00006 if figure_id == "S21" else 1e-9):
            raise ValueError(f"Saved fit R2 does not match {target}: {r2}")
        observed["PREDICTED_EV"] = predicted
        cases[target] = observed
        metrics[target] = (r2, rmse, float(model.LOMO_RMSE_EV))
        equations[target] = equation
    return cases, metrics, equations


def colored_equation(expression):
    lines = []
    for line in expression.split("\n"):
        pieces = []
        for token in re.split(r"([xyzw])", line):
            if not token:
                continue
            math = rf"$\mathbf{{{token}}}$" if token in VARIABLE_COLORS else "$" + token + "$"
            pieces.append(TextArea(math, textprops={
                "fontsize": 10.0,
                "color": VARIABLE_COLORS.get(token, "#202020"),
                "fontfamily": "DejaVu Sans", "math_fontfamily": "dejavusans",
            }))
        lines.append(HPacker(children=pieces, align="baseline", pad=0, sep=0.20))
    return lines[0] if len(lines) == 1 else VPacker(children=lines, align="center", pad=0, sep=0.6)


def linear_display(expression, descriptors, target):
    symbols = sp.symbols("X0:4")
    polynomial = sp.Poly(expression, *symbols)
    out = ""
    for symbol, variable in zip(symbols, "xyzw"):
        value = float(polynomial.coeff_monomial(symbol))
        if abs(value) < 1e-9:
            continue
        sign = "-" if value < 0 else ("+" if out else "")
        out += f"{sign}{abs(value):.2f}{variable}"
        if target == "O--CO" and variable == "y":
            out += "\n"
    intercept = float(polynomial.coeff_monomial(1))
    out += ("-" if intercept < 0 else "+") + f"{abs(intercept):.2f}"
    return out


def legends(fig, descriptors, *, top_metal, top_charge, descriptor_y, metal_label_x):
    color_for = dict(zip(CHARGES, COLORS))
    metal_handles = [Line2D([], [], ls="none", marker=shape, markersize=7,
                            markerfacecolor="white", markeredgecolor="#282828", label=metal)
                     for metal, shape in MARKERS.items()]
    charge_handles = [Line2D([], [], ls="none", marker="o", markersize=7,
                             markerfacecolor=color_for[charge], markeredgecolor="#282828", label=label)
                      for charge, label in zip(CHARGES, CHARGE_LABELS)]
    fig.legend(handles=metal_handles, ncol=5, loc="center",
               bbox_to_anchor=(0.53, top_metal), frameon=False,
               fontsize=11.2, columnspacing=2.2, handletextpad=0.4)
    fig.legend(handles=charge_handles, ncol=5, loc="center",
               bbox_to_anchor=(0.53, top_charge), frameon=False,
               fontsize=11.2, columnspacing=2.0, handletextpad=0.4)
    fig.text(metal_label_x, top_metal, "Metals", ha="left", va="center", fontsize=11.2)
    fig.text(metal_label_x, top_charge, "Charge", ha="left", va="center", fontsize=11.2)
    for x, descriptor, variable in zip((0.15, 0.38, 0.62, 0.85), descriptors, "xyzw"):
        fig.text(x, descriptor_y,
                 r"$\mathbf{" + variable + r"}=" + LABEL_SYMBOLS[descriptor] + "$",
                 color=VARIABLE_COLORS[variable], ha="center", va="center", fontsize=11.7)
    return color_for


def plot_s21():
    cases, metrics, equations = load_cases("S21")
    descriptors = CONFIG["S21"]["descriptors"]
    plt.rcParams.update({"font.family": "DejaVu Sans", "mathtext.fontset": "dejavusans"})
    fig = plt.figure(figsize=(14.0, 7.7))
    colors = legends(fig, descriptors, top_metal=0.975, top_charge=0.931,
                     descriptor_y=0.855, metal_label_x=0.065)
    for index, target in enumerate(LINEAR_TARGETS):
        row = 0 if index < 5 else 1
        col = index if row == 0 else index - 5
        x = 0.070 + 0.190 * col if row == 0 else 0.095 + 0.220 * col
        y = 0.520 if row == 0 else 0.130
        ax = fig.add_axes([x, y, 0.150, 0.285])
        sub = cases[target]
        values = np.r_[sub.PREDICTED_EV.to_numpy(float), sub[target].to_numpy(float)]
        span = float(np.ptp(values))
        upper_pad = 0.45 if target == "H--CO2" else 0.23 if target == "OC--OH" else 0.12
        pad = max(0.04, upper_pad * span)
        lo, hi = float(values.min() - pad), float(values.max() + pad)
        ax.plot([lo, hi], [lo, hi], ls="--", lw=1.0, color="#696969", zorder=1)
        for _, point in sub.iterrows():
            metal = str(point.Metal).split("(")[0]
            charge = min(CHARGES, key=lambda q: abs(q - point.Charge))
            ax.scatter(point.PREDICTED_EV, point[target], marker=MARKERS[metal], s=42,
                       facecolor=colors[charge], edgecolor="#282828", linewidth=0.75, zorder=2)
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_aspect("equal", adjustable="box")
        ax.xaxis.set_major_locator(MaxNLocator(nbins=3))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=3))
        ax.xaxis.set_major_formatter(FormatStrFormatter("%.2f"))
        ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
        ax.tick_params(labelsize=9.3, length=3.3, pad=2)
        for spine in ax.spines.values():
            spine.set_linewidth(0.8)
        ax.set_title(TARGET_LABELS[target], fontsize=12.8, fontweight="bold", pad=4)
        r2, rmse, lomo = metrics[target]
        ax.text(0.035, 0.97, rf"$R^2={r2:.2f}$" + "\n" +
                rf"RMSE$={rmse:.2f}$ eV" + "\n" +
                rf"LOMO RMSE$={lomo:.2f}$ eV",
                transform=ax.transAxes, va="top", ha="left", fontsize=8.7, linespacing=1.03)
        ax.text(0.97, 0.035, f"({chr(ord('a') + index)})",
                transform=ax.transAxes, va="bottom", ha="right", fontsize=11.0, fontweight="bold")
        ax.add_artist(AnnotationBbox(
            colored_equation(linear_display(equations[target], descriptors, target)),
            (0.5, -0.155), xycoords=ax.transAxes, frameon=False,
            box_alignment=(0.5, 0.5), annotation_clip=False))
    fig.text(0.022, 0.662, r"$\Delta E^{\mathrm{ad}}_{s*}$ (eV)", rotation=90,
             ha="center", va="center", fontsize=13.0)
    fig.text(0.022, 0.272, r"$\Delta E^{\mathrm{act}}_i$ (eV)", rotation=90,
             ha="center", va="center", fontsize=13.0)
    fig.text(0.51, 0.010, "Scaling-relation prediction (eV)",
             ha="center", va="bottom", fontsize=13.0)
    output = ROOT / "FIGURE_S21/FIGURE_S21.PNG"
    output.parent.mkdir(exist_ok=True)
    fig.savefig(output, dpi=1200, facecolor="white")
    plt.close(fig)
    return output


def plot_s23():
    cases, metrics, _ = load_cases("S23")
    descriptors = CONFIG["S23"]["descriptors"]
    plt.rcParams.update({"font.family": "DejaVu Sans", "mathtext.fontset": "dejavusans",
                         "axes.spines.top": True, "axes.spines.right": True})
    fig = plt.figure(figsize=(13.0, 8.0), facecolor="white")
    colors = legends(fig, descriptors, top_metal=0.980, top_charge=0.940,
                     descriptor_y=0.887, metal_label_x=0.065)
    width, height = 0.170, 0.310
    positions = ((0.075, 0.530), (0.305, 0.530), (0.535, 0.530),
                 (0.765, 0.530), (0.180, 0.105), (0.415, 0.105),
                 (0.650, 0.105))
    for index, (target, (left, bottom)) in enumerate(zip(NONLINEAR_TARGETS, positions)):
        ax = fig.add_axes([left, bottom, width, height])
        sub = cases[target]
        values = np.r_[sub.PREDICTED_EV.to_numpy(float), sub[target].to_numpy(float)]
        pad = max(0.05, float(np.ptp(values)) * 0.12)
        lo, hi = float(values.min() - pad), float(values.max() + pad)
        ax.plot((lo, hi), (lo, hi), ls="--", lw=1.0, color="#696969", zorder=1)
        for _, point in sub.iterrows():
            metal = str(point.Metal).split("(")[0]
            charge = min(CHARGES, key=lambda q: abs(q - point.Charge))
            ax.scatter(point.PREDICTED_EV, point[target], marker=MARKERS[metal], s=43,
                       facecolor=colors[charge], edgecolor="#282828", linewidth=0.75, zorder=2)
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_aspect("equal", adjustable="box")
        ax.xaxis.set_major_locator(MaxNLocator(nbins=3))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=3))
        ax.xaxis.set_major_formatter(FormatStrFormatter("%.2f"))
        ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
        ax.tick_params(labelsize=9.3, length=3.3, pad=2)
        for spine in ax.spines.values():
            spine.set_linewidth(0.8)
        ax.set_title(TARGET_LABELS[target], fontsize=12.8, fontweight="bold", pad=4)
        r2, rmse, lomo = metrics[target]
        ax.text(0.035, 0.97, rf"$R^2={r2:.2f}$" + "\n" +
                rf"RMSE$={rmse:.2f}$ eV" + "\n" +
                rf"LOMO RMSE$={lomo:.2f}$ eV",
                transform=ax.transAxes, va="top", ha="left", fontsize=8.8, linespacing=1.03)
        ax.text(0.97, 0.035, f"({chr(ord('a') + index)})",
                transform=ax.transAxes, va="bottom", ha="right", fontsize=11.0, fontweight="bold")
        ax.add_artist(AnnotationBbox(
            colored_equation(DISPLAY_NONLINEAR[target]), (0.5, -0.16),
            xycoords=ax.transAxes, frameon=False,
            box_alignment=(0.5, 0.5), annotation_clip=False))
    fig.text(0.027, 0.685, r"$\Delta E^{\mathrm{ad}}_{s*}$ (eV)", rotation=90,
             ha="center", va="center", fontsize=13.0)
    fig.text(0.027, 0.260, r"$\Delta E^{\mathrm{act}}_i$ (eV)", rotation=90,
             ha="center", va="center", fontsize=13.0)
    fig.text(0.505, 0.006, "Scaling-relation prediction (eV)",
             ha="center", va="bottom", fontsize=13.0)
    output = ROOT / "FIGURE_S23/FIGURE_S23.PNG"
    output.parent.mkdir(exist_ok=True)
    fig.savefig(output, dpi=1200, facecolor="white")
    plt.close(fig)
    return output
