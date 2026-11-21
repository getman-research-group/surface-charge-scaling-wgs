"""Recreate manuscript Figure 5 from the repository energy and SR result pickles."""

from pathlib import Path
import re
from runpy import run_path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, FormatStrFormatter, MultipleLocator
from matplotlib.offsetbox import AnnotationBbox, HPacker, TextArea, VPacker
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
parse_equation = run_path(str(ROOT / "FOUR_DESCRIPTOR_PARITY.PY"))["parse_equation"]
DESCRIPTORS = ("CH*", "HCO*", "H2O*", "O*")
EXTRA_PAIRS = (
    ("H*", "O*"), ("H*", "CO* FCC"), ("H*", "CO2*"),
    ("H*", "H2O*"), ("O*", "H2O*"), ("O*", "CO* FCC"),
    ("O*", "CO2*"), ("CO* FCC", "H2O*"),
    ("CO* FCC", "CO2*"), ("H2O*", "CO2*"),
)
TS = ("HO--H", "O--CO", "H--CO2", "OC--OH")
CHARGES = (-0.11, -0.05, 0.0, 0.05, 0.11)
CHARGE_COLORS = ("#748eba", "#b4c7e7", "#FFFFFF", "#f8cbad", "#d69163")
MARKERS = {"Pt": "o", "Cu": "s", "Ru": "^", "Ni": "D", "Au": "v"}
VARIABLE_COLORS = {"x": "#7041a6", "y": "#157a52",
                   "z": "#007c91", "w": "#a45116"}
FONT_SCALE = 1.30  # Match the effective printed type size of Figure 3.


def fs(points):
    return points * FONT_SCALE
ENERGY_SYMBOLS = {
    "H*": r"\Delta E_{\mathrm{H}^{*}}^{\mathrm{ad}}",
    "HO--H": r"\Delta E_{\mathrm{HO{-}H}}^{\mathrm{act}}",
}


def colored_equation(expression):
    lines = []
    for line in expression.splitlines():
        pieces = []
        for token in re.split(r"([xyzw])", line):
            if not token:
                continue
            math = rf"$\mathbf{{{token}}}$" if token in VARIABLE_COLORS else f"${token}$"
            pieces.append(TextArea(math, textprops={
                "fontsize": fs(11.5), "color": VARIABLE_COLORS.get(token, "#202020"),
                "fontfamily": "DejaVu Sans", "math_fontfamily": "dejavusans"}))
        lines.append(HPacker(children=pieces, align="baseline", pad=0, sep=0.4))
    return lines[0] if len(lines) == 1 else VPacker(children=lines, align="center", pad=0, sep=0.8)


def matched_data():
    def result_path(family, pair):
        names = ("CH", "HCO", *(name.replace("*", "").replace(" ", "_") for name in pair))
        return ROOT / family / ("RESULTS_" + "_".join(names) + ".PKL")

    def pooled_rmse(data, column):
        return float(np.sqrt(np.sum(data.N * data[column] ** 2) / data.N.sum()))

    records = []
    for number, pair in enumerate(EXTRA_PAIRS, start=1):
        linear = pd.read_pickle(result_path("LINEAR_FOUR_DESCRIPTORS", pair))
        nonlinear = pd.read_pickle(result_path("NONLINEAR_FOUR_DESCRIPTORS", pair))
        if len(linear) != 9 or len(nonlinear) != 9:
            raise ValueError(f"Unexpected number of targets for {pair}")
        for scope in ("pooled", "transition_states"):
            l = linear if scope == "pooled" else linear.loc[linear.TARGET.isin(TS)]
            n = nonlinear if scope == "pooled" else nonlinear.loc[nonlinear.TARGET.isin(TS)]
            if l.N.sum() != n.N.sum() or l.N.sum() != (199 if scope == "pooled" else 74):
                raise ValueError(f"Mismatched observation count for {pair}, {scope}")
            records.append({
                "set_number": number,
                "set_key": "+".join(("CH*", "HCO*", *pair)),
                "scope": scope,
                "linear_fit_RMSE_eV": pooled_rmse(l, "FIT_RMSE_EV"),
                "linear_nested_RMSE_eV": pooled_rmse(l, "LOMO_RMSE_EV"),
                "nonlinear_fit_RMSE_eV": pooled_rmse(n, "FIT_RMSE_EV"),
                "nonlinear_nested_RMSE_eV": pooled_rmse(n, "LOMO_RMSE_EV"),
            })
    return pd.DataFrame.from_records(records)


def parity_data():
    from sympy import lambdify, symbols

    d = pd.read_pickle(ROOT / "ADSORPTION_TS_ENERGIES.PKL")
    models = pd.read_pickle(
        ROOT / "NONLINEAR_FOUR_DESCRIPTORS/RESULTS_CH_HCO_O_H2O.PKL"
    ).set_index("TARGET")
    frames = []
    for species in ("H*", "HO--H"):
        observed = d.loc[d[species].notna() & d[list(DESCRIPTORS)].notna().all(axis=1)]
        model = models.loc[species]
        expr = parse_equation(model.EQUATION, DESCRIPTORS)
        function = lambdify(symbols("X0:4"), expr, modules="numpy")
        predicted = np.asarray(function(*[observed[name].to_numpy(float)
                                          for name in DESCRIPTORS]), dtype=float)
        predicted = np.broadcast_to(predicted, (len(observed),))
        actual = observed[species].to_numpy(float)
        rmse = float(np.sqrt(np.mean((predicted - actual) ** 2)))
        if len(observed) != int(model.N) or abs(rmse - float(model.FIT_RMSE_EV)) > 1e-8:
            raise ValueError(f"Figure 5 parity values disagree with saved model for {species}")
        frames.append(pd.DataFrame({
            "species": species, "metal": observed.Metal.to_numpy(),
            "charge": observed.Charge.to_numpy(), "actual_eV": actual,
            "predicted_eV": predicted,
        }))
    return pd.concat(frames, ignore_index=True)


def main():
    matched = matched_data()
    parity = parity_data()
    metrics = pd.read_pickle(ROOT / "NONLINEAR_FOUR_DESCRIPTORS/RESULTS_CH_HCO_O_H2O.PKL")
    metrics = metrics.set_index("TARGET")

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": fs(11),
                         "axes.linewidth": 0.8, "axes.grid": False})
    fig, axes = plt.subplots(2, 2, figsize=(7.4, 8.8))
    plot_width = 0.35
    plot_height = plot_width * 7.4 / 8.8
    left_edges = (0.14, 0.60)
    bottom_edges = (0.535, 0.11)
    for row in range(2):
        for col in range(2):
            axes[row, col].set_position([
                left_edges[col], bottom_edges[row], plot_width, plot_height])
    color = dict(zip(CHARGES, CHARGE_COLORS))

    for j, (species, title, equation) in enumerate([
        ("H*", r"$\mathbf{H}^{*}$", "0.26x+0.77w+0.10xw+1.38"),
        ("HO--H", r"$\mathbf{HO{-}H}$ TS",
         "-0.09x/w-0.65z\n-0.22yw+2.00"),
    ]):
        ax = axes[0, j]
        subset = parity.loc[parity.species.eq(species)]
        values = np.r_[subset.actual_eV, subset.predicted_eV]
        lo = np.floor((values.min() - 0.10) * 5) / 5
        hi = np.ceil((values.max() + 0.10) * 5) / 5
        ax.plot([lo, hi], [lo, hi], "--", lw=1, color="0.45", zorder=1)
        for point in subset.itertuples(index=False):
            charge = min(CHARGES, key=lambda q: abs(q - point.charge))
            ax.scatter(point.predicted_eV, point.actual_eV,
                       s=40, marker=MARKERS[str(point.metal).split("(")[0]],
                       facecolor=color[charge], edgecolor="0.2", linewidth=0.75,
                       zorder=2)
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_aspect("equal", adjustable="box")
        tick_spacing = 0.2 if species == "H*" else 0.4
        ax.xaxis.set_major_locator(MultipleLocator(tick_spacing))
        ax.yaxis.set_major_locator(MultipleLocator(tick_spacing))
        ax.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))
        ax.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))
        ax.set_title(title, fontsize=fs(15), fontweight="bold", pad=7)
        ax.set_ylabel(rf"${ENERGY_SYMBOLS[species]}$ (eV)", fontsize=fs(12))
        m = metrics.loc[species]
        ax.text(0.03, 0.97,
                rf"$R^2$ = {m.FIT_R2:.2f}" + "\n" +
                f"RMSE = {m.FIT_RMSE_EV:.2f} eV\n" +
                f"LOMO = {m.LOMO_RMSE_EV:.2f} eV\n" +
                rf"$C_x$ = {m.COMPLEXITY:g}",
                transform=ax.transAxes, va="top", ha="left", fontsize=fs(10.5))
        ax.text(0.98, 0.04, f"({chr(ord('a')+j)})", transform=ax.transAxes,
                ha="right", va="bottom", fontsize=fs(13), fontweight="bold")
        ax.add_artist(AnnotationBbox(
            colored_equation(equation), (0.5, -0.19), xycoords=ax.transAxes,
            frameon=False, box_alignment=(0.5, 0.5), annotation_clip=False))

    for j, (scope, title) in enumerate([
        ("pooled", "All predicted species"),
        ("transition_states", "Transition states only"),
    ]):
        ax = axes[1, j]
        subset = matched.loc[matched.scope.eq(scope)].sort_values("set_number")
        assert len(subset) == 10
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_box_aspect(1)
        ax.set_xlim(0.12, 14)
        ax.set_ylim(0.12, 14)
        for axis in (ax.xaxis, ax.yaxis):
            axis.set_major_locator(FixedLocator([1, 10]))
            axis.set_major_formatter(FormatStrFormatter("%.1f"))
        ax.plot([0.12, 14], [0.12, 14], "--", lw=1, color="0.55", zorder=0)
        ax.scatter(subset.linear_fit_RMSE_eV, subset.nonlinear_fit_RMSE_eV,
                   marker="o", s=45, facecolor="white", edgecolor="black",
                   linewidth=1, label="Fit", zorder=2)
        ax.scatter(subset.linear_nested_RMSE_eV,
                   subset.nonlinear_nested_RMSE_eV,
                   marker="s", s=45, facecolor="0.45", edgecolor="black",
                   linewidth=0.8, label="LOMO", zorder=2)
        ax.set_title(title, fontsize=fs(12.5), fontweight="bold", pad=7)
        ax.set_xlabel("Linear RMSE (eV)", fontsize=fs(11))
        ax.set_ylabel("Nonlinear RMSE (eV)", fontsize=fs(11))
        ax.text(0.98, 0.04, f"({chr(ord('c')+j)})", transform=ax.transAxes,
                ha="right", va="bottom", fontsize=fs(13), fontweight="bold")
        if j == 0:
            ax.legend(frameon=False, fontsize=fs(8.5), loc="upper left")
        label_specs = {
            "pooled": {
                "CH*+HCO*+O*+H2O*": (r"$\mathrm{O}^{*}+\mathrm{H_2O}^{*}$", (9, -1)),
                "CH*+HCO*+O*+CO2*": (r"$\mathrm{O}^{*}+\mathrm{CO_2}^{*}$", (27, -19)),
                "CH*+HCO*+H*+H2O*": (r"$\mathrm{H}^{*}+\mathrm{H_2O}^{*}$", (-55, 34)),
            },
            "transition_states": {
                "CH*+HCO*+O*+H2O*": (r"$\mathrm{O}^{*}+\mathrm{H_2O}^{*}$", (10, -3)),
                "CH*+HCO*+O*+CO2*": (r"$\mathrm{O}^{*}+\mathrm{CO_2}^{*}$", (25, 16)),
                "CH*+HCO*+CO* FCC+H2O*": (r"$\mathrm{CO}_{\mathrm{FCC}}^{*}+\mathrm{H_2O}^{*}$", (-70, 40)),
            },
        }
        for point in subset.itertuples(index=False):
            if point.nonlinear_nested_RMSE_eV >= point.linear_nested_RMSE_eV:
                continue
            label, offset = label_specs[scope][point.set_key]
            ax.annotate(label,
                        (point.linear_nested_RMSE_eV,
                         point.nonlinear_nested_RMSE_eV),
                        xytext=offset, textcoords="offset points",
                        fontsize=fs(8.2), ha="left", va="center",
                        arrowprops={"arrowstyle": "-", "color": "0.4", "lw": 0.6})

    metal_handles = [Line2D([], [], marker=mark, linestyle="none",
                            markerfacecolor="white", markeredgecolor="0.2",
                            markersize=7, label=metal)
                     for metal, mark in MARKERS.items()]
    charge_handles = [Line2D([], [], marker="o", linestyle="none",
                             markerfacecolor=color[q], markeredgecolor="0.2",
                             markersize=7, label=label)
                      for q, label in zip(CHARGES,
                          (r"$-\frac{2}{18}$", r"$-\frac{1}{18}$", "0",
                           r"$+\frac{1}{18}$", r"$+\frac{2}{18}$"))]
    for y, label, handles, positions in (
        (0.978, "Metals", metal_handles, (0.19, 0.33, 0.47, 0.61, 0.75)),
        (0.940, "Charge", charge_handles,
         (0.19, 0.35, 0.51, 0.67, 0.83))):
        fig.text(0.035, y, label, ha="left", va="center", fontsize=fs(9.5))
        for x, handle in zip(positions, handles):
            fig.add_artist(Line2D([x], [y], linestyle="none",
                                  marker=handle.get_marker(), markersize=6,
                                  markerfacecolor=handle.get_markerfacecolor(),
                                  markeredgecolor="0.2",
                                  transform=fig.transFigure))
            fig.text(x + 0.018, y, handle.get_label(), ha="left", va="center",
                     fontsize=fs(9.5))
    for x, var, energy in (
        (0.10, "x", r"\Delta E_{\mathrm{CH}^{*}}^{\mathrm{ad}}"),
        (0.32, "y", r"\Delta E_{\mathrm{HCO}^{*}}^{\mathrm{ad}}"),
        (0.56, "z", r"\Delta E_{\mathrm{H_2O}^{*}}^{\mathrm{ad}}"),
        (0.80, "w", r"\Delta E_{\mathrm{O}^{*}}^{\mathrm{ad}}"),
    ):
        fig.text(x, 0.900, rf"$\mathbf{{{var}}}={energy}$",
                 color=VARIABLE_COLORS[var], ha="center", va="center", fontsize=fs(10.5))
    fig.savefig(HERE / "FIGURE_5.PNG", dpi=1200)
    plt.close(fig)


if __name__ == "__main__":
    main()
