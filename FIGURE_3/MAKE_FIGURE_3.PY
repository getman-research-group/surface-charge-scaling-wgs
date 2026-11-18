"""Recreate manuscript Figure 3 from repository equations and DFT energies.

Run from anywhere with: python FIGURE_3/MAKE_FIGURE_3.PY
The neutral equations are applied to charged surfaces without refitting.
The output is FIGURE_3.PNG at 1200 DPI.
"""

from pathlib import Path
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
from matplotlib.offsetbox import AnnotationBbox, HPacker, TextArea, VPacker
import numpy as np
import pandas as pd
import sympy as sp
from sympy.parsing.sympy_parser import (
    parse_expr, standard_transformations, implicit_multiplication_application,
)

OUT = Path(__file__).resolve().parent
REPO = OUT.parent
TARGETS = ("H*", "OC--OH")
MODELS = (
    ("Neutral linear", ("CH*", "HCO*"),
     REPO / "LINEAR_NEUTRAL" / "RESULTS_CH_HCO.PKL", "catmap"),
    ("Neutral nonlinear", ("OH*", "COOH*"),
     REPO / "NONLINEAR_NEUTRAL" / "RESULTS_OH_COOH.PKL", "symantic"),
)
MARKERS = {"Pt": "o", "Cu": "s", "Ru": "^", "Ni": "D", "Au": "v"}
CHARGES = (-0.11, -0.05, 0.0, 0.05, 0.11)
PLOT_CHARGES = (-0.11, -0.05, 0.05, 0.11)
COLORS = ("#748eba", "#b4c7e7", "#FFFFFF", "#f8cbad", "#d69163")
DESCRIPTOR_COLORS = ("#7041a6", "#157a52")
ROW_NAMES = {
    "Neutral linear": "Linear SRs\nTrain: neutral\nTest: charged data",
    "Neutral nonlinear": "Nonlinear SRs\nTrain: neutral\nTest: charged data",
}
DESCRIPTOR_LABELS = {
    "CH*": r"\mathbf{CH}^{*}",
    "HCO*": r"\mathbf{HCO}^{*}",
    "OH*": r"\mathbf{OH}^{*}",
    "COOH*": r"\mathbf{COOH}^{*}",
}
ENERGY_TERMS = {
    "CH*": r"\Delta E_{\mathrm{CH}^{*}}^{\mathrm{ad}}",
    "HCO*": r"\Delta E_{\mathrm{HCO}^{*}}^{\mathrm{ad}}",
    "OH*": r"\Delta E_{\mathrm{OH}^{*}}^{\mathrm{ad}}",
    "COOH*": r"\Delta E_{\mathrm{COOH}^{*}}^{\mathrm{ad}}",
}

def descriptor_label(name, variable):
    return (rf"$\Delta E_{{{DESCRIPTOR_LABELS[name]}}}^{{\mathrm{{ad}}}}"
            rf"\, (\mathbf{{{variable}}})$")


def equation_row(path, pair, target):
    saved = pd.read_pickle(path)
    selection = saved.loc[saved.TARGET.eq(target)]
    if len(selection) != 1:
        raise ValueError(f"Expected one equation for {target} in {path}")
    return selection.iloc[0]


def equation_as_xy(equation, pair):
    expression = equation
    for descriptor, name in zip(pair, ("X0", "X1")):
        term = ENERGY_TERMS[descriptor]
        if term not in expression:
            raise ValueError(f"{descriptor} does not occur in {equation}")
        expression = expression.replace(term, name)
    expression = re.sub(
        r"\\frac\{(X[01])\}\{(X[01])\}",
        lambda match: f"({match.group(1)}/{match.group(2)})", expression,
    )
    expression = expression.replace(r"\left(", "(").replace(r"\right)", ")")
    if "\\" in expression:
        raise ValueError(f"Unsupported saved equation: {expression}")
    return expression


def evaluate(row, pair, kind, frame):
    X0, X1 = sp.symbols("X0 X1")
    expression = parse_expr(
        equation_as_xy(row.EQUATION, pair),
        local_dict={"X0": X0, "X1": X1},
        transformations=standard_transformations + (implicit_multiplication_application,),
    )
    function = sp.lambdify((X0, X1), expression, modules="numpy")
    values = function(frame[pair[0]].to_numpy(float), frame[pair[1]].to_numpy(float))
    values = np.asarray(values, dtype=float)
    if values.shape != (len(frame),) or not np.isfinite(values).all():
        raise ValueError("Equation did not yield finite predictions")
    return values


def rmse(frame):
    return float(np.sqrt(np.mean((frame.predicted_eV-frame.actual_eV)**2)))


def r2(frame):
    y=frame.actual_eV.to_numpy(float)
    p=frame.predicted_eV.to_numpy(float)
    return float(1-np.sum((p-y)**2)/np.sum((y-y.mean())**2))


def split_terms(expression):
    terms=[]
    start=0
    depth=0
    for index,char in enumerate(expression):
        if char=="(":
            depth+=1
        elif char==")":
            depth-=1
        elif char in "+-" and depth==0 and index>start:
            terms.append(expression[start:index].strip())
            start=index
    terms.append(expression[start:].strip())
    if depth:
        raise ValueError("Unbalanced equation parentheses")
    return terms


def energy_symbol(name):
    if name == "OC--OH":
        return r"\Delta E_{\mathrm{OC{-}OH}}^{\mathrm{act}}"
    species = {"CH*": r"\mathrm{CH}^{*}", "HCO*": r"\mathrm{HCO}^{*}",
               "CO2*": r"\mathrm{CO}_{2}^{*}", "COOH*": r"\mathrm{COOH}^{*}",
               "H*": r"\mathrm{H}^{*}", "CO* FCC": r"\mathrm{CO}_{\rm FCC}^{*}",
               "H2O*": r"\mathrm{H}_{2}\mathrm{O}^{*}",
               "O*": r"\mathrm{O}^{*}", "OH*": r"\mathrm{OH}^{*}"}
    return rf"\Delta E_{{{species[name]}}}^{{\mathrm{{ad}}}}"


def equation_label(model, pair, target, selected, kind):
    """Round saved coefficients to two decimals and color the descriptor symbols."""
    expression = equation_as_xy(selected.EQUATION, pair)
    names = {
        "X0": r"\textcolor[HTML]{7041A6}{\boldsymbol{x}}",
        "X1": r"\textcolor[HTML]{157A52}{\boldsymbol{y}}",
    }
    rendered=[]
    for term in split_terms(expression):
        sign="-" if term.startswith("-") else "+"
        unsigned=term[1:].strip() if term.startswith(("+","-")) else term
        match=re.fullmatch(r"(\d+(?:\.\d+)?)(.*)",unsigned)
        if not match:
            raise ValueError(f"Cannot format equation term {term}")
        coefficient,basis=match.groups()
        basis=re.sub(r"\s+","",basis)
        if basis.startswith("(") and basis.endswith(")") and "/" in basis:
            basis=basis[1:-1]
        basis=basis.replace("X0",names["X0"]).replace("X1",names["X1"])
        rendered.append((sign,f"{float(coefficient):.2f}"+basis))
    label_lines=[]
    for start in range(0,len(rendered),2):
        chunk=rendered[start:start+2]
        content="".join((sign if start or n else ("-" if sign=="-" else ""))+body
                        for n,(sign,body) in enumerate(chunk))
        label_lines.append("$" + content + "$")
    return "\n".join(label_lines)


def fit_annotation(charged_r2, charged_rmse):
    return (rf"$R^2$ = {charged_r2:.2f}" + "\n"
            + f"RMSE = {charged_rmse:.3f} eV")


def equation_rhs_box(label):
    """Build a colored RHS using the same DejaVu math renderer as the axes."""
    x_token = r"\textcolor[HTML]{7041A6}{\boldsymbol{x}}"
    y_token = r"\textcolor[HTML]{157A52}{\boldsymbol{y}}"
    splitter = re.compile("(" + re.escape(x_token) + "|" + re.escape(y_token) + ")")
    lines = []
    for raw_line in label.splitlines():
        expression = raw_line[1:-1]
        expression = expression.replace(r"\frac{" + y_token + "}{" + x_token + "}",
                                        y_token + "/" + x_token)
        expression = expression.replace(r"\left(", "(").replace(r"\right)", ")")
        fragments = []
        for token in splitter.split(expression):
            if not token.strip():
                continue
            if token in (x_token, y_token):
                symbol = "x" if token == x_token else "y"
                color = DESCRIPTOR_COLORS[0 if symbol == "x" else 1]
                tex = rf"$\mathbf{{{symbol}}}$"
            else:
                color = "#202020"
                tex = "$" + token + "$"
            fragments.append(TextArea(tex, textprops={
                "fontsize": 16.5, "color": color,
                "fontfamily": "DejaVu Sans", "math_fontfamily": "dejavusans"}))
        lines.append(HPacker(children=fragments, align="baseline", pad=0, sep=0.8))
    return VPacker(children=lines, align="center", pad=0, sep=3.0)


def main():
    data = pd.read_pickle(REPO / "ADSORPTION_TS_ENERGIES.PKL")
    if len(data) != 25:
        raise ValueError("Expected 25 metal-charge observations")
    rows=[]
    for model,pair,path,kind in MODELS:
        for target in TARGETS:
            selected=equation_row(path,pair,target)
            available=data.loc[data[target].notna() & data[list(pair)].notna().all(axis=1)]
            prediction=evaluate(selected,pair,kind,available)
            neutral=np.isclose(available.Charge,0)
            fit_rmse=np.sqrt(np.mean(
                (prediction[neutral]-available.loc[neutral,target].to_numpy(float))**2
            ))
            if abs(fit_rmse-selected.FIT_RMSE_EV)>5e-5:
                raise ValueError(f"Saved equation does not match {target} fit RMSE")
            for (_,observation),value in zip(available.iterrows(),prediction):
                if np.isclose(observation.Charge,0):
                    continue
                rows.append({
                    "model":model,"target":target,"metal":observation.Metal,
                    "charge":float(observation.Charge),
                    "actual_eV":float(observation[target]),
                    "predicted_eV":float(value),
                })
    cases=pd.DataFrame(rows)
    if len(cases)!=72:
        raise ValueError(f"Expected 72 charged observations, got {len(cases)}")
    summary=[]
    for (model,target),subset in cases.groupby(["model","target"],sort=False):
        summary.append({
            "model":model,"target":target,
            "charged_R2":r2(subset),"charged_RMSE_eV":rmse(subset),
        })
    summary=pd.DataFrame(summary)

    plt.rcParams.update({"font.size": 17, "text.usetex":False,
                         "axes.spines.top": True, "axes.spines.right": True})
    figure_height = 10.2
    fig, axes = plt.subplots(2, 2, figsize=(8.6, figure_height), constrained_layout=False)
    color_for = dict(zip(CHARGES, COLORS))
    equation_annotations = []
    for i, (model, pair, path, kind) in enumerate(MODELS):
        for j, target in enumerate(TARGETS):
            ax = axes[i,j]
            sub = cases.loc[(cases.model == model)&(cases.target == target)]
            target_data = cases.loc[cases.target == target]
            extrema = np.r_[target_data.actual_eV.to_numpy(), target_data.predicted_eV.to_numpy()]
            lo = np.floor((extrema.min()-0.13)*2)/2
            hi = np.ceil((extrema.max()+0.13)*2)/2
            ax.plot([lo,hi], [lo,hi], color="#555555", lw=1.1, ls="--", zorder=1)
            for _, point in sub.iterrows():
                metal = str(point.metal).split("(")[0]
                charge = min(CHARGES, key=lambda c: abs(c-point.charge))
                assert abs(charge-point.charge) < 1e-8
                ax.scatter(point.predicted_eV, point.actual_eV,
                           marker=MARKERS[metal], s=47,
                           facecolor=color_for[charge], edgecolor="#303030",
                           linewidth=0.8, zorder=2)
            ax.set_xlim(lo,hi); ax.set_ylim(lo,hi)
            if target == "OC--OH":
                ax.xaxis.set_major_locator(MaxNLocator(nbins=4, integer=True))
                ax.yaxis.set_major_locator(MaxNLocator(nbins=4, integer=True))
            ax.set_aspect("equal", adjustable="box")
            for spine in ax.spines.values():
                spine.set_visible(True)
                spine.set_linewidth(0.8)
            ax.set_ylabel(rf"${energy_symbol(target)}$ (eV)", fontsize=17,
                          labelpad=2 if j == 1 else -2)
            ax.yaxis.set_label_coords(-0.18 if j == 1 else -0.30,0.5)
            selected = equation_row(path,pair,target)
            equation_box = equation_rhs_box(
                equation_label(model,pair,target,selected,kind))
            annotation = AnnotationBbox(
                equation_box, (0.5,-0.25),
                xycoords=ax.transAxes, frameon=False,
                box_alignment=(0.5,0.5), annotation_clip=False)
            ax.add_artist(annotation)
            equation_annotations.append(annotation)
            panel=summary.loc[(summary.model==model)&(summary.target==target)].iloc[0]
            metrics_label=fit_annotation(panel.charged_R2,panel.charged_RMSE_eV)
            ax.text(0.03, 0.97, metrics_label,
                    transform=ax.transAxes, ha="left", va="top", fontsize=15.0,
                    path_effects=[path_effects.withStroke(linewidth=2.5,
                                                         foreground="white")])
            ax.text(0.97,0.04,f"({chr(ord('a')+2*i+j)})",
                    transform=ax.transAxes,ha="right",va="bottom",
                    fontsize=17.5,fontweight="bold")
    metal_handles = [Line2D([], [], linestyle="none", marker=m, markersize=7,
                            markerfacecolor="white", markeredgecolor="#303030", label=metal)
                     for metal,m in MARKERS.items()]
    charge_labels = (r"$-\frac{2}{18}$", r"$-\frac{1}{18}$",
                     r"$+\frac{1}{18}$", r"$+\frac{2}{18}$")
    charge_handles = [Line2D([], [], linestyle="none", marker="o", markersize=7,
                             markerfacecolor=color_for[c], markeredgecolor="#303030",
                             label=label)
                      for c, label in zip(PLOT_CHARGES, charge_labels)]
    # Use the blank left margin for legend labels; keep each marker row compact.
    legend_x = (0.225, 0.365, 0.505, 0.645, 0.785)
    charge_legend_x = (0.225, 0.405, 0.585, 0.765)
    legend_rows = ((0.060, "Metals", metal_handles),
                   (0.038, "Charge", charge_handles))
    legend_artists = {}
    for legend_y, legend_label, handles in legend_rows:
        heading = ("Charge" if legend_label == "Charge"
                   else legend_label)
        artists = [fig.text(0.055, legend_y, heading,
                            ha="left", va="center",
                            fontsize=12.0 if legend_label == "Charge" else 16.5)]
        marker_x = charge_legend_x if legend_label == "Charge" else legend_x
        for x, handle in zip(marker_x, handles):
            marker = Line2D([x], [legend_y], linestyle="none",
                            marker=handle.get_marker(), markersize=7.5,
                            markerfacecolor=handle.get_markerfacecolor(),
                            markeredgecolor="#303030",
                            transform=fig.transFigure)
            fig.add_artist(marker)
            artists.append(marker)
            artists.append(fig.text(x + 0.024, legend_y, handle.get_label(),
                                    ha="left", va="center", fontsize=16.5))
        legend_artists[legend_label] = artists
    plot_width=0.32
    plot_height=plot_width*8.6/figure_height
    for i,center_y in enumerate(value / figure_height for value in (7.9,3.9)):
        for j,left_x in enumerate((0.21,0.625)):
            axes[i,j].set_position([left_x,center_y-plot_height/2,
                                    plot_width,plot_height])
    for j,target in enumerate(TARGETS):
        position=axes[0,j].get_position()
        title=target if target == "H*" else "OC--OH TS"
        fig.text((position.x0+position.x1)/2,(7.9 + plot_width*8.6/2 + 0.45)/figure_height,title,
                 ha="center",va="center",fontsize=23,fontweight="bold")
    pair_texts=[]
    for i,(model,pair,_,_) in enumerate(MODELS):
        position=axes[i,0].get_position()
        middle=(position.y0+position.y1)/2
        fig.text(-0.005,middle,ROW_NAMES[model],ha="center",va="center",
                 rotation=90,fontsize=13.8,fontweight="bold",
                 linespacing=1.05)
        parts=[]
        # With 90-degree rotation, the pair reads from bottom to top.
        for label,color in (("(","#222222"),
                            (descriptor_label(pair[0], "x"),DESCRIPTOR_COLORS[0]),
                            (",","#222222"),
                            (descriptor_label(pair[1], "y"),DESCRIPTOR_COLORS[1]),
                            (")","#222222")):
            parts.append(fig.text(0.058,middle,label,ha="center",va="center",
                                  rotation=90,fontsize=16.4,fontweight="bold",
                                  fontfamily="DejaVu Sans",usetex=False,color=color))
        pair_texts.append((middle,parts))
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    for middle,parts in pair_texts:
        heights=[part.get_window_extent(renderer).height/fig.bbox.height
                 for part in parts]
        gap=0.0035
        bottom=middle-(sum(heights)+(len(parts)-1)*gap)/2
        for part,height in zip(parts,heights):
            part.set_y(bottom+height/2)
            bottom+=height+gap
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    lowest_equation = min(annotation.get_window_extent(renderer).y0
                          for annotation in equation_annotations[-2:])
    metal_text_height = max(artist.get_window_extent(renderer).height
                            for artist in legend_artists["Metals"]
                            if not isinstance(artist, Line2D))
    metal_y = (lowest_equation - 0.12*fig.dpi - metal_text_height/2) / fig.bbox.height
    charge_y = metal_y - 0.40/figure_height
    for label, new_y in (("Metals", metal_y), ("Charge", charge_y)):
        for artist in legend_artists[label]:
            if isinstance(artist, Line2D):
                artist.set_ydata([new_y])
            else:
                artist.set_y(new_y)
    fig.savefig(OUT / "FIGURE_3.PNG", dpi=1200, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
