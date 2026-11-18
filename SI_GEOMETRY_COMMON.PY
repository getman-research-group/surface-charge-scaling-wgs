"""Recreate SI geometry schematics from the measured values in this repository.

The labels give tabulated DFT distances and angles. Drawn atom positions are
illustrative and are not DFT atomic coordinates.
"""

from __future__ import annotations

from dataclasses import dataclass
import csv
from math import cos, radians, sin
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Arc, Ellipse, FancyArrowPatch


HERE = Path(__file__).resolve().parent
MEASUREMENTS_CSV = HERE / "GEOMETRY_MEASUREMENTS.CSV"

plt.rcParams.update({
    "font.family": "STIXGeneral",
    "mathtext.fontset": "stix",
    "font.size": 11,
    "axes.linewidth": 0.8,
    "savefig.facecolor": "white",
    "text.latex.preamble": r"\usepackage{amsmath}",
})

CHARGES = [r"$-\tfrac{1}{18}$", r"$0$", r"$+\tfrac{1}{18}$"]
CHARGE_COLORS = ["#222222"] * 3
METAL_COLORS = {
    "Pt": ("#d7d7e2", "#8a8b9b"),
    "Cu": ("#cc8a3e", "#916025"),
    "Ru": ("#2d999a", "#1f6f70"),
    "Ni": ("#58cf54", "#328f33"),
    "Au": ("#ffd42d", "#bc9b13"),
}
ATOM_COLORS = {"O": "#e52525", "C": "#777777", "H": "white"}
ATOM_EDGES = {"O": "#b81d1d", "C": "#555555", "H": "#777777"}


@dataclass(frozen=True)
class Spec:
    number: int
    slug: str
    title: str
    key: str
    metals: tuple[str, ...]
    anchor: str
    angle_label: str | None = None


SPECS = [
    Spec(1, "oh", r"OH$^*$", "OH", ("Cu", "Ru", "Ni", "Au"), "O", r"$\angle M\!\!\!-O-H$"),
    Spec(2, "hco", r"HCO$^*$", "HCO", ("Pt", "Cu", "Ru", "Ni", "Au"), "C", r"$\angle M\!\!\!-C-O$"),
    Spec(3, "co_fcc", r"CO$^*_{\mathrm{FCC}}$", "CO", ("Pt", "Cu", "Ru", "Ni", "Au"), "C", r"$\angle M\!\!\!-C-O$"),
    Spec(4, "h2o", r"H$_2$O$^*$", "H2O", ("Pt", "Cu", "Ru", "Ni", "Au"), "O", r"$\angle M\!\!\!-O-H$"),
    Spec(5, "co2", r"CO$_2^*$", "CO2", ("Pt", "Cu", "Ru", "Ni", "Au"), "C", r"$\angle M\!\!\!-C-O$"),
    Spec(6, "h", r"H$^*$", "H", ("Pt", "Cu", "Ru", "Ni", "Au"), "H"),
    Spec(7, "o", r"O$^*$", "O", ("Pt", "Cu", "Ni", "Au"), "O"),
    Spec(8, "ch", r"CH$^*$", "CH", ("Pt", "Cu", "Ru", "Ni", "Au"), "C", r"$\angle M\!\!\!-C-H$"),
    Spec(9, "cooh", r"COOH$^*$", "COOH", ("Pt", "Cu", "Ru", "Ni", "Au"), "C", r"$\angle M\!\!\!-C-O$"),
    Spec(10, "h_co2_ts", "H--CO$_2$ TS", "H_CO2", ("Pt", "Cu", "Ni", "Au"), "C", r"$\angle M\!\!\!-C-O$"),
    Spec(11, "o_co_ts", "O--CO TS", "O_CO", ("Pt", "Cu", "Ni", "Au"), "C", r"$\angle M\!\!\!-C-O$"),
    Spec(12, "ho_h_ts", "HO--H TS", "HO_H", ("Cu", "Ni", "Au"), "O", r"$\angle M\!\!\!-O-H$"),
    Spec(13, "oc_oh_ts", "OC--OH TS", "OC_OH", ("Cu", "Ni", "Au"), "O", r"$\angle M\!\!\!-O-H$"),
]


def read_measurements() -> dict[tuple[str, str], tuple[list[str], list[str | None]]]:
    measurements = {}
    with MEASUREMENTS_CSV.open(newline="") as handle:
        for row in csv.DictReader(handle):
            distances = [row[f"DISTANCE_{state}_A"] for state in ("NEG", "ZERO", "POS")]
            angles = [row[f"ANGLE_{state}_DEG"] or None for state in ("NEG", "ZERO", "POS")]
            measurements[(row["METAL"], row["SPECIES"])] = (distances, angles)
    return measurements


def physical_ratio(ax) -> float:
    # A data-space x displacement must be divided by this to draw a circle.
    return ax.bbox.width / ax.bbox.height


def ball(ax, x: float, y: float, radius: float, element: str, ratio: float,
         *, z: int = 5):
    ax.add_patch(Ellipse((x, y), 2 * radius / ratio, 2 * radius,
                         facecolor=ATOM_COLORS[element], edgecolor=ATOM_EDGES[element],
                         linewidth=0.82, zorder=z))


def bond(ax, a, b, *, dashed=False, width=1.5, z=3):
    ax.plot((a[0], b[0]), (a[1], b[1]),
            color="#777777" if dashed else "#555555", linewidth=width,
            linestyle=(0, (3, 2.5)) if dashed else "-", zorder=z)


def polar(origin, length, degrees_ccw, ratio):
    angle = radians(degrees_ccw)
    return origin[0] + length * cos(angle) / ratio, origin[1] + length * sin(angle)


def arc(ax, center, begin, end, ratio, radius=0.088):
    while end < begin:
        end += 360
    if end - begin > 180:
        begin, end = end, begin + 360
    ax.add_patch(Arc(center, 2 * radius / ratio, 2 * radius,
                     theta1=begin, theta2=end, color="#666666",
                     linewidth=0.85, zorder=7))


def surface(ax, metal: str, ratio: float):
    fill, edge = METAL_COLORS[metal]
    for x in (0.12, 0.27, 0.42, 0.57, 0.72, 0.87):
        ax.add_patch(Ellipse((x, 0.34), 0.242 / ratio, 0.242,
                             facecolor=fill, edgecolor=edge, linewidth=0.75,
                             zorder=2))
    ax.plot((0, 1), (0.461, 0.461), color="#a0a0a0", linewidth=0.65,
            linestyle=(0, (3, 3)), zorder=1)


def main_atom_y(distance: float) -> float:
    # Height is schematic but increases with the tabulated metal--atom distance.
    return min(0.765, 0.565 + 0.071 * (distance - 1.5))


def distance_annotation(ax, x_atom: float, y_atom: float, element: str,
                        distance: str, metal: str, row: int):
    atom_radius = 0.032 if element == "H" else 0.053
    bottom = y_atom - atom_radius
    arrow_x = 0.285
    ax.plot((0.42, arrow_x), (0.461, 0.461), color="#777777", lw=0.66, zorder=7)
    ax.plot((x_atom, arrow_x), (bottom, bottom), color="#777777", lw=0.66, zorder=7)
    ax.add_patch(FancyArrowPatch((arrow_x, 0.461), (arrow_x, bottom),
                                 arrowstyle="<->", mutation_scale=6.4,
                                 shrinkA=0, shrinkB=0, color="#3e3e3e",
                                 linewidth=0.8, zorder=8))
    label = rf"$d_{{\mathrm{{{metal}-{element}}}}}={distance}\,\mathrm{{\AA}}$"
    ax.text(0.268, (0.461 + bottom) / 2, label, ha="right", va="center",
            fontsize=8.8, color=CHARGE_COLORS[row], zorder=10)


def angle_annotation(ax, value: str | None, row: int):
    if value is None:
        return
    ax.text(0.98, 0.505, rf"$\theta={value}^\circ$", ha="right", va="center",
            fontsize=10.0, color=CHARGE_COLORS[row], zorder=10)


def draw_molecule(ax, spec: Spec, center, theta: str | None, ratio: float):
    """Draw one consistent side-view motif, with a measured angle where known."""
    x, y = center
    t = float(theta) if theta is not None else None
    key = spec.key
    if key in ("H", "O"):
        ball(ax, x, y, 0.033 if key == "H" else 0.056, key, ratio)
        bond(ax, (x, 0.461), (x, y - (0.033 if key == "H" else 0.056)), dashed=True, width=0.85, z=3)
        return

    if key == "CO":
        # Keep the FCC CO construction used in the approved manuscript draft.
        oxygen = (x, y + 0.131)
        bond(ax, center, oxygen, width=1.55)
        ball(ax, x, y, 0.053, "C", ratio)
        ball(ax, *oxygen, 0.052, "O", ratio)
        left_pt = polar(center, 0.22, 90 + t, ratio)
        bond(ax, (left_pt[0], 0.461), center, dashed=True, width=0.9, z=3)
        arc(ax, center, 90, 90 + t, ratio)
        return

    if key == "OC_OH":
        # The O--H angle and two-fragment arrangement follow main Fig. 1.
        ox = x - 0.11
        oh_center = (ox, y)
        h = polar(oh_center, 0.14, 270 + t, ratio)
        carbon = (x + 0.17, y + 0.005)
        oxygen = (carbon[0], carbon[1] + 0.13)
        bond(ax, oh_center, h)
        bond(ax, carbon, oxygen)
        bond(ax, (oh_center[0] + 0.052 / ratio, y + 0.012),
             (carbon[0] - 0.05 / ratio, carbon[1] + 0.012), dashed=True, width=1.0)
        bond(ax, (ox, 0.461), oh_center, dashed=True, width=0.85)
        ball(ax, *oh_center, 0.053, "O", ratio)
        ball(ax, *h, 0.027, "H", ratio)
        ball(ax, *carbon, 0.052, "C", ratio)
        ball(ax, *oxygen, 0.050, "O", ratio)
        arc(ax, oh_center, 270, 270 + t, ratio)
        return

    anchor = spec.anchor
    ball(ax, x, y, 0.053, anchor, ratio)
    bond(ax, (x, 0.461), center, dashed=True, width=0.85)

    if key in ("OH", "CH"):
        h = polar(center, 0.145, 270 + t, ratio)
        bond(ax, center, h)
        ball(ax, *h, 0.028, "H", ratio)
        arc(ax, center, 270, 270 + t, ratio)
    elif key == "HCO":
        o = polar(center, 0.15, 270 + t, ratio)
        h = polar(center, 0.13, 135, ratio)
        bond(ax, center, o)
        bond(ax, center, h)
        ball(ax, *o, 0.052, "O", ratio)
        ball(ax, *h, 0.027, "H", ratio)
        arc(ax, center, 270, 270 + t, ratio)
    elif key == "H2O":
        h1 = polar(center, 0.127, 270 + t, ratio)
        h2 = polar(center, 0.125, 150, ratio)
        bond(ax, center, h1)
        bond(ax, center, h2)
        ball(ax, *h1, 0.027, "H", ratio)
        ball(ax, *h2, 0.027, "H", ratio)
        arc(ax, center, 270, 270 + t, ratio)
    elif key in ("CO2", "H_CO2"):
        o1 = polar(center, 0.15, 270 + t, ratio)
        o2 = polar(center, 0.15, 90 + t, ratio)
        bond(ax, center, o1)
        bond(ax, center, o2)
        ball(ax, *o1, 0.051, "O", ratio)
        ball(ax, *o2, 0.051, "O", ratio)
        if key == "H_CO2":
            h = (o1[0] + 0.052 / ratio, o1[1] - 0.046)
            bond(ax, o1, h, dashed=True, width=1.0)
            ball(ax, *h, 0.027, "H", ratio)
        arc(ax, center, 270, 270 + t, ratio)
    elif key == "COOH":
        oxygen = polar(center, 0.15, 270 + t, ratio)
        hydroxyl_o = polar(center, 0.14, 145, ratio)
        h = polar(hydroxyl_o, 0.105, 145, ratio)
        for p in (oxygen, hydroxyl_o):
            bond(ax, center, p)
            ball(ax, *p, 0.050, "O", ratio)
        bond(ax, hydroxyl_o, h)
        ball(ax, *h, 0.026, "H", ratio)
        arc(ax, center, 270, 270 + t, ratio)
    elif key == "O_CO":
        oxygen = polar(center, 0.15, 270 + t, ratio)
        incoming = polar(center, 0.18, 145, ratio)
        bond(ax, center, oxygen)
        bond(ax, center, incoming, dashed=True, width=1.0)
        ball(ax, *oxygen, 0.051, "O", ratio)
        ball(ax, *incoming, 0.051, "O", ratio)
        arc(ax, center, 270, 270 + t, ratio)
    elif key == "HO_H":
        attached_h = polar(center, 0.13, 270 + t, ratio)
        incoming_h = polar(center, 0.18, 155, ratio)
        bond(ax, center, attached_h)
        bond(ax, center, incoming_h, dashed=True, width=1.0)
        ball(ax, *attached_h, 0.027, "H", ratio)
        ball(ax, *incoming_h, 0.027, "H", ratio)
        arc(ax, center, 270, 270 + t, ratio)
    else:
        raise ValueError(key)


def draw_one(spec: Spec, measurements, output_path: Path):
    cols = len(spec.metals)
    fig_width = 1.0 + cols * 2.86
    fig, axs = plt.subplots(3, cols, figsize=(fig_width, 6.25), dpi=220,
                            squeeze=False)
    fig.subplots_adjust(left=0.066, right=0.997, top=0.855, bottom=0.045,
                        wspace=0.027, hspace=0.045)
    fig.canvas.draw()
    for col, metal in enumerate(spec.metals):
        if (metal, spec.key) not in measurements:
            raise ValueError(f"No SI geometry table row for {metal}, {spec.key}")
        ds, angles = measurements[(metal, spec.key)]
        axs[0, col].set_title(metal + ("(0001)" if metal == "Ru" else "(111)"),
                              fontsize=16, fontweight="bold", pad=8)
        for row in range(3):
            ax = axs[row, col]
            ax.set_xlim(0, 1)
            ax.set_ylim(0.09, 1)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_facecolor("white")
            for spine in ax.spines.values():
                spine.set_edgecolor("#444444")
                spine.set_linewidth(0.72)
            ratio = physical_ratio(ax) * (1 - 0.09)
            surface(ax, metal, ratio)
            anchor_x = 0.55 if spec.key != "OC_OH" else 0.66
            anchor_y = main_atom_y(float(ds[row]))
            draw_molecule(ax, spec, (anchor_x, anchor_y), angles[row], ratio)
            distance_annotation(ax, anchor_x if spec.key != "OC_OH" else anchor_x - 0.11,
                                anchor_y, spec.anchor, ds[row], metal, row)
            angle_annotation(ax, angles[row], row)
            ax.text(0.025, 0.975, f"({chr(97 + col * 3 + row)})",
                    transform=ax.transAxes, va="top", ha="left",
                    fontsize=13.5, color="#111111", zorder=10)
            if row == 0 and spec.angle_label:
                ax.text(0.98, 0.97, spec.angle_label,
                        transform=ax.transAxes, va="top", ha="right",
                        fontsize=9.7, color="#555555", zorder=10)
    for row, charge in enumerate(CHARGES):
        axs[row, 0].text(-0.12, 0.51, charge, transform=axs[row, 0].transAxes,
                         ha="center", va="center", rotation=90, fontsize=15,
                         fontweight="bold", color=CHARGE_COLORS[row], usetex=True)
    fig.suptitle(spec.title, fontsize=19, fontweight="bold", y=0.989)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=1200, bbox_inches="tight", pad_inches=0.045)
    plt.close(fig)
    return output_path


if __name__ == "__main__":
    data = read_measurements()
    generated = [draw_one(spec, data, HERE / f"FIGURE_S{spec.number:02d}" / f"FIGURE_S{spec.number:02d}.PNG") for spec in SPECS]
    print(f"Generated {len(generated)} SI geometry PNGs at 1200 DPI")
