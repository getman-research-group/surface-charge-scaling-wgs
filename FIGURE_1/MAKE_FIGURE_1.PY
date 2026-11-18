"""Recreate manuscript Figure 1 as a 1200-DPI geometry schematic.

The OH* and OC--OH TS values are read from manuscript Fig. 1; the CO* values
are from the Pt rows of the SI geometry table. Schematic atom positions are
illustrative; numeric distances and angles are from those sources.
"""

from pathlib import Path
from math import cos, radians, sin, tan

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse, FancyArrowPatch


OUT = Path(__file__).resolve().parent

plt.rcParams.update(
    {
        "font.family": "STIXGeneral",
        "mathtext.fontset": "stix",
        "font.size": 20,
        "axes.linewidth": 0.8,
        "savefig.facecolor": "white",
    }
)

CHARGES = [r"$-\frac{1}{18}$", r"$0$", r"$+\frac{1}{18}$"]
ROW_COLORS = ["#222222"] * 3
PLATINUM = "#d7d7e2"
OXYGEN = "#e52525"
CARBON = "#777777"
EDGE = "#555555"
FIGURE_HEIGHT = 7.2
Y_MIN = 0.20
AXES_ASPECT = ((11.0 * (0.995 - 0.085)) /
               (FIGURE_HEIGHT * (0.85 - 0.055)) * (1.0 - Y_MIN))

DATA = [
    {
        "name": r"OH$^*$",
        "atom": "O",
        "distance": ["2.01", "1.98", "1.95"],
        "angle": ["103.1", "106.4", "112.5"],
        "distance_label": r"$d_{\mathrm{Pt-O}}$",
        "angle_label": r"$\angle\mathrm{Pt-O-H}$",
    },
    {
        "name": r"CO$^*_{\mathrm{FCC}}$",
        "atom": "C",
        "distance": ["2.1", "2.1", "2.1"],
        "angle": ["128.4", "128.2", "128.5"],
        "distance_label": r"$d_{\mathrm{Pt-C}}$",
        "angle_label": r"$\angle\mathrm{Pt-C-O}$",
    },
    {
        "name": "OC–OH TS",
        "atom": "O",
        "distance": ["2.13", "2.1", "2.07"],
        "angle": ["102", "106", "113"],
        "distance_label": r"$d_{\mathrm{Pt-O}}$",
        "angle_label": r"$\angle\mathrm{Pt-O-H}$",
    },
]


def ball(ax, x, y, r, color, edge=EDGE, lw=0.8, z=3):
    # The panels are wider than tall, so correct x radii to keep atoms round.
    ax.add_patch(Ellipse((x, y), 2 * r / AXES_ASPECT, 2 * r,
                         facecolor=color, edgecolor=edge, lw=lw, zorder=z))


def angle_arc(ax, x, y, first_degree, span_degree, radius=0.100):
    a = np.deg2rad(np.linspace(first_degree, first_degree + span_degree, 90))
    ax.plot(x + radius * np.cos(a) / AXES_ASPECT,
            y + radius * np.sin(a), color="#666666", lw=0.85, zorder=4)


def surface(ax):
    # Slightly overlapping spheres match the manuscript's side-view shorthand.
    for x in [0.125, 0.275, 0.425, 0.575, 0.725, 0.875]:
        ball(ax, x, 0.40, 0.145, PLATINUM, edge="#8a8b9b", lw=0.75, z=2)
    ax.plot([0, 1], [0.54, 0.54], color="#a0a0a0", lw=0.65, ls=(0, (3, 3)), zorder=1)


def distance_arrow(ax, arrow_x, pt_x, atom_x, atom_bottom, text, color, label_side):
    pt_top = 0.54
    guide_color = "#777777"
    ax.plot([pt_x, arrow_x], [pt_top, pt_top], color=guide_color, lw=0.75, zorder=4)
    ax.plot([atom_x, arrow_x], [atom_bottom, atom_bottom], color=guide_color, lw=0.75, zorder=4)
    ax.add_patch(FancyArrowPatch((arrow_x, pt_top), (arrow_x, atom_bottom),
                                 arrowstyle="<->", mutation_scale=9,
                                 color="#3e3e3e", lw=0.8, zorder=5))
    offset = -0.015 if label_side == "left" else 0.015
    label, value = text.split(" = ", 1)
    for y, line in ((atom_bottom + 0.09, label), (pt_top + 0.075, value)):
        ax.text(arrow_x + offset, y, line,
                ha="right" if label_side == "left" else "left", va="center",
                fontsize=19, color=color, zorder=6)


def diagonal_bond_arrow(ax, pt_x, atom_x, atom_center, text, color):
    """Show the angled Pt--C bond and dimension it to the carbon edge."""
    pt_top = 0.54
    ax.plot([pt_x, atom_x], [pt_top, atom_center],
            color="#777777", lw=0.9, ls="--", zorder=2.5)
    right_pt_x = 0.575
    carbon_edge_x = atom_x
    carbon_edge_y = atom_center - 0.060
    ax.add_patch(FancyArrowPatch((right_pt_x, pt_top),
                                 (carbon_edge_x, carbon_edge_y),
                                 arrowstyle="<->", mutation_scale=9,
                                 shrinkA=0, shrinkB=0,
                                 color="#3e3e3e", lw=0.8, zorder=5))
    label, value = text.split(" = ", 1)
    for y, line in ((carbon_edge_y + 0.09, label), (pt_top + 0.075, value)):
        ax.text(0.61, y, line, ha="left", va="center",
                fontsize=19, color=color, zorder=6)


def geometry(ax, col, row):
    surface(ax)
    angle_text = r"$\theta=" + DATA[col]["angle"][row] + r"^\circ$"
    # Atom coordinates vary gently with charge. The drawing is schematic;
    # tabulated numbers, rather than these coordinates, define the geometry.
    if col == 0:  # OH* atop Pt
        ox = 0.50
        oy = 0.73 + [0.010, 0.005, -0.005][row]
        angle = float(DATA[col]["angle"][row])
        hx = ox + 0.13
        hy = oy + 0.13 * AXES_ASPECT * tan(radians(angle - 90))
        ax.plot([ox, hx], [oy, hy], color="#555555", lw=1.55, zorder=3)
        ball(ax, ox, oy, 0.064, OXYGEN, edge="#bd1e1e")
        ball(ax, hx, hy, 0.030, "white", edge="#777777")
        ax.plot([0.50, ox], [0.52, oy - 0.07], color="#777777", lw=0.9, ls="--", zorder=2)
        angle_arc(ax, ox, oy, 270, angle)
        distance_arrow(ax, 0.31, 0.425, ox, oy - 0.064,
                       DATA[col]["distance_label"] + " = " + DATA[col]["distance"][row] + r" $\mathrm{\AA}$",
                       ROW_COLORS[row], "left")
        ax.text(0.62, 0.655, angle_text, ha="left", va="center",
                fontsize=20, color=ROW_COLORS[row], zorder=6)
    elif col == 1:  # FCC CO*; C is the surface-binding atom
        cx = 0.50
        cy = 0.72 + [0, 0.003, 0.008][row]
        oy = cy + 0.145
        theta = float(DATA[col]["angle"][row])
        pt_to_c_direction = radians(90 + theta)
        pt_x = cx + ((0.54 - cy) * cos(pt_to_c_direction) /
                     sin(pt_to_c_direction) / AXES_ASPECT)
        ax.plot([cx, cx], [cy, oy], color="#555555", lw=1.8, zorder=3)
        ball(ax, cx, cy, 0.060, CARBON)
        ball(ax, cx, oy, 0.055, OXYGEN, edge="#bd1e1e")
        angle_arc(ax, cx, cy, 90, theta)
        diagonal_bond_arrow(ax, pt_x, cx, cy,
                            DATA[col]["distance_label"] + " = " + DATA[col]["distance"][row] + r" $\mathrm{\AA}$",
                            ROW_COLORS[row])
        ax.text(0.08, 0.795, angle_text, ha="left", va="center",
                fontsize=20, color=ROW_COLORS[row], zorder=6)
    else:  # OC--OH / CO--OH transition state
        ox = 0.35
        oy = 0.73 + [0.008, 0.004, 0][row]
        angle = float(DATA[col]["angle"][row])
        hx = ox + 0.10
        hy = oy + 0.10 * AXES_ASPECT * tan(radians(angle - 90))
        cx, cy = 0.68, 0.72
        co_oy = 0.865
        for start, end in [((ox, oy), (hx, hy)), ((cx, cy), (cx, co_oy))]:
            ax.plot([start[0], end[0]], [start[1], end[1]], color="#555555", lw=1.5, zorder=3)
        ax.plot([ox + 0.06, cx - 0.06], [oy + 0.01, cy + 0.015],
                color="#666666", lw=1.2, ls=(0, (2.5, 2)), zorder=3)
        ball(ax, ox, oy, 0.060, OXYGEN, edge="#bd1e1e")
        ball(ax, hx, hy, 0.028, "white", edge="#777777")
        ball(ax, cx, cy, 0.058, CARBON)
        ball(ax, cx, co_oy, 0.055, OXYGEN, edge="#bd1e1e")
        ax.plot([0.36, ox], [0.52, oy - 0.06], color="#777777", lw=0.9, ls="--", zorder=2)
        angle_arc(ax, ox, oy, 270, angle)
        distance_arrow(ax, 0.32, 0.275, ox, oy - 0.060,
                       DATA[col]["distance_label"] + " = " + DATA[col]["distance"][row] + r" $\mathrm{\AA}$",
                       ROW_COLORS[row], "left")
        ax.text(0.44, 0.61, angle_text, ha="left", va="center",
                fontsize=20, color=ROW_COLORS[row], zorder=6)


fig, axs = plt.subplots(3, 3, figsize=(11.0, FIGURE_HEIGHT), dpi=220)
fig.subplots_adjust(left=0.085, right=0.995, top=0.85, bottom=0.055,
                    wspace=0.025, hspace=0.035)

for col, spec in enumerate(DATA):
    top_panel = axs[0, col].get_position()
    column_center = (top_panel.x0 + top_panel.x1) / 2
    fig.text(column_center, 0.965, spec["name"], ha="center", va="center",
             fontsize=24, fontweight="bold")
    fig.text(column_center, 0.895,
             r"$\theta=" + spec["angle_label"][1:-1] + r"$",
             ha="center", va="center", fontsize=17, color="#555555")
    for row in range(3):
        ax = axs[row, col]
        ax.set_xlim(0, 1)
        ax.set_ylim(Y_MIN, 1)
        ax.set_aspect("auto")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_facecolor("white")
        for spine in ax.spines.values():
            spine.set_edgecolor("#444444")
            spine.set_linewidth(0.72)
        geometry(ax, col, row)
        ax.text(0.025, 0.96, f"({chr(97 + 3 * col + row)})", va="top", ha="left",
                fontsize=22, color="#111111", transform=ax.transAxes)

for row, charge in enumerate(CHARGES):
    axs[row, 0].text(-0.12, 0.51, charge, ha="center", va="center",
                     fontsize=22, fontweight="bold", color=ROW_COLORS[row],
                     rotation=90, transform=axs[row, 0].transAxes)

fig.text(0.016, 0.455, r"Surface charge ($h^+$ per surface Pt atom)",
         ha="center", va="center", rotation=90, fontsize=21, color="#222222")

fig.savefig(OUT / "FIGURE_1.PNG", dpi=1200,
            bbox_inches="tight", pad_inches=0.045)
plt.close(fig)
