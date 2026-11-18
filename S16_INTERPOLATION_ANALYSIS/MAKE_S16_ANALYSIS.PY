"""Calculate the Pt descriptor interpolation results used in SI Figure S16.

Run from the repository root:
    python S16_INTERPOLATION_ANALYSIS/MAKE_S16_ANALYSIS.PY

Writes the nine SI-table metrics and the node/curve values for Figure S16.
The fourth-degree polynomial through five nodes is the unique fourth-order
Lagrange interpolant; the quadratic is fitted by least squares.
"""

from pathlib import Path

import numpy as np
from numpy.polynomial import Polynomial
import pandas as pd


HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "ADSORPTION_TS_ENERGIES.PKL"
METRICS = HERE / "INTERPOLATION_METRICS.PKL"
NODES = HERE / "INTERPOLATION_NODES.PKL"
CURVES = HERE / "INTERPOLATION_CURVES.PKL"
DESCRIPTORS = (
    "H*", "O*", "OH*", "CH*", "CO* FCC",
    "CO2*", "H2O*", "HCO*", "COOH*",
)
FINE_CHARGE = np.linspace(-0.11, 0.11, 501)


def main():
    energies = pd.read_pickle(DATA)
    pt = energies.loc[energies["Metal"].eq("Pt")].sort_values("Charge")
    charges = pt["Charge"].to_numpy(float)
    if len(charges) != 5 or not np.allclose(
        charges, [-0.11, -0.05, 0.0, 0.05, 0.11]
    ):
        raise ValueError("Expected five Pt adsorption energies at the stored charges")
    scaled_nodes = charges / 0.11
    scaled_curve = FINE_CHARGE / 0.11
    metric_rows = []
    node_rows = []
    curve_rows = []

    for descriptor in DESCRIPTORS:
        values = pt[descriptor].to_numpy(float)
        if not np.isfinite(values).all():
            raise ValueError(f"Missing Pt energy for {descriptor}")
        quadratic = Polynomial.fit(scaled_nodes, values, deg=2)
        lagrange = Polynomial.fit(scaled_nodes, values, deg=4)
        quadratic_at_nodes = quadratic(scaled_nodes)
        lagrange_at_nodes = lagrange(scaled_nodes)
        residuals = quadratic_at_nodes - values
        lagrange_error = float(np.max(np.abs(lagrange_at_nodes - values)))
        if lagrange_error > 1e-10:
            raise ValueError(f"Interpolant misses a DFT node for {descriptor}")

        metric_rows.append({
            "DESCRIPTOR": descriptor,
            "QUADRATIC_FIT_RMSE_EV": float(np.sqrt(np.mean(residuals ** 2))),
            "QUADRATIC_MAX_ABS_ERROR_EV": float(np.max(np.abs(residuals))),
            "LAGRANGE_MAX_ABS_NODE_ERROR_EV": lagrange_error,
        })
        for charge, dft, quadratic_value, lagrange_value in zip(
            charges, values, quadratic_at_nodes, lagrange_at_nodes
        ):
            node_rows.append({
                "DESCRIPTOR": descriptor,
                "CHARGE": float(charge),
                "DFT_EV": float(dft),
                "QUADRATIC_EV": float(quadratic_value),
                "LAGRANGE_EV": float(lagrange_value),
            })
        for charge, quadratic_value, lagrange_value in zip(
            FINE_CHARGE, quadratic(scaled_curve), lagrange(scaled_curve)
        ):
            curve_rows.append({
                "DESCRIPTOR": descriptor,
                "CHARGE": float(charge),
                "QUADRATIC_EV": float(quadratic_value),
                "LAGRANGE_EV": float(lagrange_value),
            })

    metrics = pd.DataFrame(metric_rows)
    nodes = pd.DataFrame(node_rows)
    curves = pd.DataFrame(curve_rows)
    if metrics.shape[0] != 9 or nodes.shape[0] != 45 or curves.shape[0] != 4509:
        raise ValueError("Unexpected interpolation output size")
    metrics.to_pickle(METRICS)
    nodes.to_pickle(NODES)
    curves.to_pickle(CURVES)

    print(metrics[[
        "DESCRIPTOR", "QUADRATIC_FIT_RMSE_EV",
        "QUADRATIC_MAX_ABS_ERROR_EV",
    ]].to_string(index=False, float_format=lambda value: f"{value:.4f}"))
    worst = metrics.loc[metrics["QUADRATIC_MAX_ABS_ERROR_EV"].idxmax()]
    print(
        f"Largest quadratic error: {worst['QUADRATIC_MAX_ABS_ERROR_EV']:.4f} eV "
        f"for {worst['DESCRIPTOR']}; Lagrange error at nodes < 1e-10 eV."
    )


if __name__ == "__main__":
    main()
