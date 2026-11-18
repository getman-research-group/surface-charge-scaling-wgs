"""Reaction-energy construction and BEP fitting for the four WGS transitions.

The electronic gas references and reaction-energy directions reproduce the
original BEP analysis in SR/ts_reaction_energies/calculate.py. HO--H uses the
reverse water-formation reaction energy, as in the SI BEP equation.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
DATA_FILE = ROOT / "ADSORPTION_TS_ENERGIES.PKL"
GAS_ENERGY_EV = {
    "H2": -6.726328,
    "O2": -9.853529,
    "OH": -7.549490,
    "CO": -14.795800,
    "CO2": -22.985849,
    "H2O": -14.218838,
    "COOH": -24.464805,
}
ADSORPTION_REFERENCE_EV = {
    "H*": GAS_ENERGY_EV["H2"] / 2,
    "O*": GAS_ENERGY_EV["O2"] / 2,
    "OH*": GAS_ENERGY_EV["OH"],
    "CO* FCC": GAS_ENERGY_EV["CO"],
    "CO2*": GAS_ENERGY_EV["CO2"],
    "H2O*": GAS_ENERGY_EV["H2O"],
    "COOH*": GAS_ENERGY_EV["COOH"],
}

# Each tuple is (reactant adsorbates, product adsorbates) in the forward direction.
REACTIONS = {
    "HO--H": (("OH*", "H*"), ("H2O*",)),
    "O--CO": (("CO* FCC", "O*"), ("CO2*",)),
    "H--CO2": (("COOH*",), ("CO2*", "H*")),
    "OC--OH": (("OH*", "CO* FCC"), ("COOH*",)),
}
TARGET_TEX = {
    "HO--H": r"\mathrm{HO{-}H}",
    "O--CO": r"\mathrm{O{-}CO}",
    "H--CO2": r"\mathrm{H{-}CO}_{2}",
    "OC--OH": r"\mathrm{OC{-}OH}",
}


def reaction_energies(frame: pd.DataFrame, target: str) -> np.ndarray:
    """Compute the reaction-energy predictor used by the SI BEP equations."""
    reactants, products = REACTIONS[target]
    required = set(reactants + products)
    if frame[list(required)].isna().any().any():
        raise ValueError(f"Missing adsorption energy required for {target}")
    gas_offset = (
        sum(ADSORPTION_REFERENCE_EV[name] for name in products)
        - sum(ADSORPTION_REFERENCE_EV[name] for name in reactants)
    )
    forward_energy = (
        frame[list(products)].sum(axis=1).to_numpy(float)
        - frame[list(reactants)].sum(axis=1).to_numpy(float)
        + gas_offset
    )
    # The SI fitted HO--H against water dissociation energy, the negative of
    # the forward water-formation energy shown in the manuscript.
    return -forward_energy if target == "HO--H" else forward_energy


def fit_line(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    design = np.column_stack((x, np.ones(len(x))))
    if len(x) < 3 or np.linalg.matrix_rank(design) != 2:
        raise ValueError("BEP fitting needs at least three distinct reaction energies")
    slope, intercept = np.linalg.lstsq(design, y, rcond=None)[0]
    return float(slope), float(intercept)


def rmse(actual: np.ndarray, predicted: np.ndarray) -> float:
    return float(np.sqrt(np.mean((actual - predicted) ** 2)))


def r2(actual: np.ndarray, predicted: np.ndarray) -> float:
    denominator = float(np.sum((actual - actual.mean()) ** 2))
    if denominator == 0:
        raise ValueError("R2 is undefined for a constant activation-energy target")
    return float(1 - np.sum((actual - predicted) ** 2) / denominator)


def equation(target: str, slope: float, intercept: float) -> str:
    subscript = TARGET_TEX[target]
    sign = "+" if intercept >= 0 else "-"
    return (
        rf"\Delta E^{{\mathrm{{act}}}}_{{{subscript}}} = "
        rf"{slope:.17g}\Delta E^{{\mathrm{{rxn}}}}_{{{subscript}}} "
        rf"{sign} {abs(intercept):.17g}"
    )


def fit_target(frame: pd.DataFrame, target: str) -> dict[str, str]:
    observed = frame.loc[frame[target].notna()].copy()
    if observed.Metal.nunique() != 4:
        raise ValueError(f"Expected four metals with {target} transition states")
    x = reaction_energies(observed, target)
    y = observed[target].to_numpy(float)
    metals = observed.Metal.to_numpy()
    slope, intercept = fit_line(x, y)
    fitted = slope * x + intercept
    heldout = np.empty(len(y), dtype=float)
    for metal in sorted(np.unique(metals)):
        test = metals == metal
        fold_slope, fold_intercept = fit_line(x[~test], y[~test])
        heldout[test] = fold_slope * x[test] + fold_intercept
    return {
        "TARGET": target,
        "N": str(len(observed)),
        "EQUATION": equation(target, slope, intercept),
        "R2": format(r2(y, fitted), ".17g"),
        "RMSE_EV": format(rmse(y, fitted), ".17g"),
        "LOMO_RMSE_EV": format(rmse(y, heldout), ".17g"),
    }


def run_bep(neutral_only: bool, output_folder: Path,
            argv: list[str] | None = None) -> Path:
    parser = argparse.ArgumentParser(
        description="Fit a BEP relation from DFT reaction and activation energies"
    )
    parser.add_argument(
        "--target", choices=(*REACTIONS, "ALL"),
        help="Transition-state label; omit to be prompted, or use ALL",
    )
    parser.add_argument("--output", type=Path, help="Output CSV path")
    args = parser.parse_args(argv)
    if args.target is None:
        print("Available TS targets: " + ", ".join(REACTIONS))
        target = input("Enter one TS target (or ALL): ").strip().upper()
    else:
        target = args.target
    lookup = {name.upper(): name for name in REACTIONS}
    if target != "ALL" and target not in lookup:
        parser.error("Choose one of " + ", ".join(REACTIONS) + ", or ALL")
    targets = list(REACTIONS) if target == "ALL" else [lookup[target]]

    frame = pd.read_pickle(DATA_FILE)
    if len(frame) != 25 or frame.Metal.nunique() != 5:
        raise ValueError("Expected five metals and five charge states per metal")
    if neutral_only:
        frame = frame.loc[np.isclose(frame.Charge, 0)]
        if len(frame) != 5:
            raise ValueError("Expected one neutral observation per metal")
    results = [fit_target(frame, name) for name in targets]
    safe_target = "ALL" if target == "ALL" else lookup[target].replace("--", "_")
    output = args.output or output_folder / f"RESULTS_{safe_target}.CSV"
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=("TARGET", "N", "EQUATION", "R2", "RMSE_EV", "LOMO_RMSE_EV"),
        )
        writer.writeheader()
        writer.writerows(results)
    print(f"Wrote {output}")
    return output
