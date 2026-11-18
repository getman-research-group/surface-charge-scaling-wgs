"""Archive the full-precision sources for SI S4.2.1--S4.2.4 and S5.

This exporter reads the saved analysis CSVs used to build the SI tables. Run
it from the repository root while the original SR work folder is available:

    python3 EXPORT_SI_RESULT_PKLS.PY --source-root /path/to/SR

It writes 36 PKLs in each two-descriptor model folder and one PKL in each BEP
folder. The two-descriptor LOMO values are nested metal-out scores; the BEP
LOMO values refit the same linear form without one metal. Equations retain
full-precision coefficients and use the adsorbate-energy notation.
"""

from __future__ import annotations

import argparse
import json
import re
from importlib.machinery import SourceFileLoader
from importlib.util import module_from_spec, spec_from_loader
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import sympy as sp


ROOT = Path(__file__).resolve().parent
ADSORBATES = (
    "H*", "O*", "OH*", "CH*", "CO* FCC", "CO2*", "H2O*", "HCO*", "COOH*"
)
TRANSITION_STATES = ("HO--H", "O--CO", "H--CO2", "OC--OH")
TARGETS = ADSORBATES + TRANSITION_STATES
PAIRS = tuple(combinations(ADSORBATES, 2))
NAMES = {
    "H*": r"\mathrm{H}^{*}", "O*": r"\mathrm{O}^{*}",
    "OH*": r"\mathrm{OH}^{*}", "CH*": r"\mathrm{CH}^{*}",
    "CO* FCC": r"\mathrm{CO}^{*}_{\mathrm{FCC}}",
    "CO2*": r"\mathrm{CO}_{2}^{*}",
    "H2O*": r"\mathrm{H}_{2}\mathrm{O}^{*}",
    "HCO*": r"\mathrm{HCO}^{*}", "COOH*": r"\mathrm{COOH}^{*}",
    "HO--H": r"\mathrm{HO{-}H}", "O--CO": r"\mathrm{O{-}CO}",
    "H--CO2": r"\mathrm{H{-}CO}_{2}",
    "OC--OH": r"\mathrm{OC{-}OH}",
}
X0, X1 = sp.symbols("X0 X1")
KEY = ["descriptor1", "descriptor2", "target"]
COLUMNS = [
    "MODEL", "DATASET", "DESCRIPTOR_1", "DESCRIPTOR_2", "TARGET",
    "TARGET_TYPE", "N", "EQUATION", "COMPLEXITY", "FIT_R2",
    "FIT_RMSE_EV", "LOMO_RMSE_EV", "PREFERRED_CRITERIA_MET",
]


def energy(name: str) -> str:
    kind = "act" if name in TRANSITION_STATES else "ad"
    return rf"\Delta E^{{\mathrm{{{kind}}}}}_{{{NAMES[name]}}}"


def equation(target: str, pair: tuple[str, str], expression: sp.Expr) -> str:
    rhs = sp.latex(expression, symbol_names={
        X0: energy(pair[0]), X1: energy(pair[1])
    })
    return f"{energy(target)} = {rhs}"


def pair_filename(pair: tuple[str, str]) -> str:
    slug = lambda name: re.sub(r"[^A-Z0-9]+", "_", name.upper()).strip("_")
    return f"RESULTS_{slug(pair[0])}_{slug(pair[1])}.PKL"


def as_float(value: float) -> sp.Float:
    return sp.Float(format(float(value), ".17g"))


def linear_expression(row: pd.Series) -> sp.Expr:
    return (
        as_float(row.coefficient1) * X0
        + as_float(row.coefficient2) * X1
        + as_float(row.intercept)
    )


def neutral_symbolic_expression(row: pd.Series) -> sp.Expr:
    terms = json.loads(row.terms)
    coefficients = json.loads(row.coefficients)
    if len(terms) != len(coefficients):
        raise ValueError("SyMANTIC term and coefficient counts differ")
    expression = as_float(row.intercept)
    for coefficient, term in zip(coefficients, terms):
        expression += as_float(coefficient) * sp.sympify(
            term, locals={"X0": X0, "X1": X1}
        )
    return expression


def full_symbolic_expression(row: pd.Series) -> sp.Expr:
    return sp.sympify(
        row.equation_X0_X1, locals={"X0": X0, "X1": X1}
    )


def load_family(source: Path, family: str) -> pd.DataFrame:
    if family == "LINEAR_NEUTRAL":
        data = pd.read_csv(source / "catmap_neutral_sparse_pairs_results/all_fits.csv")
        data["N"] = data.n_metals
        data["FIT_R2"] = data.fit_R2
        data["FIT_RMSE_EV"] = data.fit_RMSE_eV
        data["LOMO_RMSE_EV"] = data.LOMO_RMSE_eV
        data["EXPRESSION"] = data.apply(linear_expression, axis=1)
    elif family == "LINEAR_ALL_CHARGES":
        data = pd.read_csv(source / "catmap_all_charges_sparse_pairs_results/all_fits.csv")
        data["N"] = data.n_observations
        data["FIT_R2"] = data.fit_R2
        data["FIT_RMSE_EV"] = data.fit_RMSE_eV
        data["LOMO_RMSE_EV"] = data.LOMO_RMSE_eV
        data["EXPRESSION"] = data.apply(linear_expression, axis=1)
    elif family == "NONLINEAR_NEUTRAL":
        data = pd.read_csv(source / "neutral_symantic_36_pairs_results/selected_equations.csv")
        nested = pd.read_csv(
            source / "neutral_symantic_nested_LOMO_results/nested_LOMO_target_metrics.csv"
        )
        if not nested.complete.all():
            raise ValueError("Neutral SyMANTIC nested LOMO results are incomplete")
        # The selected-equation CSV's LOMO column is from an earlier fixed-form
        # analysis. Replace it with the SI's nested metal-out score.
        data = data.drop(columns="LOMO_RMSE_eV").merge(
            nested[KEY + ["target_type", "nested_LOMO_RMSE_eV"]],
            on=KEY, validate="one_to_one"
        )
        data["N"] = data.n_neutral
        data["FIT_R2"] = data.fit_R2
        data["FIT_RMSE_EV"] = data.fit_RMSE_eV
        data["LOMO_RMSE_EV"] = data.nested_LOMO_RMSE_eV
        data["EXPRESSION"] = data.apply(neutral_symbolic_expression, axis=1)
    elif family == "NONLINEAR_ALL_CHARGES":
        data = pd.read_csv(
            source / "symantic_all_charges_36_pair_tables/selected_equations_with_LOMO.csv"
        )
        if not data.equation_verified.all():
            raise ValueError("Unverified full-dataset SyMANTIC equation")
        data["N"] = data.n_observations
        data["FIT_R2"] = data.direct_fit_R2
        data["FIT_RMSE_EV"] = data.direct_fit_RMSE_eV
        data["LOMO_RMSE_EV"] = data.LOMO_RMSE_eV
        data["EXPRESSION"] = data.apply(full_symbolic_expression, axis=1)
    else:
        raise ValueError(f"Unknown family: {family}")
    if len(data) != 396 or data[KEY].duplicated().any():
        raise ValueError(f"Expected 396 unique equations for {family}")
    if set(map(tuple, data[["descriptor1", "descriptor2"]].drop_duplicates().values)) != set(PAIRS):
        raise ValueError(f"Expected all 36 adsorbate-energy pairs for {family}")
    return data


def pooled_rows(rows: pd.DataFrame) -> pd.DataFrame:
    summaries = []
    for label, subset in (
        ("POOLED_ALL", rows),
        ("POOLED_ADSORBATES", rows.loc[rows.TARGET_TYPE == "adsorbate"]),
        ("POOLED_TS", rows.loc[rows.TARGET_TYPE == "TS"]),
    ):
        n = int(subset.N.sum())
        base = {column: np.nan for column in COLUMNS}
        base.update({
            "MODEL": rows.MODEL.iloc[0], "DATASET": rows.DATASET.iloc[0],
            "DESCRIPTOR_1": rows.DESCRIPTOR_1.iloc[0],
            "DESCRIPTOR_2": rows.DESCRIPTOR_2.iloc[0],
            "TARGET": label, "TARGET_TYPE": "pooled", "N": n,
            "FIT_RMSE_EV": float(np.sqrt(
                np.dot(subset.N, subset.FIT_RMSE_EV ** 2) / n
            )),
            "LOMO_RMSE_EV": float(np.sqrt(
                np.dot(subset.N, subset.LOMO_RMSE_EV ** 2) / n
            )),
        })
        summaries.append(base)
    return pd.DataFrame(summaries, columns=COLUMNS)


def export_family(source: Path, family: str) -> list[Path]:
    data = load_family(source, family)
    outputs = []
    for pair in PAIRS:
        subset = data.loc[
            data.descriptor1.eq(pair[0]) & data.descriptor2.eq(pair[1])
        ].copy()
        expected = [target for target in TARGETS if target not in pair]
        if len(subset) != 11 or set(subset.target) != set(expected):
            raise ValueError(f"Wrong target coverage for {family}: {pair}")
        subset["ORDER"] = subset.target.map({name: i for i, name in enumerate(expected)})
        subset = subset.sort_values("ORDER")
        rows = []
        for _, item in subset.iterrows():
            preferred = (
                bool(item.passes_criteria) if family == "NONLINEAR_ALL_CHARGES"
                else np.nan
            )
            rows.append({
                "MODEL": "CATMAP linear" if family.startswith("LINEAR") else "SyMANTIC nonlinear",
                "DATASET": "neutral-only" if family.endswith("NEUTRAL") else "all charges",
                "DESCRIPTOR_1": pair[0], "DESCRIPTOR_2": pair[1],
                "TARGET": item.target, "TARGET_TYPE": item.target_type,
                "N": int(item.N),
                "EQUATION": equation(item.target, pair, item.EXPRESSION),
                "COMPLEXITY": float(item.complexity),
                "FIT_R2": float(item.FIT_R2),
                "FIT_RMSE_EV": float(item.FIT_RMSE_EV),
                "LOMO_RMSE_EV": float(item.LOMO_RMSE_EV),
                "PREFERRED_CRITERIA_MET": preferred,
            })
        frame = pd.DataFrame(rows, columns=COLUMNS)
        frame = pd.concat([frame, pooled_rows(frame)], ignore_index=True)
        path = ROOT / family / pair_filename(pair)
        frame.to_pickle(path)
        outputs.append(path)
    return outputs


def load_bep_common():
    path = ROOT / "BEP_COMMON.PY"
    spec = spec_from_loader("BEP_COMMON", SourceFileLoader("BEP_COMMON", str(path)))
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def export_bep(neutral_only: bool) -> Path:
    bep = load_bep_common()
    frame = pd.read_pickle(ROOT / "ADSORPTION_TS_ENERGIES.PKL")
    if neutral_only:
        frame = frame.loc[np.isclose(frame.Charge, 0)]
    records = []
    for target in TRANSITION_STATES:
        item = bep.fit_target(frame, target)
        records.append({
            "MODEL": "BEP linear",
            "DATASET": "neutral-only" if neutral_only else "all charges",
            "TARGET": target, "N": int(item["N"]),
            "EQUATION": item["EQUATION"],
            "FIT_R2": float(item["R2"]),
            "FIT_RMSE_EV": float(item["RMSE_EV"]),
            "LOMO_RMSE_EV": float(item["LOMO_RMSE_EV"]),
        })
    result = pd.DataFrame(records)
    n = int(result.N.sum())
    result.loc[len(result)] = {
        "MODEL": "BEP linear", "DATASET": result.DATASET.iloc[0],
        "TARGET": "POOLED_TS", "N": n, "EQUATION": "",
        "FIT_R2": np.nan,
        "FIT_RMSE_EV": float(np.sqrt(np.dot(result.N, result.FIT_RMSE_EV ** 2) / n)),
        "LOMO_RMSE_EV": float(np.sqrt(np.dot(result.N, result.LOMO_RMSE_EV ** 2) / n)),
    }
    folder = "BEP_NEUTRAL" if neutral_only else "BEP_ALL_CHARGES"
    path = ROOT / folder / "RESULTS_ALL.PKL"
    result.to_pickle(path)
    return path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source-root", type=Path, default=ROOT.parent,
        help="Original SR work folder with the saved SI analysis CSVs",
    )
    args = parser.parse_args()
    outputs = []
    for family in (
        "LINEAR_NEUTRAL", "LINEAR_ALL_CHARGES",
        "NONLINEAR_NEUTRAL", "NONLINEAR_ALL_CHARGES",
    ):
        outputs.extend(export_family(args.source_root, family))
    outputs.extend((export_bep(True), export_bep(False)))
    if len(outputs) != 146 or len(set(outputs)) != 146:
        raise ValueError("Expected exactly 146 distinct result PKLs")
    print(f"Wrote {len(outputs)} PKLs: 36 per two-descriptor family and 2 BEP files")


if __name__ == "__main__":
    main()
