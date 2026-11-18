"""Shared data, scoring, and fitting logic for the four two-descriptor models.

Each model has its own runnable script and folder. This module keeps the
descriptor input, target list, reported metrics, and equation notation
identical across those scripts.
"""

from __future__ import annotations

import argparse
import csv
import io
import os
import random
import re
import warnings
from contextlib import redirect_stdout
from dataclasses import dataclass
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
DATA_FILE = ROOT / "ADSORPTION_TS_ENERGIES.PKL"
ADSORBATES = (
    "H*", "O*", "OH*", "CH*", "CO* FCC", "CO2*", "H2O*", "HCO*", "COOH*"
)
TRANSITION_STATES = ("HO--H", "O--CO", "H--CO2", "OC--OH")
SPECIES_TEX = {
    "H*": r"\mathrm{H}^{*}",
    "O*": r"\mathrm{O}^{*}",
    "OH*": r"\mathrm{OH}^{*}",
    "CH*": r"\mathrm{CH}^{*}",
    "CO* FCC": r"\mathrm{CO}^{*}_{\mathrm{FCC}}",
    "CO2*": r"\mathrm{CO}_{2}^{*}",
    "H2O*": r"\mathrm{H}_{2}\mathrm{O}^{*}",
    "HCO*": r"\mathrm{HCO}^{*}",
    "COOH*": r"\mathrm{COOH}^{*}",
    "HO--H": r"\mathrm{HO{-}H}",
    "O--CO": r"\mathrm{O{-}CO}",
    "H--CO2": r"\mathrm{H{-}CO}_{2}",
    "OC--OH": r"\mathrm{OC{-}OH}",
}


@dataclass
class TargetResult:
    target: str
    equation: str
    actual: np.ndarray
    fitted: np.ndarray
    heldout: np.ndarray


def energy_tex(name: str) -> str:
    kind = "act" if name in TRANSITION_STATES else "ad"
    return rf"\Delta E^{{\mathrm{{{kind}}}}}_{{{SPECIES_TEX[name]}}}"


def score_r2(actual: np.ndarray, predicted: np.ndarray) -> float:
    actual = np.asarray(actual, dtype=float)
    predicted = np.asarray(predicted, dtype=float)
    denominator = np.sum((actual - actual.mean()) ** 2)
    numerator = np.sum((actual - predicted) ** 2)
    if denominator == 0:
        return 1.0 if numerator == 0 else float("-inf")
    return float(1 - numerator / denominator)


def score_rmse(actual: np.ndarray, predicted: np.ndarray) -> float:
    actual = np.asarray(actual, dtype=float)
    predicted = np.asarray(predicted, dtype=float)
    if actual.shape != predicted.shape or not np.isfinite(predicted).all():
        raise ValueError("Predictions are missing, nonfinite, or incorrectly shaped")
    return float(np.sqrt(np.mean((actual - predicted) ** 2)))


def load_data(neutral_only: bool) -> pd.DataFrame:
    frame = pd.read_pickle(DATA_FILE)
    if len(frame) != 25 or frame.Metal.nunique() != 5:
        raise ValueError("Expected 25 observations from five metals")
    if frame[list(ADSORBATES)].isna().any().any():
        raise ValueError("The nine candidate adsorption energies must be complete")
    if frame[list(ADSORBATES) + list(TRANSITION_STATES)].notna().sum().sum() != 299:
        raise ValueError("Unexpected adsorption-energy or TS-energy availability")
    if neutral_only:
        frame = frame.loc[np.isclose(frame.Charge, 0)]
        if len(frame) != 5:
            raise ValueError("Expected one neutral observation per metal")
    frame = frame.sort_values(["Metal", "Charge"]).reset_index(drop=True)
    frame["SURFACE_ID"] = (
        frame.Metal if neutral_only
        else frame.Metal + "|" + frame.Charge.map(lambda value: f"{value:+.2f}")
    )
    if not frame.SURFACE_ID.is_unique:
        raise ValueError("Metal-charge surface identifiers are not unique")
    return frame


def read_pair(arguments: argparse.Namespace) -> tuple[str, str]:
    if arguments.descriptors is None:
        print("Available adsorption-energy inputs: " + ", ".join(ADSORBATES))
        raw = input("Enter two species labels, separated by a comma: ")
        entered = [part.strip() for part in raw.split(",")]
    else:
        entered = list(arguments.descriptors)
    lookup = {name.upper(): name for name in ADSORBATES}
    if len(entered) != 2 or any(value.upper() not in lookup for value in entered):
        raise ValueError("Choose two labels from: " + ", ".join(ADSORBATES))
    pair = tuple(sorted((lookup[value.upper()] for value in entered),
                        key=ADSORBATES.index))
    if pair[0] == pair[1]:
        raise ValueError("The two descriptor adsorption energies must differ")
    return pair


def arguments_for(model_name: str, output_folder: Path,
                  argv: list[str] | None = None):
    parser = argparse.ArgumentParser(description=model_name)
    parser.add_argument(
        "--descriptors", nargs=2, metavar=("FIRST", "SECOND"),
        help="Supply two adsorbate labels instead of answering the prompt",
    )
    parser.add_argument(
        "--output", type=Path,
        help="CSV path; default is this model's folder",
    )
    parser.add_argument(
        "--target", choices=ADSORBATES + TRANSITION_STATES,
        help="Fit one target for a short diagnostic run; default fits all 11",
    )
    args = parser.parse_args(argv)
    pair = read_pair(args)
    targets = [name for name in ADSORBATES + TRANSITION_STATES if name not in pair]
    if args.target:
        if args.target in pair:
            raise ValueError("A supplied descriptor cannot also be a target")
        targets = [args.target]
    slug = lambda value: re.sub(r"[^A-Z0-9]+", "_", value.upper()).strip("_")
    suffix = ("_" + slug(args.target)) if args.target else ""
    output = args.output or output_folder / (
        f"RESULTS_{slug(pair[0])}_{slug(pair[1])}{suffix}.CSV"
    )
    return pair, targets, output


def write_results(model: str, pair: tuple[str, str], results: list[TargetResult],
                  output: Path) -> None:
    fields = ("MODEL", "INPUT_1", "INPUT_2", "SPECIES", "EQUATION",
              "R2", "RMSE_EV", "LOMO_RMSE_EV")
    rows = []
    for result in results:
        rows.append({
            "MODEL": model, "INPUT_1": energy_tex(pair[0]),
            "INPUT_2": energy_tex(pair[1]), "SPECIES": result.target,
            "EQUATION": result.equation,
            "R2": format(score_r2(result.actual, result.fitted), ".17g"),
            "RMSE_EV": format(score_rmse(result.actual, result.fitted), ".17g"),
            "LOMO_RMSE_EV": format(
                score_rmse(result.actual, result.heldout), ".17g"
            ),
        })

    if len(results) == 11:
        scopes = (
            ("POOLED_ALL", results),
            ("POOLED_ADSORBATES", [
                result for result in results if result.target in ADSORBATES
            ]),
            ("POOLED_TS", [
                result for result in results
                if result.target in TRANSITION_STATES
            ]),
        )
        for label, scope in scopes:
            actual = np.concatenate([item.actual for item in scope])
            fitted = np.concatenate([item.fitted for item in scope])
            heldout = np.concatenate([item.heldout for item in scope])
            rows.append({
                "MODEL": model, "INPUT_1": energy_tex(pair[0]),
                "INPUT_2": energy_tex(pair[1]), "SPECIES": label,
                "EQUATION": "", "R2": "",
                "RMSE_EV": format(score_rmse(actual, fitted), ".17g"),
                "LOMO_RMSE_EV": format(
                    score_rmse(actual, heldout), ".17g"
                ),
            })

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {output}")


def linear_equation(target: str, subset: tuple[str, ...],
                    coefficients: np.ndarray) -> str:
    pieces = []
    for index, descriptor in enumerate(subset):
        coefficient = float(coefficients[index])
        sign = "-" if coefficient < 0 else "+"
        term = f"{format(abs(coefficient), '.17g')}{energy_tex(descriptor)}"
        pieces.append(("-" if sign == "-" else "") + term if index == 0
                      else f" {sign} {term}")
    intercept = float(coefficients[-1])
    sign = "-" if intercept < 0 else "+"
    pieces.append(f" {sign} {format(abs(intercept), '.17g')}")
    return f"{energy_tex(target)} = {''.join(pieces)}"


def symbolic_equation(target: str, expression, pair: tuple[str, str]) -> str:
    import sympy as sp
    x0, x1 = sp.symbols("X0 X1")
    rhs = sp.latex(
        expression,
        symbol_names={x0: energy_tex(pair[0]), x1: energy_tex(pair[1])},
    )
    return f"{energy_tex(target)} = {rhs}"


def _catmap_candidates(frame: pd.DataFrame, pair: tuple[str, str],
                       targets: list[str]) -> dict[str, list[dict]]:
    try:
        from catmap.functions import parse_constraint, scaling_coefficient_matrix
    except ImportError as exc:
        raise RuntimeError(
            "CATMAP is required. Install python-catmap 0.3.1 or expose its "
            "installation through PYTHONPATH."
        ) from exc
    found = {target: [] for target in targets}
    surfaces = frame.SURFACE_ID.tolist()
    parameters = {
        target: [None if pd.isna(value) else float(value)
                 for value in frame[target]]
        for target in targets
    }
    for subset in ((pair[0],), (pair[1],), pair):
        descriptors = {
            str(row.SURFACE_ID): [float(row[name]) for name in subset]
            for _, row in frame.iterrows()
        }
        bounds = [
            parse_constraint([None] * (len(subset) + 1), target)
            for target in targets
        ]
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", message="Number of energies specified is less than"
            )
            matrix = scaling_coefficient_matrix(
                parameters, descriptors, surfaces, parameter_names=targets,
                coeff_mins=[low for low, _ in bounds],
                coeff_maxs=[high for _, high in bounds],
            )
        coefficients = np.asarray(matrix.T, dtype=float)
        expected = (len(targets), len(subset) + 1)
        if coefficients.shape != expected:
            raise ValueError(
                f"CATMAP returned shape {coefficients.shape}, expected {expected}"
            )
        for target, coefficient in zip(targets, coefficients):
            observed = frame.loc[frame[target].notna()]
            x = np.column_stack((
                observed[list(subset)].to_numpy(float),
                np.ones(len(observed)),
            ))
            actual = observed[target].to_numpy(float)
            predicted = x @ coefficient
            found[target].append({
                "subset": subset, "coefficients": coefficient,
                "r2": score_r2(actual, predicted),
                "rmse": score_rmse(actual, predicted),
            })
    return found


def _choose_catmap(candidates: list[dict]) -> dict:
    best_r2 = max(item["r2"] for item in candidates)
    eligible = [
        item for item in candidates
        if item["r2"] >= best_r2 - 0.02 - 1e-12
    ]
    return min(eligible, key=lambda item: (
        len(item["subset"]), -item["r2"], item["rmse"], item["subset"]
    ))


def _predict_catmap(candidate: dict, frame: pd.DataFrame) -> np.ndarray:
    features = np.column_stack((
        frame[list(candidate["subset"])].to_numpy(float),
        np.ones(len(frame)),
    ))
    return features @ candidate["coefficients"]


def run_catmap(neutral_only: bool, model_name: str, output_folder: Path,
               argv: list[str] | None = None) -> Path:
    pair, targets, output = arguments_for(model_name, output_folder, argv)
    frame = load_data(neutral_only)
    full = {
        target: _choose_catmap(candidates)
        for target, candidates in _catmap_candidates(frame, pair, targets).items()
    }
    heldout = {target: np.full(len(frame), np.nan) for target in targets}
    for metal in sorted(frame.Metal.unique()):
        train = frame.loc[frame.Metal != metal]
        test = frame.loc[frame.Metal == metal]
        fold = {
            target: _choose_catmap(candidates)
            for target, candidates in
            _catmap_candidates(train, pair, targets).items()
        }
        for target in targets:
            observed = test.loc[test[target].notna()]
            if len(observed):
                heldout[target][observed.index] = _predict_catmap(
                    fold[target], observed
                )

    results = []
    for target in targets:
        observed = frame.loc[frame[target].notna()]
        selected = full[target]
        results.append(TargetResult(
            target=target,
            equation=linear_equation(
                target, selected["subset"], selected["coefficients"]
            ),
            actual=observed[target].to_numpy(float),
            fitted=_predict_catmap(selected, observed),
            heldout=heldout[target][observed.index],
        ))
    write_results(model_name, pair, results, output)
    return output


def _predict_symbolic(expression, frame: pd.DataFrame,
                      pair: tuple[str, str]) -> np.ndarray:
    import sympy as sp
    x0, x1 = sp.symbols("X0 X1")
    function = sp.lambdify((x0, x1), expression, modules="numpy")
    with np.errstate(all="ignore"):
        predicted = np.asarray(
            function(
                frame[pair[0]].to_numpy(float),
                frame[pair[1]].to_numpy(float),
            ),
            dtype=float,
        )
    if predicted.ndim == 0:
        predicted = np.full(len(frame), float(predicted))
    if predicted.shape != (len(frame),) or not np.isfinite(predicted).all():
        raise ValueError("SyMANTIC equation gives nonfinite predictions")
    return predicted


def _symbolic_candidates(observed: pd.DataFrame, pair: tuple[str, str],
                         target: str, neutral_only: bool) -> list[dict]:
    import sympy as sp
    try:
        import torch
        from symantic import SymanticModel
    except ImportError as exc:
        raise RuntimeError(
            "SyMANTIC 1.1.0.6 and PyTorch are required. Install them or "
            "expose their installation through PYTHONPATH."
        ) from exc

    data = observed[[target, *pair]].rename(columns={
        target: "Target", pair[0]: "X0", pair[1]: "X1"
    }).reset_index(drop=True)
    random.seed(0)
    np.random.seed(0)
    torch.manual_seed(0)
    model = SymanticModel(
        df=data, operators=["+", "-", "*", "/"],
        n_expansion=None, n_term=2 if neutral_only else 3,
        sis_features=10, metrics=[0.0, 1.1], device="cpu",
        disp=False, pareto=False,
    )
    with warnings.catch_warnings(), redirect_stdout(io.StringIO()), patch(
        "builtins.input", return_value="no"
    ):
        warnings.simplefilter("ignore", category=RuntimeWarning)
        _, pareto = model.fit()
    selected_rows = (
        model.final_df if neutral_only else pareto
    )
    selected_rows = selected_rows.loc[
        selected_rows.Complexity <= (10 if neutral_only else 50) + 1e-8
    ]
    x0, x1 = sp.symbols("X0 X1")
    rows = []
    for index, summary in selected_rows.iterrows():
        raw = model.final_df.loc[index]
        if isinstance(raw, pd.DataFrame):
            raise ValueError("SyMANTIC returned duplicate candidate indices")
        terms = raw.Equations
        if not isinstance(terms, (list, tuple, np.ndarray)):
            terms = [str(terms)]
        terms = [str(term) for term in terms]
        coefficients = [float(value) for value in raw.Coefficients[:len(terms)]]
        intercept = float(raw.Intercepts)
        expression = sp.Float(format(intercept, ".17g"))
        for coefficient, term in zip(coefficients, terms):
            basis = sp.sympify(term, locals={"X0": x0, "X1": x1})
            expression += sp.Float(format(coefficient, ".17g")) * basis
        try:
            predicted = _predict_symbolic(expression, data, ("X0", "X1"))
            actual = data.Target.to_numpy(float)
            rmse = score_rmse(actual, predicted)
            if not np.isclose(rmse, float(summary.Loss),
                              rtol=1e-5, atol=1e-6):
                continue
            rows.append({
                "expression": expression, "terms": terms,
                "complexity": float(summary.Complexity),
                "r2": score_r2(actual, predicted), "rmse": rmse,
            })
        except (ValueError, TypeError, ZeroDivisionError, OverflowError):
            continue
    if not rows:
        raise ValueError(
            f"No verified SyMANTIC candidate for {target} using {pair}"
        )
    return rows


def _neutral_order(candidates: list[dict], target: str) -> list[dict]:
    r2_floor = 0.65 if target in TRANSITION_STATES else 0.70
    rmse_ceiling = 0.46 if target in TRANSITION_STATES else 0.60
    eligible = [
        item for item in candidates
        if item["r2"] >= r2_floor and item["rmse"] < rmse_ceiling
    ]
    pool = eligible if eligible else candidates
    best_r2 = max(item["r2"] for item in pool)
    best_rmse = min(item["rmse"] for item in pool)
    near = [
        item for item in pool
        if best_r2 - item["r2"] <= 0.02 + 1e-12
        or item["rmse"] - best_rmse <= 0.01 + 1e-12
    ]
    rest = [item for item in pool if item not in near]
    return (
        sorted(near, key=lambda item: (
            item["complexity"], item["rmse"]
        ))
        + sorted(rest, key=lambda item: (
            -item["r2"], item["complexity"]
        ))
    )


def _full_order(candidates: list[dict], target: str) -> list[dict]:
    is_ts = target in TRANSITION_STATES
    r2_floor = 0.65 if is_ts else 0.70
    complexity_ceiling = 25.0 if is_ts else 14.0
    rmse_ceiling = 0.46 if is_ts else 0.60
    tiers = (
        [item for item in candidates
         if item["r2"] > r2_floor
         and item["complexity"] < complexity_ceiling
         and item["rmse"] < rmse_ceiling],
        [item for item in candidates
         if item["complexity"] < complexity_ceiling
         and item["rmse"] < rmse_ceiling],
        [item for item in candidates
         if item["complexity"] < complexity_ceiling],
        candidates,
    )
    pool = next((tier for tier in tiers if tier), None)
    if pool is None:
        raise ValueError("No verified SyMANTIC candidate")
    best_r2 = max(item["r2"] for item in pool)
    near = [
        item for item in pool
        if item["r2"] >= best_r2 - 0.02 - 1e-12
    ]
    return sorted(near, key=lambda item: (
        item["complexity"], item["rmse"], -item["r2"]
    ))


def _neutral_identifiable(candidate: dict, frame: pd.DataFrame,
                          pair: tuple[str, str], all_folds: bool) -> bool:
    import sympy as sp
    x0, x1 = sp.symbols("X0 X1")
    columns = []
    for term in candidate["terms"]:
        expression = sp.sympify(term, locals={"X0": x0, "X1": x1})
        columns.append(_predict_symbolic(expression, frame, pair))
    design = np.column_stack([*columns, np.ones(len(frame))])
    if design.shape[1] > len(frame):
        return False
    if np.linalg.matrix_rank(design) != design.shape[1]:
        return False
    if all_folds:
        for metal in frame.Metal.unique():
            train = frame.Metal.ne(metal).to_numpy()
            if train.sum() < design.shape[1]:
                return False
            if np.linalg.matrix_rank(design[train]) != design.shape[1]:
                return False
    return True


def _select_symbolic(observed: pd.DataFrame, pair: tuple[str, str],
                     target: str, neutral_only: bool,
                     final_fit: bool) -> dict:
    candidates = _symbolic_candidates(observed, pair, target, neutral_only)
    ordered = (
        _neutral_order(candidates, target) if neutral_only
        else _full_order(candidates, target)
    )
    if not neutral_only:
        return ordered[0]
    for candidate in ordered:
        try:
            if _neutral_identifiable(
                candidate, observed, pair, all_folds=final_fit
            ):
                return candidate
        except (ValueError, TypeError, ZeroDivisionError, OverflowError):
            pass
    raise ValueError(
        f"No identifiable neutral SyMANTIC candidate for {target}"
    )


def run_symantic(neutral_only: bool, model_name: str, output_folder: Path,
                 argv: list[str] | None = None) -> Path:
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")
    try:
        import torch
    except ImportError as exc:
        raise RuntimeError("PyTorch is required for SyMANTIC") from exc
    torch.set_num_threads(1)
    try:
        torch.set_num_interop_threads(1)
    except RuntimeError:
        pass

    pair, targets, output = arguments_for(model_name, output_folder, argv)
    frame = load_data(neutral_only)
    results = []
    for number, target in enumerate(targets, 1):
        observed = frame.loc[frame[target].notna()].copy()
        selected = _select_symbolic(
            observed, pair, target, neutral_only, final_fit=True
        )
        fitted = _predict_symbolic(
            selected["expression"], observed, pair
        )
        heldout = np.full(len(observed), np.nan)
        for metal in sorted(observed.Metal.unique()):
            train = observed.loc[observed.Metal != metal]
            test = observed.loc[observed.Metal == metal]
            fold = _select_symbolic(
                train, pair, target, neutral_only, final_fit=False
            )
            heldout[observed.Metal.eq(metal).to_numpy()] = (
                _predict_symbolic(fold["expression"], test, pair)
            )
        results.append(TargetResult(
            target=target,
            equation=symbolic_equation(
                target, selected["expression"], pair
            ),
            actual=observed[target].to_numpy(float),
            fitted=fitted, heldout=heldout,
        ))
        print(f"[{number}/{len(targets)}] {target} complete", flush=True)
    write_results(model_name, pair, results, output)
    return output
