"""Shared full-dataset fitting and nested metal-out validation for four inputs."""

from __future__ import annotations

import argparse
import csv
import io
import os
import random
import re
import warnings
from contextlib import redirect_stdout
from itertools import combinations
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
DATA = ROOT / "ADSORPTION_TS_ENERGIES.PKL"
ADSORBATES = ("H*", "O*", "OH*", "CH*", "CO* FCC", "CO2*", "H2O*", "HCO*", "COOH*")
TS = ("HO--H", "O--CO", "H--CO2", "OC--OH")
SPECIES_TEX = {
    "H*": r"\mathrm{H}^{*}", "O*": r"\mathrm{O}^{*}",
    "OH*": r"\mathrm{OH}^{*}", "CH*": r"\mathrm{CH}^{*}",
    "CO* FCC": r"\mathrm{CO}^{*}_{\mathrm{FCC}}",
    "CO2*": r"\mathrm{CO}_{2}^{*}",
    "H2O*": r"\mathrm{H}_{2}\mathrm{O}^{*}",
    "HCO*": r"\mathrm{HCO}^{*}", "COOH*": r"\mathrm{COOH}^{*}",
    "HO--H": r"\mathrm{HO{-}H}", "O--CO": r"\mathrm{O{-}CO}",
    "H--CO2": r"\mathrm{H{-}CO}_{2}", "OC--OH": r"\mathrm{OC{-}OH}",
}


def energy_tex(name):
    kind = "act" if name in TS else "ad"
    return rf"\Delta E^{{\mathrm{{{kind}}}}}_{{{SPECIES_TEX[name]}}}"


def rmse(actual, predicted):
    actual = np.asarray(actual, dtype=float)
    predicted = np.asarray(predicted, dtype=float)
    if actual.shape != predicted.shape or not np.isfinite(predicted).all():
        raise ValueError("Missing, nonfinite, or incorrectly shaped predictions")
    return float(np.sqrt(np.mean(np.square(actual - predicted))))


def r2(actual, predicted):
    actual = np.asarray(actual, dtype=float)
    predicted = np.asarray(predicted, dtype=float)
    denominator = float(np.square(actual - actual.mean()).sum())
    numerator = float(np.square(actual - predicted).sum())
    if denominator == 0:
        return 1.0 if numerator == 0 else float("-inf")
    return 1.0 - numerator / denominator


def read_data():
    data = pd.read_pickle(DATA).sort_values(["Metal", "Charge"]).reset_index(drop=True)
    if len(data) != 25 or data.Metal.nunique() != 5:
        raise ValueError("Expected 25 metal-charge observations for five metals")
    if data[list(ADSORBATES)].isna().any().any():
        raise ValueError("All nine candidate descriptor adsorption energies must exist")
    if int(data[list(ADSORBATES + TS)].notna().sum().sum()) != 299:
        raise ValueError("Unexpected target availability in the energy pickle")
    data["SURFACE_ID"] = data.Metal + "|" + data.Charge.map(lambda value: f"{value:+.2f}")
    if not data.SURFACE_ID.is_unique:
        raise ValueError("Metal-charge identifiers are not unique")
    return data


def parse_arguments(model_name, output_folder, argv=None):
    parser = argparse.ArgumentParser(description=model_name)
    parser.add_argument("--descriptors", nargs=4, metavar=("FIRST", "SECOND", "THIRD", "FOURTH"),
                        help="Four adsorbate labels instead of the prompt")
    parser.add_argument("--target", choices=ADSORBATES + TS,
                        help="Fit one target for a diagnostic check; default fits all nine")
    parser.add_argument("--output", type=Path, help="Output CSV path")
    args = parser.parse_args(argv)
    if args.descriptors is None:
        print("Available adsorption-energy inputs: " + ", ".join(ADSORBATES))
        entered = [part.strip() for part in input("Enter four species labels, separated by commas: ").split(",")]
    else:
        entered = list(args.descriptors)
    lookup = {name.upper(): name for name in ADSORBATES}
    if len(entered) != 4 or any(value.upper() not in lookup for value in entered):
        parser.error("Choose four labels from: " + ", ".join(ADSORBATES))
    # Preserve input order: SyMANTIC assigns X0--X3 in this order, as in the
    # original four-descriptor screen with CH* and HCO* first.
    descriptors = tuple(lookup[value.upper()] for value in entered)
    if len(set(descriptors)) != 4:
        parser.error("The four descriptor adsorption energies must differ")
    if args.target in descriptors:
        parser.error("A supplied descriptor cannot also be a target")
    targets = [name for name in ADSORBATES + TS if name not in descriptors]
    if args.target:
        targets = [args.target]

    def slug(value):
        return re.sub(r"[^A-Z0-9]+", "_", value.upper()).strip("_")

    suffix = "_" + slug(args.target) if args.target else ""
    output = args.output or output_folder / (
        "RESULTS_" + "_".join(slug(item) for item in descriptors) + suffix + ".CSV"
    )
    return descriptors, targets, output


def write_results(model_name, descriptors, results, output):
    fields = ("MODEL", "INPUT_1", "INPUT_2", "INPUT_3", "INPUT_4", "SPECIES",
              "EQUATION", "COMPLEXITY", "R2", "RMSE_EV", "LOMO_RMSE_EV")
    rows = []
    for result in results:
        row = {"MODEL": model_name, "SPECIES": result["target"],
               "EQUATION": result["equation"],
               "COMPLEXITY": format(result["complexity"], ".17g"),
               "R2": format(r2(result["actual"], result["fitted"]), ".17g"),
               "RMSE_EV": format(rmse(result["actual"], result["fitted"]), ".17g"),
               "LOMO_RMSE_EV": format(rmse(result["actual"], result["heldout"]), ".17g")}
        row.update({f"INPUT_{i}": energy_tex(descriptor)
                    for i, descriptor in enumerate(descriptors, 1)})
        rows.append(row)
    if len(results) == 9:
        scopes = (
            ("POOLED_ALL", results),
            ("POOLED_ADSORBATES", [item for item in results if item["target"] in ADSORBATES]),
            ("POOLED_TS", [item for item in results if item["target"] in TS]),
        )
        for label, subset in scopes:
            actual = np.concatenate([item["actual"] for item in subset])
            fitted = np.concatenate([item["fitted"] for item in subset])
            heldout = np.concatenate([item["heldout"] for item in subset])
            row = {"MODEL": model_name, "SPECIES": label,
                   "RMSE_EV": format(rmse(actual, fitted), ".17g"),
                   "LOMO_RMSE_EV": format(rmse(actual, heldout), ".17g")}
            row.update({f"INPUT_{i}": energy_tex(descriptor)
                        for i, descriptor in enumerate(descriptors, 1)})
            rows.append(row)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {output}")
    return output


def catmap_candidates(frame, descriptors, targets):
    try:
        from catmap.functions import parse_constraint, scaling_coefficient_matrix
    except ImportError as exc:
        raise RuntimeError("CATMAP is required (python-catmap==0.3.1)") from exc
    found = {target: [] for target in targets}
    surfaces = frame.SURFACE_ID.tolist()
    parameters = {target: [None if pd.isna(value) else float(value) for value in frame[target]]
                  for target in targets}
    for size in range(1, 5):
        for subset in combinations(descriptors, size):
            descriptor_values = {
                str(row.SURFACE_ID): [float(row[name]) for name in subset]
                for _, row in frame.iterrows()
            }
            bounds = [parse_constraint([None] * (len(subset) + 1), target) for target in targets]
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", message="Number of energies specified is less than")
                matrix = scaling_coefficient_matrix(
                    parameters, descriptor_values, surfaces, parameter_names=targets,
                    coeff_mins=[low for low, _ in bounds],
                    coeff_maxs=[high for _, high in bounds],
                )
            coefficients = np.asarray(matrix.T, dtype=float)
            if coefficients.shape != (len(targets), len(subset) + 1):
                raise ValueError("CATMAP returned an unexpected coefficient matrix")
            for target, coefficient in zip(targets, coefficients):
                observed = frame.loc[frame[target].notna()]
                features = np.column_stack((observed[list(subset)].to_numpy(float),
                                            np.ones(len(observed))))
                actual = observed[target].to_numpy(float)
                predicted = features @ coefficient
                found[target].append({
                    "subset": subset, "coefficients": coefficient,
                    "r2": r2(actual, predicted), "rmse": rmse(actual, predicted),
                    "complexity": float(len(subset)),
                })
    return found


def choose_catmap(candidates):
    best = max(item["r2"] for item in candidates)
    near = [item for item in candidates if item["r2"] >= best - 0.02 - 1e-12]
    return min(near, key=lambda item: (len(item["subset"]), -item["r2"],
                                       item["rmse"], item["subset"]))


def predict_catmap(model, frame):
    features = np.column_stack((frame[list(model["subset"])].to_numpy(float),
                                np.ones(len(frame))))
    return features @ model["coefficients"]


def catmap_equation(target, model):
    terms = []
    for coefficient, name in zip(model["coefficients"][:-1], model["subset"]):
        terms.append((float(coefficient), energy_tex(name)))
    terms.append((float(model["coefficients"][-1]), ""))
    rhs = ""
    for value, term in terms:
        sign = "-" if value < 0 else (" + " if rhs else "")
        if value < 0 and rhs:
            sign = " - "
        rhs += sign + repr(abs(value)) + term
    return f"{energy_tex(target)} = {rhs}"


def run_catmap(argv=None):
    descriptors, targets, output = parse_arguments(
        "ALL_CHARGES_LINEAR_CATMAP_4D", ROOT / "LINEAR_FOUR_DESCRIPTORS", argv
    )
    frame = read_data()
    selected = {target: choose_catmap(candidates)
                for target, candidates in catmap_candidates(frame, descriptors, targets).items()}
    heldout = {target: np.full(len(frame), np.nan) for target in targets}
    for metal in sorted(frame.Metal.unique()):
        train = frame.loc[frame.Metal.ne(metal)]
        test = frame.loc[frame.Metal.eq(metal)]
        fold = {target: choose_catmap(candidates)
                for target, candidates in catmap_candidates(train, descriptors, targets).items()}
        for target in targets:
            observed = test.loc[test[target].notna()]
            if len(observed):
                heldout[target][observed.index] = predict_catmap(fold[target], observed)
    results = []
    for target in targets:
        observed = frame.loc[frame[target].notna()]
        model = selected[target]
        results.append({
            "target": target, "equation": catmap_equation(target, model),
            "complexity": model["complexity"],
            "actual": observed[target].to_numpy(float),
            "fitted": predict_catmap(model, observed),
            "heldout": heldout[target][observed.index],
        })
    return write_results("ALL_CHARGES_LINEAR_CATMAP_4D", descriptors, results, output)


def symbolic_basis(terms, frame, descriptors):
    import sympy as sp
    symbols = sp.symbols("X0:4")
    x = frame[list(descriptors)].to_numpy(float)
    columns = []
    for term in terms:
        expression = sp.sympify(term, locals={str(symbol): symbol for symbol in symbols})
        function = sp.lambdify(symbols, expression, modules="numpy")
        with np.errstate(all="ignore"):
            values = np.asarray(function(*x.T), dtype=float)
        if values.ndim == 0:
            values = np.full(len(x), float(values))
        if values.shape != (len(x),) or not np.isfinite(values).all():
            raise ValueError(f"Nonfinite SyMANTIC basis term: {term}")
        columns.append(values)
    return np.column_stack(columns)


def symbolic_candidates(observed, descriptors, target):
    try:
        import torch
        from symantic import SymanticModel
    except ImportError as exc:
        raise RuntimeError("SyMANTIC 1.1.0.6 and PyTorch are required") from exc
    data = observed[[target, *descriptors]].rename(columns={
        target: "Target", **{descriptor: f"X{i}" for i, descriptor in enumerate(descriptors)}
    }).reset_index(drop=True)
    random.seed(0)
    np.random.seed(0)
    torch.manual_seed(0)
    model = SymanticModel(
        df=data, operators=["+", "-", "*", "/"], n_expansion=None,
        n_term=3, sis_features=10, metrics=[0.0, 1.1], device="cpu",
        disp=False, pareto=False,
    )
    with warnings.catch_warnings(), redirect_stdout(io.StringIO()), patch("builtins.input", return_value="no"):
        warnings.simplefilter("ignore", category=RuntimeWarning)
        model.fit()
    rows = []
    for _, candidate in model.final_df.iterrows():
        complexity = float(candidate.Complexity)
        if complexity > 15 + 1e-8:
            continue
        terms = candidate.Equations
        if not isinstance(terms, (list, tuple, np.ndarray)):
            terms = [str(terms)]
        terms = [str(term) for term in terms]
        coefficients = np.asarray(candidate.Coefficients[:len(terms)], dtype=float)
        intercept = float(candidate.Intercepts)
        try:
            matrix = symbolic_basis(terms, data, ("X0", "X1", "X2", "X3"))
            actual = data.Target.to_numpy(float)
            fitted = matrix @ coefficients + intercept
            error = rmse(actual, fitted)
            if not np.isclose(error, float(candidate.Loss), rtol=1e-5, atol=1e-6):
                continue
            rows.append({
                "terms": terms, "coefficients": coefficients, "intercept": intercept,
                "complexity": complexity, "r2": r2(actual, fitted), "rmse": error,
            })
        except (ValueError, TypeError, ZeroDivisionError, OverflowError):
            continue
    if not rows:
        raise ValueError(f"No verified SyMANTIC candidate within complexity 15 for {target}")
    return rows


def choose_symbolic(candidates, target):
    is_ts = target in TS
    r2_floor = 0.65 if is_ts else 0.70
    rmse_ceiling = 0.46 if is_ts else 0.60
    complexity_ok = lambda item: item["complexity"] < (15 + 1e-8 if is_ts else 14)
    tiers = (
        [item for item in candidates if item["r2"] > r2_floor
         and item["rmse"] < rmse_ceiling and complexity_ok(item)],
        [item for item in candidates if item["rmse"] < rmse_ceiling and complexity_ok(item)],
        [item for item in candidates if complexity_ok(item)],
        candidates,
    )
    pool = next(tier for tier in tiers if tier)
    best_r2 = max(item["r2"] for item in pool)
    near = [item for item in pool if item["r2"] >= best_r2 - 0.02 - 1e-12]
    return min(near, key=lambda item: (item["complexity"], item["rmse"], -item["r2"]))


def predict_symbolic(model, frame, descriptors):
    matrix = symbolic_basis(model["terms"], frame, descriptors)
    prediction = matrix @ model["coefficients"] + model["intercept"]
    if not np.isfinite(prediction).all():
        raise ValueError("Nonfinite SyMANTIC predictions")
    return prediction


def symbolic_equation(target, model, descriptors):
    import sympy as sp
    symbols = sp.symbols("X0:4")
    expression = sp.Float(format(model["intercept"], ".17g"))
    for coefficient, term in zip(model["coefficients"], model["terms"]):
        expression += sp.Float(format(coefficient, ".17g")) * sp.sympify(
            term, locals={str(symbol): symbol for symbol in symbols}
        )
    rhs = sp.latex(expression, symbol_names={
        symbol: energy_tex(name) for symbol, name in zip(symbols, descriptors)
    })
    return f"{energy_tex(target)} = {rhs}"


def run_symantic(argv=None):
    descriptors, targets, output = parse_arguments(
        "ALL_CHARGES_NONLINEAR_SYMANTIC_4D", ROOT / "NONLINEAR_FOUR_DESCRIPTORS", argv
    )
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
    frame = read_data()
    results = []
    for index, target in enumerate(targets, 1):
        observed = frame.loc[frame[target].notna()].copy()
        model = choose_symbolic(symbolic_candidates(observed, descriptors, target), target)
        fitted = predict_symbolic(model, observed, descriptors)
        heldout = np.full(len(observed), np.nan)
        for metal in sorted(observed.Metal.unique()):
            train = observed.loc[observed.Metal.ne(metal)]
            test = observed.loc[observed.Metal.eq(metal)]
            fold = choose_symbolic(symbolic_candidates(train, descriptors, target), target)
            heldout[observed.Metal.eq(metal).to_numpy()] = predict_symbolic(fold, test, descriptors)
        results.append({
            "target": target, "equation": symbolic_equation(target, model, descriptors),
            "complexity": model["complexity"],
            "actual": observed[target].to_numpy(float),
            "fitted": fitted, "heldout": heldout,
        })
        print(f"[{index}/{len(targets)}] {target} complete", flush=True)
    return write_results("ALL_CHARGES_NONLINEAR_SYMANTIC_4D", descriptors, results, output)
