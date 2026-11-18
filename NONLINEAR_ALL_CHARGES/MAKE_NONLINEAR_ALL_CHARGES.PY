"""Fit all-charge two-descriptor symbolic scaling relations with SyMANTIC.

Run `python MAKE_NONLINEAR_ALL_CHARGES.PY` and enter two adsorption species
when prompted, or supply `--descriptors "CO* FCC" "H2O*"`. Nested LOMO
repeats equation search with every charge of the held-out metal omitted.
"""

from pathlib import Path
from importlib.machinery import SourceFileLoader
from importlib.util import module_from_spec, spec_from_loader
import sys

common = Path(__file__).resolve().parents[1] / "TWO_DESCRIPTOR_COMMON.PY"
spec = spec_from_loader("TWO_DESCRIPTOR_COMMON", SourceFileLoader("TWO_DESCRIPTOR_COMMON", str(common)))
module = module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)


if __name__ == "__main__":
    module.run_symantic(False, "ALL_CHARGES_NONLINEAR_SYMANTIC", Path(__file__).resolve().parent)
