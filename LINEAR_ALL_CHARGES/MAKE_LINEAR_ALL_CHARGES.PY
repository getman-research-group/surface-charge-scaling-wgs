"""Fit all-charge two-descriptor linear scaling relations with CATMAP.

Run `python MAKE_LINEAR_ALL_CHARGES.PY` and enter two adsorption species
when prompted, or supply `--descriptors "CO2*" "COOH*"`. All charge states
of a metal are held out together during nested LOMO validation.
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
    module.run_catmap(False, "ALL_CHARGES_LINEAR_CATMAP", Path(__file__).resolve().parent)
