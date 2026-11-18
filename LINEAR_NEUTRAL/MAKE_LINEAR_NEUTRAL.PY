"""Fit neutral-only two-descriptor linear scaling relations with CATMAP.

Run `python MAKE_LINEAR_NEUTRAL.PY` and enter two adsorption species when
prompted, or supply `--descriptors "CH*" "HCO*"` on the command line.
The script writes all 11 equations and paper metrics to one CSV file here.
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
    module.run_catmap(True, "NEUTRAL_LINEAR_CATMAP", Path(__file__).resolve().parent)
