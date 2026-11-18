"""Fit one neutral-surface BEP relation and report its paper metrics.

Run `python3 BEP_NEUTRAL/MAKE_BEP_NEUTRAL.PY` from the repository root.
The script prompts for a TS target; use `--target ALL` for all four.
"""

from importlib.machinery import SourceFileLoader
from importlib.util import module_from_spec, spec_from_loader
from pathlib import Path
import sys


common = Path(__file__).resolve().parents[1] / "BEP_COMMON.PY"
spec = spec_from_loader("BEP_COMMON", SourceFileLoader("BEP_COMMON", str(common)))
module = module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)


if __name__ == "__main__":
    module.run_bep(True, Path(__file__).resolve().parent)
