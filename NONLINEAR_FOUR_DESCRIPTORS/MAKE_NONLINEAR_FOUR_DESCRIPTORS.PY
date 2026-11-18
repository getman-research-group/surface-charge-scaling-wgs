"""Fit full-dataset four-input nonlinear SRs with SyMANTIC and nested CV.

Run without arguments to select four inputs, or use
    python NONLINEAR_FOUR_DESCRIPTORS/MAKE_NONLINEAR_FOUR_DESCRIPTORS.PY \
        --descriptors 'CH*' 'HCO*' 'O*' 'H2O*'
The search uses +, -, *, /, up to three additive terms and complexity <= 15.
"""

from importlib.machinery import SourceFileLoader
from importlib.util import module_from_spec, spec_from_loader
from pathlib import Path
import sys

source = Path(__file__).resolve().parents[1] / "FOUR_DESCRIPTOR_COMMON.PY"
spec = spec_from_loader("FOUR_DESCRIPTOR_COMMON", SourceFileLoader("FOUR_DESCRIPTOR_COMMON", str(source)))
module = module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)

if __name__ == "__main__":
    module.run_symantic()
