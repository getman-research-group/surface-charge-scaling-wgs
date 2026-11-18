"""Render SI Figure S23 at 1200 DPI from the nonlinear four-descriptor PKL."""

from pathlib import Path
from runpy import run_path

source = Path(__file__).resolve().parent.parent / "FOUR_DESCRIPTOR_PARITY.PY"
plot = run_path(str(source))["plot_s23"]

if __name__ == "__main__":
    print(plot())
