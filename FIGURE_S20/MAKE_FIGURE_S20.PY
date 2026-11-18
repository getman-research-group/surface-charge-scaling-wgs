"""Rebuild SI Figure S20 at 1200 DPI from the full nonlinear winner PKL."""

from pathlib import Path
from runpy import run_path

COMMON = Path(__file__).resolve().parent.parent / "TWO_DESCRIPTOR_PARITY.PY"
make_figure = run_path(str(COMMON))["make_figure"]

if __name__ == "__main__":
    print(make_figure("S20"))
