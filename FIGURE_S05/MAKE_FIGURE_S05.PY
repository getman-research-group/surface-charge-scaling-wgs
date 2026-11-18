"""Render SI Figure S05 from the repository geometry measurements at 1200 DPI."""

from pathlib import Path
from runpy import run_path

root = Path(__file__).resolve().parent.parent
module = run_path(str(root / "SI_GEOMETRY_COMMON.PY"))
number = 5
spec = next(item for item in module["SPECS"] if item.number == number)
output = root / "FIGURE_S05" / "FIGURE_S05.PNG"

if __name__ == "__main__":
    print(module["draw_one"](spec, module["read_measurements"](), output))
