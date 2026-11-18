# Two-descriptor scaling-relation scripts

The four folders each contain one runnable model. Every script reads the
full-precision `ADSORPTION_TS_ENERGIES.PKL` in the repository root. Run it
without arguments to be prompted for two adsorbate labels, or use
`--descriptors` as in the examples below. Quote labels containing `*` or a
space.

| Folder | Model | Example command |
| --- | --- | --- |
| `LINEAR_NEUTRAL` | CATMAP, neutral-only | `python3 LINEAR_NEUTRAL/MAKE_LINEAR_NEUTRAL.PY --descriptors 'CH*' 'HCO*'` |
| `LINEAR_ALL_CHARGES` | CATMAP, all charges | `python3 LINEAR_ALL_CHARGES/MAKE_LINEAR_ALL_CHARGES.PY --descriptors 'CO2*' 'COOH*'` |
| `NONLINEAR_NEUTRAL` | SyMANTIC, neutral-only | `python3 NONLINEAR_NEUTRAL/MAKE_NONLINEAR_NEUTRAL.PY --descriptors 'OH*' 'COOH*'` |
| `NONLINEAR_ALL_CHARGES` | SyMANTIC, all charges | `python3 NONLINEAR_ALL_CHARGES/MAKE_NONLINEAR_ALL_CHARGES.PY --descriptors 'CO* FCC' 'H2O*'` |

Run commands from the repository root. The script writes `RESULTS_<PAIR>.CSV`
to its own folder. Use `--output PATH` to choose another CSV path. The
optional `--target 'H*'` runs one target for a diagnostic check. A complete
run reports 11 equations: seven other adsorption energies and four activation
energies. `CO* ATOP` and transition-state energies cannot be descriptors.
Charge is not an equation input.

The CSV contains each full-precision fitted equation, its fit coefficient of
determination (`R2`), fit `RMSE_EV`, and nested `LOMO_RMSE_EV`. Three final rows
give pooled fit and nested LOMO RMSE across all predicted species, adsorption
energies, and transition states. All reported errors are in eV. The equations
use LaTeX adsorption-energy and activation-energy notation. The metrics are
calculated from unrounded equations. The source pickle is never modified.

For each nested LOMO fold, all observations for one metal are held out.
Equation selection and coefficient fitting are repeated using the other
metals. The descriptor pair supplied by the user remains fixed. CATMAP tests
each input alone and both together. SyMANTIC tests `+`, `-`, `*`, and `/`
expressions and applies the original model-specific complexity and selection
rules. The shared implementation is `TWO_DESCRIPTOR_COMMON.PY`.

Python dependencies are `numpy`, `pandas`, `sympy`, and either
`python-catmap==0.3.1` or `symantic==1.1.0.6` with PyTorch. On the OSC
environment used for this paper, CATMAP was loaded with
`PYTHONPATH=/fs/ess/PAS2536/sifat149/SR/.catmap_install`, while SyMANTIC
was loaded with `module load pytorch/2.7.1` followed by
`PYTHONPATH=/fs/ess/PAS2536/sifat149/SR/.symantic_install`. A separate
environment with those packages installed also works; `PYTHONPATH` is only
needed for the paper's local package installations.
