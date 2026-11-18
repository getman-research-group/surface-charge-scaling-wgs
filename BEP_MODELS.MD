# BEP relation scripts

The two scripts fit a Brønsted–Evans–Polanyi (BEP) line for one of the four
transition-state targets in `ADSORPTION_TS_ENERGIES.PKL`:

| Dataset | Run from the repository root |
| --- | --- |
| Neutral surfaces | `python3 BEP_NEUTRAL/MAKE_BEP_NEUTRAL.PY` |
| All charges | `python3 BEP_ALL_CHARGES/MAKE_BEP_ALL_CHARGES.PY` |

Each script prompts for `HO--H`, `O--CO`, `H--CO2`, or `OC--OH`. Enter `ALL` to
fit all four targets. For noninteractive use, add `--target O--CO` or
`--target ALL`. The CSV is written to the script's folder as
`RESULTS_<TARGET>.CSV`; `--output PATH` overrides the destination. Each row
contains the fitted equation, fit $R^2$, fit RMSE, nested leave-one-metal-out
(LOMO) RMSE, and observation count. RMSEs are in eV. The scripts require
Python, NumPy, and pandas. Their shared formulas are in `BEP_COMMON.PY`.

For each forward reaction as written in the manuscript, the predictor is

$$
\Delta E_i^{\mathrm{rxn}}=
\sum_{p}\left(\Delta E_p^{\mathrm{ad}}+G_p\right)
-\sum_{r}\left(\Delta E_r^{\mathrm{ad}}+G_r\right),
$$

where $p$ and $r$ are the adsorbed products and reactants. The $G$ values are
the neutral gas-phase electronic-energy references used to define the
adsorption energies. For H* and O*, the references are respectively half the
H₂ and O₂ electronic energies. Their numerical values are recorded in
`BEP_COMMON.PY` and come from the original BEP analysis. A least-squares line
with intercept then predicts $\Delta E_i^{\mathrm{act}}$ from
$\Delta E_i^{\mathrm{rxn}}$. Each LOMO fold fits that same linear form using
the other three metals. For all-charge fitting, it omits every charge state of
the held-out metal together. No missing transition-state value is imputed.

**HO–H sign correction.** The original `ts_reaction_energies/calculate.py`
computed its HO–H reaction energy as reactants minus products, opposite to the
forward water-formation reaction printed in the manuscript. These scripts use
products minus reactants. The corrected HO–H BEP slope therefore has the
opposite sign from the old SI equation, while its fitted predictions, $R^2$,
fit RMSE, and LOMO RMSE are unchanged. The other three equations retain their
original slope signs.
