# Energy dataset data dictionary

This document describes `adsorption_TS_energies.pkl` **as stored**. The file is a
pickled pandas DataFrame with 25 rows and 17 columns: one row for each of five
metal surfaces at each of five charge states. Its energy values are DFT
**electronic-energy differences**, in eV, rather than absolute total energies or
Gibbs free energies. All numeric columns are stored as `float64`; no rounding is
applied to the original pickle.

## Row identifiers and charge

| Column | Meaning |
| --- | --- |
| `Metal` | Element symbol: `Pt`, `Cu`, `Ni`, `Au`, or `Ru`. |
| `Metal_name` | Full element name: Platinum, Copper, Nickel, Gold, or Ruthenium, respectively. |
| `Charge` | Nominal charge coordinate, in $h^+$ per surface metal atom. A negative value denotes electron addition relative to the neutral slab; a positive value denotes electron removal. |

The surfaces are Pt(111), Cu(111), Ni(111), Au(111), and Ru(0001). Each has one
row at every stored charge. The stored decimal charge coordinates are rounded
labels for the charge states used in the manuscript; they are **not** numerically
equal to the exact fractions:

| Stored `Charge` | Charge state ($h^+$ per surface metal atom) |
| ---: | ---: |
| `-0.11` | $-2/18$ |
| `-0.05` | $-1/18$ |
| `0.00` | $0$ |
| `+0.05` | $+1/18$ |
| `+0.11` | $+2/18$ |

Here $h^+$ is a unit for electron deficiency. The charge label does not specify
where a hole is localized in the slab.

## Adsorption-energy columns

An asterisk (`*`) denotes a surface-bound species; a lone `*` denotes a vacant
site. Each adsorption-energy column contains

$$
\Delta E^{\mathrm{ad}}_{s*}(q)
=E_{s*/\mathrm{surf}}(q)-E_{\mathrm{surf}}(q)-E_{s(g)},
$$

where the adsorbate/slab and clean-slab electronic energies are calculated at
the **same charge** $q$. The gas-phase reference is calculated in a neutral
cell. For `H*`, replace $E_{s(g)}$ with $\tfrac12 E_{\mathrm{H_2}(g)}$; for
`O*`, use $\tfrac12 E_{\mathrm{O_2}(g)}$. For the other species, the reference
is the corresponding neutral gas-phase species. Negative values indicate
adsorption favorable relative to the stated reference. The pickle contains the
resulting differences, not the underlying slab or gas-phase total energies.

| Pickle column | Adsorbed species/site | Gas-phase reference | Available rows |
| --- | --- | --- | ---: |
| `H*` | Atomic H | $\tfrac12\mathrm{H_2}(g)$ | 25 |
| `O*` | Atomic O | $\tfrac12\mathrm{O_2}(g)$ | 25 |
| `OH*` | OH | $\mathrm{OH}(g)$ | 25 |
| `CH*` | CH | $\mathrm{CH}(g)$ | 25 |
| `CO* FCC` | CO at an FCC hollow site | $\mathrm{CO}(g)$ | 25 |
| `CO* ATOP` | CO at an atop site | $\mathrm{CO}(g)$ | 5, all on Pt(111) |
| `CO2*` | CO$_2$ | $\mathrm{CO_2}(g)$ | 25 |
| `H2O*` | H$_2$O | $\mathrm{H_2O}(g)$ | 25 |
| `HCO*` | HCO | $\mathrm{HCO}(g)$ | 25 |
| `COOH*` | COOH | $\mathrm{COOH}(g)$ | 25 |

`FCC` means face-centered-cubic hollow. `CO* ATOP` is retained in the source
dataset but **excluded from the scaling-relation fits and descriptor search**:
only Pt has atop-CO values. The nine other adsorption-energy columns supply
the candidate adsorption-energy descriptors and adsorption-energy targets.

## Transition-state columns

The double hyphen in each column name is a transition-state label, not a
subtraction operator. These columns contain **forward activation energies**
$\Delta E^{\mathrm{act}}_i(q)$ for the reactions as written below, not
absolute TS electronic energies:

$$
\Delta E^{\mathrm{act}}_i(q)
=E_{\mathrm{TS},i}(q)-E_{\mathrm{reactants},i}(q).
$$

Reactant and TS slab energies are evaluated at the same charge. For a reaction
with two separately adsorbed reactants, the reactant energy is their sum and a
clean-slab energy is added to the TS side where needed to balance the number of
slabs. These conventions follow the manuscript's energy definitions.

| Pickle column | Forward elementary reaction | Available rows |
| --- | --- | ---: |
| `HO--H` | OH* + H* $\rightarrow$ H$_2$O* + * | 20 |
| `O--CO` | CO* at FCC + O* $\rightarrow$ CO$_2$* + * | 20 |
| `H--CO2` | COOH* + * $\rightarrow$ CO$_2$* + H* | 14 |
| `OC--OH` | OH* + CO* at FCC $\rightarrow$ COOH* + * | 20 |

## Unavailable values

Missing values are stored as `NaN`, not zero. They should be excluded from
fits and metrics unless a separate, documented estimation method is applied.

| Column | Unavailable metal/charge combinations |
| --- | --- |
| `CO* ATOP` | All five charges on Cu, Ni, Au, and Ru (20 rows). |
| `HO--H` | All five charges on Ru (5 rows). |
| `O--CO` | All five charges on Ru (5 rows). |
| `OC--OH` | All five charges on Ru (5 rows). |
| `H--CO2` | All five charges on Ru; `+0.11` on Pt, Cu, Ni, and Au; `-0.05` on Ni and Au (11 rows total). |

The four TS searches were not included on Ru. The pickle itself does not record
a reason for each additional unavailable `H--CO2` value; it records only the
absence of a numerical barrier. The nine analyzed adsorption-energy columns
contain 225 values, and the four TS columns contain 74 values. Thus the
scaling-relation dataset has **299 available energy values**. The five Pt atop-CO
values bring the pickled file's total to 304 available energy values.

The original pickle is 4,026 bytes and has SHA-256 checksum
`885475f77c9c118f29e2dcf8c009077cfb75282978cfed5191c2bbb13c370d17`.
