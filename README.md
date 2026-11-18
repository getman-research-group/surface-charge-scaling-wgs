# Surface Charge Reshapes Scaling Relationships

Code and processed data for the manuscript **“Surface Charge Reshapes Scaling Relationships: From Linear to Nonlinear, Multi-Descriptor Models”** by Venkata Rohit Punyapu, Sifat Hossain, Matthew Neurock, Lars Grabow, and Rachel B. Getman.

This repository contains the DFT-derived adsorption and activation energies, saved scaling-relation (SR) and Brønsted–Evans–Polanyi (BEP) results, analysis scripts, the Pt water-gas-shift (WGS) microkinetic model, and figure-generation scripts. The saved result files and plot data allow the manuscript and supporting-information figures to be drawn without rerunning the equation searches.

## Repository map

| Location | Contents |
| --- | --- |
| [ADSORPTION_TS_ENERGIES.PKL](ADSORPTION_TS_ENERGIES.PKL), [DATA_DICTIONARY.MD](DATA_DICTIONARY.MD) | Original full-precision DFT-derived energy differences, units, references, and missing-value conventions. |
| `LINEAR_NEUTRAL/`, `LINEAR_ALL_CHARGES/` | CATMAP linear fits for all 36 two-descriptor pairs in each dataset. |
| `NONLINEAR_NEUTRAL/`, `NONLINEAR_ALL_CHARGES/` | SyMANTIC fits for all 36 two-descriptor pairs in each dataset. |
| `LINEAR_FOUR_DESCRIPTORS/`, `NONLINEAR_FOUR_DESCRIPTORS/` | Linear and nonlinear fits for each of the ten tested four-descriptor sets. |
| `BEP_NEUTRAL/`, `BEP_ALL_CHARGES/` | Neutral-only and full-dataset BEP equations and metrics for the four studied transition states. |
| `FIGURE_1/`–`FIGURE_6/` | Main-text figure scripts and 1200-DPI PNGs. |
| `FIGURE_S01/`–`FIGURE_S23/` | Supporting-information figure scripts and 1200-DPI PNGs. |
| `PT_DESCRIPTOR_INTERPOLATION/` | Pt adsorption-energy interpolation analysis and saved curves. |
| [GEOMETRY_MEASUREMENTS.CSV](GEOMETRY_MEASUREMENTS.CSV) | Measurements used in the geometry schematics. |
| `MKM/` | Julia WGS microkinetic model, neutral-case runner, and Julia project files. |

The shared fitting and parity-plot implementations are in the repository root. [TWO_DESCRIPTOR_MODELS.MD](TWO_DESCRIPTOR_MODELS.MD), [BEP_MODELS.MD](BEP_MODELS.MD), and [MKM/README.MD](MKM/README.MD) give model-specific commands and conventions.

## Data and model results

`ADSORPTION_TS_ENERGIES.PKL` is a pandas DataFrame with one row for each of five metal surfaces and five nominal charge states. It contains adsorption energies for nine analyzed adsorbates and activation energies for four transition states, in eV. The separate atop-CO column is not used in the SR search. The pickle stores electronic-energy differences, rather than the underlying DFT total energies. See the [data dictionary](DATA_DICTIONARY.MD) before using the values.

Each `RESULTS_*.PKL` file contains the equations and manuscript performance metrics for one descriptor set or BEP dataset. The four two-descriptor families contain 36 files each; the two four-descriptor families contain ten files each. The result metrics include fit coefficient of determination, fit RMSE, and leave-one-metal-out (LOMO) RMSE. For nonlinear models, the nested LOMO procedure repeats equation discovery and selection in every metal holdout fold. The result pickles are the archived values used for the manuscript tables and comparison figures.

The adsorption energies are the model inputs. A species listed as an input is not scored as a prediction target for that descriptor set. Charge is not an explicit input in any SR equation.

## Recreate the figures

Run plotters from the repository root with Python 3. For example:

```bash
python FIGURE_2/MAKE_FIGURE_2.PY
python FIGURE_5/MAKE_FIGURE_5.PY
python FIGURE_6/MAKE_FIGURE_6.PY
python FIGURE_S16/MAKE_FIGURE_S16.PY
```

Every figure folder follows the same naming pattern: `FIGURE_<N>/MAKE_FIGURE_<N>.PY` for main-text figures and `FIGURE_S<NN>/MAKE_FIGURE_S<NN>.PY` for SI figures. Each plotter writes its PNG beside the script. Figure 6 reads the 23-charge plot dataset in [FIGURE_6/PLOT_DATA.PKL](FIGURE_6/PLOT_DATA.PKL). The other figure scripts read the energy data, geometry measurements, interpolation results, and/or saved fit results as appropriate.

The figure PNGs are included so they can be inspected without running Python. Geometry drawings are schematic representations of the measured distances and angles; the plotted geometry values are documented in the scripts and `GEOMETRY_MEASUREMENTS.CSV`.

## Run a new fit

From the repository root, supply adsorption-energy inputs with `--descriptors`, or omit them to be prompted. For example:

```bash
python LINEAR_NEUTRAL/MAKE_LINEAR_NEUTRAL.PY --descriptors 'CH*' 'HCO*'
python NONLINEAR_ALL_CHARGES/MAKE_NONLINEAR_ALL_CHARGES.PY --descriptors 'CO* FCC' 'H2O*'
python LINEAR_FOUR_DESCRIPTORS/MAKE_LINEAR_FOUR_DESCRIPTORS.PY --descriptors 'CH*' 'HCO*' 'O*' 'H2O*'
python NONLINEAR_FOUR_DESCRIPTORS/MAKE_NONLINEAR_FOUR_DESCRIPTORS.PY --descriptors 'CH*' 'HCO*' 'O*' 'H2O*'
python BEP_NEUTRAL/MAKE_BEP_NEUTRAL.PY --target ALL
python BEP_ALL_CHARGES/MAKE_BEP_ALL_CHARGES.PY --target ALL
```

These commands write new CSV results; they do not overwrite the archived `RESULTS_*.PKL` files. The fitting implementations, equation-selection rules, and validation procedures are described in [TWO_DESCRIPTOR_MODELS.MD](TWO_DESCRIPTOR_MODELS.MD), [BEP_MODELS.MD](BEP_MODELS.MD), and the manuscript methods.

## Microkinetic model

The Julia model and its pinned package manifest are in `MKM/`. To run the neutral Pt case documented in [MKM/README.MD](MKM/README.MD):

```bash
julia --project=MKM -e 'using Pkg; Pkg.instantiate()'
julia --project=MKM MKM/RUN_NEUTRAL.JL
```

This runner writes a local `MKM/NEUTRAL_RESULT.CSV`. Figure 6 is drawn from the saved 23-charge comparison values in `FIGURE_6/PLOT_DATA.PKL`.

## Software

The figure and analysis scripts were checked with Python 3.10.13, NumPy 2.2.6, pandas 2.3.1, Matplotlib 3.10.3, and SymPy 1.14.0. The linear fitting scripts use `python-catmap==0.3.1`. The nonlinear fitting scripts use `symantic==1.1.0.6` and PyTorch; see [TWO_DESCRIPTOR_MODELS.MD](TWO_DESCRIPTOR_MODELS.MD) for the OSC module and local package paths used for the study. The MKM uses Julia 1.11.3; `MKM/Project.toml` and `MKM/Manifest.toml` specify its packages. Plotting archived results does not require CATMAP, SyMANTIC, or Julia.

## Citation and license

Please cite the associated manuscript when using this repository. Its title and authors are listed in [CITATION.cff](CITATION.cff); journal details and DOI can be added when available. Questions about the manuscript or data can be sent to the corresponding author, Rachel B. Getman, at getman.11@osu.edu.

Code is licensed under MIT. Data, figures, and documentation are licensed under CC BY 4.0. See [LICENSE.md](LICENSE.md) for the scope of each license.
