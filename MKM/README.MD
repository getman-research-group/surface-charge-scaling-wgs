# Pt WGS microkinetic model

`RUN_NEUTRAL.JL` runs the neutral-charge Pt case at 548 K using the selected full-dataset nonlinear SRs for the adsorption energies of CH*, HCO*, H2O*, and O* as inputs. `WGS_SR.JL` contains the nine selected equations with the coefficient precision stored in `NONLINEAR_FOUR_DESCRIPTORS/RESULTS_CH_HCO_O_H2O.PKL`. It uses the CO and O coverage-correction constants printed in the SI: 0.0065 and 1.8176. The original Julia kinetic implementation's zero floor on negative reverse barriers is retained in `MAKEKS.JL`.

The kinetic arrays retain the original 14-slot numbering: slot 13 is inactive, and slot 14 corresponds to manuscript reaction R13 (HCO* dissociation).

With Julia 1.11.3, run from the repository root:

```bash
julia --project=MKM -e 'using Pkg; Pkg.instantiate()'
julia --project=MKM MKM/RUN_NEUTRAL.JL
```

The run writes `NEUTRAL_RESULT.CSV`. With these inputs, the H2 TOF is **0.0028807436 min⁻¹** (4.80123933 × 10⁻⁵ s⁻¹). This differs from the manuscript's current 7.71 min⁻¹ entry. The earlier Julia run used different four-descriptor SR equations, different CO and O coverage constants, and multiplied the s⁻¹ TOF by 600 when labeling it min⁻¹. A s⁻¹ to min⁻¹ conversion uses 60.

The figure and SI should be reconciled with the model implementation before treating their reported TOF as reproduced by this code.
