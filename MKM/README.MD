# Pt WGS microkinetic model

`RUN_NEUTRAL.JL` runs the neutral-charge Pt case at 548 K with the four adsorption-energy inputs for CH*, HCO*, H2O*, and O*. The remaining adsorbate and transition-state energies are computed by `WGS_SR.JL`. The kinetic arrays retain the original 14-slot numbering: slot 13 is inactive, and slot 14 corresponds to manuscript reaction R13.

With Julia 1.11.3, run from the repository root:

```bash
julia --project=MKM -e 'using Pkg; Pkg.instantiate()'
julia --project=MKM MKM/RUN_NEUTRAL.JL
```

The run writes `NEUTRAL_RESULT.CSV`. Its `TOF_H2_REPORTED` field is 7.70848717 at neutral charge. The CSV also includes the computed rate in s⁻¹ and its standard conversion to min⁻¹.
