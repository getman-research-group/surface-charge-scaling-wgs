"""Full-dataset nonlinear four-descriptor SRs selected in the manuscript.

CH, HCO, H2O, and O are their DFT adsorption energies in eV. The equation
coefficients retain the four-decimal precision used to evaluate the SI fits.
The coverage corrections follow the values printed in the SI.
"""
function wgs_scaled_energies(descriptors; theta_CO=0.0)
    E_CH = descriptors.CH
    E_HCO = descriptors.HCO
    E_H2O = descriptors.H2O
    E_O = descriptors.O

    T = typeof(E_CH + theta_CO)
    E = Dict{Symbol, T}()

    E[:CH] = E_CH
    E[:HCO] = E_HCO
    E[:H2O] = E_H2O
    E[:O] = E_O

    E[:H] = 0.2597 * E_CH + 0.7698 * E_O + 0.1035 * E_CH * E_O + 1.3831
    E[:OH] = -0.2902 * E_CH + 0.7747 * E_O + 0.7281 * E_H2O - 3.0426
    E[:CO] = 0.07760 * E_CH / E_O - 0.3185 * E_H2O + 0.8779 * E_HCO - 0.2454
    E[:CO2] = -0.1688 * E_H2O * E_HCO - 0.0575 * E_CH - 0.2232
    E[:COOH] = 0.4432 * E_CH + 0.5026 * E_H2O + 0.7395

    E[:TS_O_CO] = 1.5452 * E_HCO * E_O + 3.4466 * E_O + 1.8019 * E_CH + 12.8399
    E[:TS_HO_H] = -0.0905 * E_CH / E_O - 0.6532 * E_H2O - 0.2211 * E_HCO * E_O + 1.999
    E[:TS_H_CO2] = -3.1810 * E_O / E_HCO - 1.3023 * (E_H2O + E_O) + 2.2359 * E_HCO + 5.6487
    E[:TS_OC_OH] = -1.3585 * E_O / E_HCO - 0.5731 * E_O + 0.2562 * E_CH + 2.4863

    E[:CO_ads] = E[:CO] + (0.0065 + 0.031135 * theta_CO) * exp(4.79 * theta_CO)
    E[:H_cov] = E[:H] + 0.4323 / (1.0 + exp(5.3619 * (0.5189 - theta_CO)))
    E[:O_cov] = E[:O] + 1.8176 / (1.0 + exp(6.6181 * (0.4745 - theta_CO)))
    E[:OH_cov] = E[:OH] + 0.0239 * exp(3.2504 * theta_CO)
    E[:COOH_cov] = E[:COOH] + 0.0712 * exp(4.4703 * theta_CO)
    E[:H_ads] = E[:H_cov]
    E[:CO_cov] = E[:CO_ads]

    return E
end
