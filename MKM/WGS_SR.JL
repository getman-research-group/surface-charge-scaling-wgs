function wgs_scaled_energies(descriptors; theta_CO=0.0)
    E_CH  = descriptors.CH
    E_O   = descriptors.O
    E_H2O = descriptors.H2O
    E_HCO = descriptors.HCO

    T = typeof(E_CH + theta_CO)
    E = Dict{Symbol, T}()

    E[:CH]  = E_CH
    E[:O]   = E_O
    E[:H2O] = E_H2O
    E[:HCO] = E_HCO

    E[:H] = 0.1254 * E_CH + 0.1232 * E_O + 0.5951
    E[:OH] = 0.7759 * E_H2O + 0.1333 * (E_CH * E_HCO) + 0.7556 * (E_HCO + E_O) - 1.4579
    E[:CO] = 16.5912 * E_HCO + 1.2231 * (E_CH * E_HCO) + 54.9842 * (E_HCO / E_CH) - 1.4895
    CO_exp_arg = 4.79 * theta_CO
    E[:CO_ads] = E[:CO] + 0.0305 * exp(CO_exp_arg) + 0.031135 * theta_CO * exp(CO_exp_arg)
    E[:CO2] = 0.1553 * (-E_H2O * E_HCO) + 0.1063
    E[:COOH] = 0.6069 * (E_CH + E_H2O) - 0.0277 * (E_O / E_H2O) - 0.2019 * E_CH + 0.6920

    E[:H_cov] = E[:H] + 0.4323 / (1.0 + exp(5.3619 * (0.5189 - theta_CO)))
    E[:O_cov] = E[:O] + 1.8276 / (1.0 + exp(6.6181 * (0.4745 - theta_CO)))
    E[:OH_cov] = E[:OH] + 0.0239 * exp(3.2504 * theta_CO)
    E[:COOH_cov] = E[:COOH] + 0.0712 * exp(4.4703 * theta_CO)
    E[:H_ads] = E[:H_cov]
    E[:CO_cov] = E[:CO_ads]

    E[:TS_HO_H] = -0.6481 * E_H2O + 0.3387 * E_HCO + 1.4815
    E[:TS_H_CO2] = -1.0278 * E_H2O + 1.3080 * E_HCO - 1.1667 * (E_O / E_HCO) + 4.2651
    E[:TS_O_CO] = 0.3269 * E_O + 2.1770 * E_CH - 2.6126 * E_HCO + 9.9781
    E[:TS_OC_OH] = 0.2491 * E_CH + 0.1719 * E_O - 0.4775 * E_HCO + 1.5356

    return E
end
