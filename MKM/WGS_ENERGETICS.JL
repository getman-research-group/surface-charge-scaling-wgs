include("WGS_SR.JL")

const WGS_KJ_PER_MOL_PER_EV = 96.48533212331002
const WGS_OMEGA_RXN5 = 0.8601
const WGS_OMEGA_RXN8 = 0.0
const WGS_OMEGA_RXN9 = 0.0107
const WGS_OMEGA_RXN10 = 0.1425
const WGS_EF0_RXN6 = 1.09
const WGS_EF0_RXN7 = 0.0
const WGS_EF0_RXN11 = 0.35
const WGS_EF0_RXN12 = 0.12
const WGS_EF0_RXN14 = 0.23
const WGS_OMEGA_RXN6 = 0.5
const WGS_OMEGA_RXN7 = 0.5
const WGS_OMEGA_RXN11 = 0.5
const WGS_OMEGA_RXN12 = 0.1696
const WGS_OMEGA_RXN14 = 0.5
const WGS_STICKING_H2O = 1.0
const WGS_STICKING_CO = 1.0
const WGS_STICKING_H2 = 1.0
const WGS_STICKING_CO2 = 1.0

function descriptor_value(descriptors, name::Symbol, default)
    return hasproperty(descriptors, name) ? getproperty(descriptors, name) : default
end

function shomate_enthalpy_increment_kjmol(T, coeffs)
    t = T / 1000.0
    A, B, C, D, E, F, G, H = coeffs
    return A * t + B * t^2 / 2.0 + C * t^3 / 3.0 + D * t^4 / 4.0 - E / t + F - H
end

function make_adsorption_sticking_constants(MW, erxn, eact, sticking, consts)
    m = MW / (1000.0 * consts.Na)
    impingement_per_bar_per_site = 1.0e5 / (consts.No * sqrt(2.0 * pi * m * consts.kb * consts.Te))

    kf = sticking * impingement_per_bar_per_site * exp(-eact / (consts.kbeV * consts.Te))
    keq = exp(-erxn / (consts.kbeV * consts.Te))
    kr = kf / keq

    return kf, kr, keq
end

function nasa7_enthalpy_ev(T, low_coeffs, high_coeffs; T_mid=1000.0)
    a = T <= T_mid ? low_coeffs : high_coeffs
    a1, a2, a3, a4, a5, a6, a7 = a
    h_over_RT = a1 + a2 * T / 2.0 + a3 * T^2 / 3.0 + a4 * T^3 / 4.0 + a5 * T^4 / 5.0 + a6 / T
    return (8.314462618 * T * h_over_RT / 1000.0) / WGS_KJ_PER_MOL_PER_EV
end

function gas_enthalpy_ev(species::Symbol, T)
    coeffs = if species == :H2O
        if T <= 1700.0
            # NIST gas-phase Shomate coefficients are tabulated from 500 K;
            # this lowest gas range is used slightly below 500 K for 473 K runs.
            (30.09200, 6.832514, 6.793435, -2.534480, 0.082139, -250.8810, 223.3967, -241.8264)
        else
            (41.96426, 8.622053, -1.499780, 0.098119, -11.15764, -272.1797, 219.7809, -241.8264)
        end
    elseif species == :OH
        if T <= 1300.0
            (32.27768, -11.36291, 13.60545, -3.846486, -0.001335, 29.75113, 225.5783, 38.98706)
        else
            (28.74701, 4.714489, -0.814725, 0.054748, -2.747829, 26.41439, 214.1166, 38.98706)
        end
    elseif species == :H
        (20.78603, 4.850638e-10, -1.582916e-10, 1.525102e-11, 3.196347e-11, 211.8020, 139.8711, 217.9994)
    elseif species == :H2
        if T <= 1000.0
            (33.066178, -11.363417, 11.432816, -2.772874, -0.158558, -9.980797, 172.707974, 0.0)
        elseif T <= 2500.0
            (18.563083, 12.257357, -2.859786, 0.268238, 1.977990, -1.147438, 156.288133, 0.0)
        else
            (43.413560, -4.293079, 1.272428, -0.096876, -20.533862, -38.515158, 162.081354, 0.0)
        end
    elseif species == :CO
        if T <= 1300.0
            (25.56759, 6.096130, 4.054656, -2.671301, 0.131021, -118.0089, 227.3665, -110.5271)
        else
            (35.15070, 1.300095, -0.205921, 0.013550, -3.282780, -127.8375, 231.7120, -110.5271)
        end
    elseif species == :CO2
        if T <= 1200.0
            (24.99735, 55.18696, -33.69137, 7.948387, -0.136638, -403.6075, 228.2431, -393.5224)
        else
            (58.16639, 2.720074, -0.492289, 0.038844, -6.447293, -425.9186, 263.6125, -393.5224)
        end
    elseif species == :HCO
        if T <= 1200.0
            (21.13803, 40.43610, -14.71337, 0.969010, 0.239639, 36.34712, 240.1695, 43.51402)
        else
            (52.79371, 2.666155, -0.392339, 0.023808, -7.457018, 11.37797, 267.2798, 43.51402)
        end
    else
        if species == :O
            return (249.18 + 2.5 * 8.314462618 * (T - 298.15) / 1000.0) / WGS_KJ_PER_MOL_PER_EV
        elseif species == :O2
            coeffs_O2 = if T <= 700.0
                (31.32234, -20.23531, 57.86644, -36.50624, -0.007374, -8.903471, 246.7945, 0.0)
            elseif T <= 2000.0
                (30.03235, 8.772972, -3.988133, 0.788313, -0.741599, -11.32468, 236.1663, 0.0)
            else
                (20.91111, 10.72071, -2.020498, 0.146449, 9.245722, 5.337651, 237.6185, 0.0)
            end
            return shomate_enthalpy_increment_kjmol(T, coeffs_O2) / WGS_KJ_PER_MOL_PER_EV
        elseif species == :COOH
            low = (
                2.92207919E+00,
                7.62453859E-03,
                3.29884437E-06,
                -1.07135205E-08,
                5.11587057E-12,
                -2.33512327E+04,
                1.12925886E+01,
            )
            high = (
                5.39206152E+00,
                4.11221455E-03,
                -1.48194900E-06,
                2.39875460E-10,
                -1.43903104E-14,
                -2.36480760E+04,
                -2.23529091E+00,
            )
            return nasa7_enthalpy_ev(T, low, high)
        end

        error("No gas enthalpy data implemented for species $species.")
    end

    Hf298_kjmol = coeffs[8]
    return (Hf298_kjmol + shomate_enthalpy_increment_kjmol(T, coeffs)) / WGS_KJ_PER_MOL_PER_EV
end

function adsorbate_gas_reference_enthalpy_ev(species::Symbol, T)
    if species == :H
        return 0.5 * gas_enthalpy_ev(:H2, T)
    elseif species == :O
        return 0.5 * gas_enthalpy_ev(:O2, T)
    end

    return gas_enthalpy_ev(species, T)
end

function adsorbate_enthalpy_ev(species::Symbol, E, consts, use_coverage_dependent_thermo)
    return adsorbate_gas_reference_enthalpy_ev(species, consts.Te) +
           surface_energy(E, species, use_coverage_dependent_thermo)
end

function wgs_reaction5_dft_energy(E, G; use_coverage_dependent_thermo=false)
    E_H = surface_energy(E, :H, use_coverage_dependent_thermo)
    E_OH = surface_energy(E, :OH, use_coverage_dependent_thermo)
    return E_OH + E_H - E[:H2O] + G[:OH] + adsorbate_gas_reference_energy(:H, G) - G[:H2O]
end

function adsorbate_gas_reference_energy(species::Symbol, G)
    if species == :H
        if !isfinite(G[:H2])
            error("H* is referenced to 1/2 H2(g), so descriptor G_H2 is required.")
        end
        return 0.5 * G[:H2]
    elseif species == :O
        if !isfinite(G[:O2])
            error("O* is referenced to 1/2 O2(g), so descriptor G_O2 is required.")
        end
        return 0.5 * G[:O2]
    end

    return G[species]
end

function wgs_reaction_dft_energy(i, E, G; use_coverage_dependent_thermo=false)
    E_H = surface_energy(E, :H, use_coverage_dependent_thermo)
    E_O = surface_energy(E, :O, use_coverage_dependent_thermo)
    E_OH = surface_energy(E, :OH, use_coverage_dependent_thermo)
    E_CO = surface_energy(E, :CO, use_coverage_dependent_thermo)
    E_COOH = surface_energy(E, :COOH, use_coverage_dependent_thermo)
    E_HCO = E[:HCO]

    if i == 5
        return E_OH + E_H - E[:H2O] +
               adsorbate_gas_reference_energy(:OH, G) +
               adsorbate_gas_reference_energy(:H, G) -
               adsorbate_gas_reference_energy(:H2O, G)
    elseif i == 6
        return E_O + E_H - E_OH +
               adsorbate_gas_reference_energy(:O, G) +
               adsorbate_gas_reference_energy(:H, G) -
               adsorbate_gas_reference_energy(:OH, G)
    elseif i == 7
        return E[:H2O] + E_O - 2.0 * E_OH +
               adsorbate_gas_reference_energy(:H2O, G) +
               adsorbate_gas_reference_energy(:O, G) -
               2.0 * adsorbate_gas_reference_energy(:OH, G)
    elseif i == 8
        return E[:CO2] - E_CO - E_O +
               adsorbate_gas_reference_energy(:CO2, G) -
               adsorbate_gas_reference_energy(:CO, G) -
               adsorbate_gas_reference_energy(:O, G)
    elseif i == 9
        return E_COOH - E_CO - E_OH +
               adsorbate_gas_reference_energy(:COOH, G) -
               adsorbate_gas_reference_energy(:CO, G) -
               adsorbate_gas_reference_energy(:OH, G)
    elseif i == 10
        return E[:CO2] + E_H - E_COOH +
               adsorbate_gas_reference_energy(:CO2, G) +
               adsorbate_gas_reference_energy(:H, G) -
               adsorbate_gas_reference_energy(:COOH, G)
    elseif i == 11
        return E[:CO2] + E_OH - E_COOH - E_O +
               adsorbate_gas_reference_energy(:CO2, G) +
               adsorbate_gas_reference_energy(:OH, G) -
               adsorbate_gas_reference_energy(:COOH, G) -
               adsorbate_gas_reference_energy(:O, G)
    elseif i == 12
        return E[:CO2] + E[:H2O] - E_COOH - E_OH +
               adsorbate_gas_reference_energy(:CO2, G) +
               adsorbate_gas_reference_energy(:H2O, G) -
               adsorbate_gas_reference_energy(:COOH, G) -
               adsorbate_gas_reference_energy(:OH, G)
    elseif i == 14
        return E_CO + E_H - E_HCO +
               adsorbate_gas_reference_energy(:CO, G) +
               adsorbate_gas_reference_energy(:H, G) -
               adsorbate_gas_reference_energy(:HCO, G)
    end

    error("No active WGS DFT reaction energy implemented for reaction $i.")
end

function wgs_reaction5_enthalpy(E, consts; use_coverage_dependent_thermo=false)
    H_OH_ads = adsorbate_enthalpy_ev(:OH, E, consts, use_coverage_dependent_thermo)
    H_H_ads = adsorbate_enthalpy_ev(:H, E, consts, use_coverage_dependent_thermo)
    H_H2O_ads = adsorbate_enthalpy_ev(:H2O, E, consts, use_coverage_dependent_thermo)
    return H_OH_ads + H_H_ads - H_H2O_ads
end

function wgs_reaction_enthalpy(i, E, consts; use_coverage_dependent_thermo=false)
    if i == 5
        return adsorbate_enthalpy_ev(:OH, E, consts, use_coverage_dependent_thermo) +
               adsorbate_enthalpy_ev(:H, E, consts, use_coverage_dependent_thermo) -
               adsorbate_enthalpy_ev(:H2O, E, consts, use_coverage_dependent_thermo)
    elseif i == 6
        return adsorbate_enthalpy_ev(:O, E, consts, use_coverage_dependent_thermo) +
               adsorbate_enthalpy_ev(:H, E, consts, use_coverage_dependent_thermo) -
               adsorbate_enthalpy_ev(:OH, E, consts, use_coverage_dependent_thermo)
    elseif i == 7
        return adsorbate_enthalpy_ev(:H2O, E, consts, use_coverage_dependent_thermo) +
               adsorbate_enthalpy_ev(:O, E, consts, use_coverage_dependent_thermo) -
               2.0 * adsorbate_enthalpy_ev(:OH, E, consts, use_coverage_dependent_thermo)
    elseif i == 8
        return adsorbate_enthalpy_ev(:CO2, E, consts, use_coverage_dependent_thermo) -
               adsorbate_enthalpy_ev(:CO, E, consts, use_coverage_dependent_thermo) -
               adsorbate_enthalpy_ev(:O, E, consts, use_coverage_dependent_thermo)
    elseif i == 9
        return adsorbate_enthalpy_ev(:COOH, E, consts, use_coverage_dependent_thermo) -
               adsorbate_enthalpy_ev(:CO, E, consts, use_coverage_dependent_thermo) -
               adsorbate_enthalpy_ev(:OH, E, consts, use_coverage_dependent_thermo)
    elseif i == 10
        return adsorbate_enthalpy_ev(:CO2, E, consts, use_coverage_dependent_thermo) +
               adsorbate_enthalpy_ev(:H, E, consts, use_coverage_dependent_thermo) -
               adsorbate_enthalpy_ev(:COOH, E, consts, use_coverage_dependent_thermo)
    elseif i == 11
        return adsorbate_enthalpy_ev(:CO2, E, consts, use_coverage_dependent_thermo) +
               adsorbate_enthalpy_ev(:OH, E, consts, use_coverage_dependent_thermo) -
               adsorbate_enthalpy_ev(:COOH, E, consts, use_coverage_dependent_thermo) -
               adsorbate_enthalpy_ev(:O, E, consts, use_coverage_dependent_thermo)
    elseif i == 12
        return adsorbate_enthalpy_ev(:CO2, E, consts, use_coverage_dependent_thermo) +
               adsorbate_enthalpy_ev(:H2O, E, consts, use_coverage_dependent_thermo) -
               adsorbate_enthalpy_ev(:COOH, E, consts, use_coverage_dependent_thermo) -
               adsorbate_enthalpy_ev(:OH, E, consts, use_coverage_dependent_thermo)
    elseif i == 14
        return adsorbate_enthalpy_ev(:CO, E, consts, use_coverage_dependent_thermo) +
               adsorbate_enthalpy_ev(:H, E, consts, use_coverage_dependent_thermo) -
               adsorbate_enthalpy_ev(:HCO, E, consts, use_coverage_dependent_thermo)
    end

    error("No active WGS reaction enthalpy implemented for reaction $i.")
end

function wgs_reaction_omega(i)
    if i == 5
        return WGS_OMEGA_RXN5
    elseif i == 6
        return WGS_OMEGA_RXN6
    elseif i == 7
        return WGS_OMEGA_RXN7
    elseif i == 8
        return WGS_OMEGA_RXN8
    elseif i == 9
        return WGS_OMEGA_RXN9
    elseif i == 10
        return WGS_OMEGA_RXN10
    elseif i == 11
        return WGS_OMEGA_RXN11
    elseif i == 12
        return WGS_OMEGA_RXN12
    elseif i == 14
        return WGS_OMEGA_RXN14
    end

    error("No omega implemented for reaction $i.")
end

function wgs_reaction_neutral_forward_barrier(i)
    if i == 6
        return WGS_EF0_RXN6
    elseif i == 7
        return WGS_EF0_RXN7
    elseif i == 11
        return WGS_EF0_RXN11
    elseif i == 12
        return WGS_EF0_RXN12
    elseif i == 14
        return WGS_EF0_RXN14
    end

    error("No neutral forward barrier implemented for reaction $i.")
end

function wgs_gas_phase_energies(descriptors)
    G_H2O = descriptor_value(descriptors, :G_H2O, 0.0)
    G_CO  = descriptor_value(descriptors, :G_CO,  0.0)
    G_H2  = descriptor_value(descriptors, :G_H2,  NaN)
    G_O2  = descriptor_value(descriptors, :G_O2,  NaN)
    G_CO2 = descriptor_value(descriptors, :G_CO2, 0.0)
    G_CH   = descriptor_value(descriptors, :G_CH,   0.0)
    G_H    = descriptor_value(descriptors, :G_H,    0.0)
    G_OH   = descriptor_value(descriptors, :G_OH,   0.0)
    G_O    = descriptor_value(descriptors, :G_O,    0.0)
    G_COOH = descriptor_value(descriptors, :G_COOH, 0.0)
    G_HCO  = descriptor_value(descriptors, :G_HCO,  0.0)

    return Dict{Symbol, Float64}(
        :H2O  => G_H2O,
        :CO   => G_CO,
        :H2   => G_H2,
        :O2   => G_O2,
        :CO2  => G_CO2,
        :CH   => G_CH,
        :H    => G_H,
        :OH   => G_OH,
        :O    => G_O,
        :COOH => G_COOH,
        :HCO  => G_HCO,
    )
end

function surface_energy(E, species::Symbol, use_coverage_dependent_thermo)
    if use_coverage_dependent_thermo
        cov_species = Symbol(string(species), "_cov")
        if haskey(E, cov_species)
            return E[cov_species]
        end
    end

    return E[species]
end

function is_neutral_base_descriptor_set(E; atol=1.0e-8)
    return abs(E[:CH] - (-6.9168)) <= atol &&
           abs(E[:O] - (-1.2550)) <= atol &&
           abs(E[:H2O] - (-0.24235)) <= atol &&
           abs(E[:HCO] - (-2.36384)) <= atol
end

function neutral_base_descriptors(descriptors)
    return (
        CH = -6.9168,
        O = -1.2550,
        H2O = -0.24235,
        HCO = -2.36384,
        G_H2O = descriptor_value(descriptors, :G_H2O, 0.0),
        G_CO = descriptor_value(descriptors, :G_CO, 0.0),
        G_H2 = descriptor_value(descriptors, :G_H2, NaN),
        G_O2 = descriptor_value(descriptors, :G_O2, NaN),
        G_CO2 = descriptor_value(descriptors, :G_CO2, 0.0),
        G_CH = descriptor_value(descriptors, :G_CH, 0.0),
        G_H = descriptor_value(descriptors, :G_H, 0.0),
        G_OH = descriptor_value(descriptors, :G_OH, 0.0),
        G_O = descriptor_value(descriptors, :G_O, 0.0),
        G_COOH = descriptor_value(descriptors, :G_COOH, 0.0),
        G_HCO = descriptor_value(descriptors, :G_HCO, 0.0),
        coverage_dependent_thermo = descriptor_value(descriptors, :coverage_dependent_thermo, false),
    )
end

function wgs_reaction_energies(E, G, consts; use_coverage_dependent_thermo=false)
    E_H = surface_energy(E, :H, use_coverage_dependent_thermo)
    E_O = surface_energy(E, :O, use_coverage_dependent_thermo)
    E_OH = surface_energy(E, :OH, use_coverage_dependent_thermo)
    E_CO = surface_energy(E, :CO, use_coverage_dependent_thermo)
    E_COOH = surface_energy(E, :COOH, use_coverage_dependent_thermo)
    T = promote_type(
        typeof(E_H),
        typeof(E_O),
        typeof(E_OH),
        typeof(E_CO),
        typeof(E_COOH),
        typeof(E[:H2O]),
        typeof(E[:CO2]),
        typeof(E[:HCO]),
    )
    Erxn = zeros(T, 14)

    Erxn[1]  = E[:H2O]
    Erxn[2]  = E[:CO]
    Erxn[3]  = 2.0 * E[:H]
    Erxn[4]  = E[:CO2]
    Erxn[5]  = wgs_reaction_enthalpy(5, E, consts; use_coverage_dependent_thermo=use_coverage_dependent_thermo)
    Erxn[6]  = wgs_reaction_enthalpy(6, E, consts; use_coverage_dependent_thermo=use_coverage_dependent_thermo)
    Erxn[7]  = wgs_reaction_enthalpy(7, E, consts; use_coverage_dependent_thermo=use_coverage_dependent_thermo)
    Erxn[8]  = wgs_reaction_enthalpy(8, E, consts; use_coverage_dependent_thermo=use_coverage_dependent_thermo)
    Erxn[9]  = wgs_reaction_enthalpy(9, E, consts; use_coverage_dependent_thermo=use_coverage_dependent_thermo)
    Erxn[10] = wgs_reaction_enthalpy(10, E, consts; use_coverage_dependent_thermo=use_coverage_dependent_thermo)
    Erxn[11] = wgs_reaction_enthalpy(11, E, consts; use_coverage_dependent_thermo=use_coverage_dependent_thermo)
    Erxn[12] = wgs_reaction_enthalpy(12, E, consts; use_coverage_dependent_thermo=use_coverage_dependent_thermo)
    Erxn[14] = wgs_reaction_enthalpy(14, E, consts; use_coverage_dependent_thermo=use_coverage_dependent_thermo)

    return Erxn
end

function wgs_forward_barriers(E, G, Erxn, consts, Erxn_neutral; use_coverage_dependent_thermo=false)
    T = promote_type(eltype(Erxn), typeof(E[:TS_HO_H]), typeof(E[:TS_O_CO]), typeof(E[:TS_OC_OH]), typeof(E[:TS_H_CO2]))
    Eact = zeros(T, 14)

    # Adsorption is treated as barrierless in the forward direction.
    Eact[1:4] .= 0.0

    # These scaling relations return DFT forward activation energies directly.
    Eact_DFT = Dict(
        5 => E[:TS_HO_H],
        8 => E[:TS_O_CO],
        9 => E[:TS_OC_OH],
        10 => E[:TS_H_CO2],
    )

    for i in (5, 8, 9, 10)
        delta_E_DFT = wgs_reaction_dft_energy(i, E, G; use_coverage_dependent_thermo=use_coverage_dependent_thermo)
        Eact[i] = Eact_DFT[i] + wgs_reaction_omega(i) * (Erxn[i] - delta_E_DFT)
    end

    for i in (6, 7, 11, 12, 14)
        omega = wgs_reaction_omega(i)
        Eact_DFT_q = wgs_reaction_neutral_forward_barrier(i) + omega * (Erxn[i] - Erxn_neutral[i])
        delta_E_DFT = wgs_reaction_dft_energy(i, E, G; use_coverage_dependent_thermo=use_coverage_dependent_thermo)
        Eact[i] = Eact_DFT_q + omega * (Erxn[i] - delta_E_DFT)
    end

    return Eact
end

function default_wgs_consts_for_energetics()
    return (
        kb   = 1.380649e-23,
        h    = 6.62607015e-34,
        Na   = 6.02214076e23,
        kbeV = 8.617333262e-5,
        Rgas = 8.314462618,
        Te   = 523.0,
        No   = 1.1e19,
    )
end

function wgs_energetics(descriptors, consts=default_wgs_consts_for_energetics(); theta_CO=0.0)
    E = wgs_scaled_energies(descriptors; theta_CO=theta_CO)
    G = wgs_gas_phase_energies(descriptors)
    use_coverage_dependent_thermo = descriptor_value(descriptors, :coverage_dependent_thermo, false)
    Erxn = wgs_reaction_energies(E, G, consts; use_coverage_dependent_thermo=use_coverage_dependent_thermo)

    neutral_descriptors = neutral_base_descriptors(descriptors)
    E_neutral = wgs_scaled_energies(neutral_descriptors; theta_CO=theta_CO)
    G_neutral = wgs_gas_phase_energies(neutral_descriptors)
    neutral_coverage_dependent_thermo = descriptor_value(neutral_descriptors, :coverage_dependent_thermo, false)
    Erxn_neutral = wgs_reaction_energies(
        E_neutral,
        G_neutral,
        consts;
        use_coverage_dependent_thermo=neutral_coverage_dependent_thermo,
    )

    Eact = wgs_forward_barriers(
        E,
        G,
        Erxn,
        consts,
        Erxn_neutral;
        use_coverage_dependent_thermo=use_coverage_dependent_thermo,
    )

    return E, Erxn, Eact
end

function make_wgs_rate_constants(descriptors, consts; theta_CO=0.0)
    E, Erxn, Eact = wgs_energetics(descriptors, consts; theta_CO=theta_CO)

    MW_H2O = 18.01528
    MW_CO  = 28.0101
    MW_H2  = 2.01588
    MW_CO2 = 44.0095

    T = promote_type(eltype(Erxn), eltype(Eact), typeof(E[:CO_ads]), typeof(E[:H_ads]))
    kf = zeros(T, 14)
    kr = zeros(T, 14)
    keq = zeros(T, 14)

    kf[1],  kr[1],  keq[1]  = make_adsorption_sticking_constants(MW_H2O, Erxn[1], Eact[1], WGS_STICKING_H2O, consts)
    kf[2],  kr[2],  keq[2]  = make_adsorption_sticking_constants(MW_CO, E[:CO_ads], Eact[2], WGS_STICKING_CO, consts)
    kf[3],  kr[3],  keq[3]  = make_adsorption_sticking_constants(MW_H2, 2.0 * E[:H_ads], Eact[3], WGS_STICKING_H2, consts)
    kf[4],  kr[4],  keq[4]  = make_adsorption_sticking_constants(MW_CO2, Erxn[4], Eact[4], WGS_STICKING_CO2, consts)

    for i in (5, 6, 7, 8, 9, 10, 11, 12, 14)
        kf[i], kr[i], keq[i] = makeks(1.0, 1.0, 1.0, Erxn[i], Eact[i], consts)
    end

    return kf, kr, keq, E, Erxn, Eact
end
