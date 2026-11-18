function wgs_net_rates(x, p)
    consts, Rg, V, q, Nsites, Cf_H2O, Cf_CO, Cf_H2, Cf_CO2, descriptors = p

    C_H2O = x[1]
    C_CO  = x[2]
    C_H2  = x[3]
    C_CO2 = x[4]

    N_H2O  = x[5]
    N_CO   = x[6]
    N_H    = x[7]
    N_CO2  = x[8]
    N_OH   = x[9]
    N_O    = x[10]
    N_COOH = x[11]
    N_HCO  = x[12]

    N_star = Nsites - N_H2O - N_CO - N_H - N_CO2 - N_OH - N_O - N_COOH - N_HCO

    th_H2O  = N_H2O  / Nsites
    th_CO   = N_CO   / Nsites
    th_H    = N_H    / Nsites
    th_CO2  = N_CO2  / Nsites
    th_OH   = N_OH   / Nsites
    th_O    = N_O    / Nsites
    th_COOH = N_COOH / Nsites
    th_HCO  = N_HCO  / Nsites
    th_star = N_star / Nsites

    pH2O = C_H2O * Rg * consts.Te
    pCO  = C_CO  * Rg * consts.Te
    pH2  = C_H2  * Rg * consts.Te
    pCO2 = C_CO2 * Rg * consts.Te

    kf, kr, _, _, _, _ = make_wgs_rate_constants(descriptors, consts; theta_CO=th_CO)
    # Legacy array slot 13 is disabled. Julia slot 14 is manuscript R13,
    # HCO* + * <=> CO* + H*.
    inactive_rate = zero(kf[1] * pH2O * th_star)

    rf = [
        kf[1]  * pH2O * th_star,
        kf[2]  * pCO  * th_star,
        kf[3]  * pH2  * th_star^2,
        kf[4]  * pCO2 * th_star,
        kf[5]  * th_H2O  * th_star,
        kf[6]  * th_OH   * th_star,
        kf[7]  * th_OH^2,
        kf[8]  * th_CO   * th_O,
        kf[9]  * th_CO   * th_OH,
        kf[10] * th_COOH * th_star,
        kf[11] * th_COOH * th_O,
        kf[12] * th_COOH * th_OH,
        inactive_rate,
        kf[14] * th_HCO  * th_star,
    ]

    rr = [
        kr[1]  * th_H2O,
        kr[2]  * th_CO,
        kr[3]  * th_H^2,
        kr[4]  * th_CO2,
        kr[5]  * th_OH   * th_H,
        kr[6]  * th_O    * th_H,
        kr[7]  * th_H2O  * th_O,
        kr[8]  * th_CO2  * th_star,
        kr[9]  * th_COOH * th_star,
        kr[10] * th_CO2  * th_H,
        kr[11] * th_CO2  * th_OH,
        kr[12] * th_CO2  * th_H2O,
        inactive_rate,
        kr[14] * th_CO   * th_H,
    ]

    rnet = rf .- rr

    return rnet
end

function wgs_surface_tof(x, p)
    rnet = wgs_net_rates(x, p)
    return rnet[8] + rnet[10] + rnet[11] + rnet[12]
end

function WGS_cstr!(xdot, x, p, t)
    consts, Rg, V, q, Nsites, Cf_H2O, Cf_CO, Cf_H2, Cf_CO2, descriptors = p

    C_H2O = x[1]
    C_CO  = x[2]
    C_H2  = x[3]
    C_CO2 = x[4]

    N_H2O  = x[5]
    N_CO   = x[6]
    N_H    = x[7]
    N_CO2  = x[8]
    N_OH   = x[9]
    N_O    = x[10]
    N_COOH = x[11]
    N_HCO  = x[12]

    rnet = wgs_net_rates(x, p)

    dN_H2O  = Nsites * ( rnet[1] - rnet[5] + rnet[7] + rnet[12] )
    dN_CO   = Nsites * ( rnet[2] - rnet[8] - rnet[9] + rnet[14] )
    dN_H    = Nsites * ( 2.0*rnet[3] + rnet[5] + rnet[6] + rnet[10] + rnet[14] )
    dN_CO2  = Nsites * ( rnet[4] + rnet[8] + rnet[10] + rnet[11] + rnet[12] )
    dN_OH   = Nsites * ( rnet[5] - rnet[6] - 2.0*rnet[7] - rnet[9] + rnet[11] - rnet[12] )
    dN_O    = Nsites * ( rnet[6] + rnet[7] - rnet[8] - rnet[11] )
    dN_COOH = Nsites * ( rnet[9] - rnet[10] - rnet[11] - rnet[12] )
    dN_HCO  = Nsites * ( -rnet[14] )

    dC_H2O = (q / V) * (Cf_H2O - C_H2O) + (Nsites / V) * (-rnet[1])
    dC_CO  = (q / V) * (Cf_CO  - C_CO)  + (Nsites / V) * (-rnet[2])
    dC_H2  = (q / V) * (Cf_H2  - C_H2)  + (Nsites / V) * (-rnet[3])
    dC_CO2 = (q / V) * (Cf_CO2 - C_CO2) + (Nsites / V) * (-rnet[4])

    xdot[1]  = dC_H2O
    xdot[2]  = dC_CO
    xdot[3]  = dC_H2
    xdot[4]  = dC_CO2
    xdot[5]  = dN_H2O
    xdot[6]  = dN_CO
    xdot[7]  = dN_H
    xdot[8]  = dN_CO2
    xdot[9]  = dN_OH
    xdot[10] = dN_O
    xdot[11] = dN_COOH
    xdot[12] = dN_HCO

    return nothing
end
