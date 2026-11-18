"""Run the Pt WGS microkinetic model at neutral charge.

The four descriptor adsorption energies use the Pt neutral values from the
original kinetic inputs. All other adsorption energies and four DFT-studied
barriers are predicted by the four-descriptor equations in WGS_SR.JL.

Run from the repository root with
    julia --project=MKM MKM/RUN_NEUTRAL.JL
"""

using DifferentialEquations
using OrdinaryDiffEqRosenbrock
using Printf

include(joinpath(@__DIR__, "MAKEKS.JL"))
include(joinpath(@__DIR__, "WGS_ENERGETICS.JL"))
include(joinpath(@__DIR__, "WGS_CSTR.JL"))

consts = (
    kb = 1.380649e-23,
    h = 6.62607015e-34,
    Na = 6.02214076e23,
    kbeV = 8.617333262e-5,
    Rgas = 8.314462618,
    Te = 548.0,
    No = 1.1e19,
)

const CHARGE = 0.0
const DESCRIPTORS = (
    CH = -6.9168,
    HCO = -2.36384,
    H2O = -0.24235,
    O = -1.2550,
    G_H = -3.363164,
    G_O = -4.9267645,
    G_OH = -7.54949,
    G_CH = -6.07495,
    G_CO = -14.7958,
    G_CO2 = -22.985849,
    G_H2O = -14.218838,
    G_HCO = -17.115986,
    G_COOH = -24.464805,
    G_H2 = 2.0 * -3.363164,
    G_O2 = 2.0 * -4.9267645,
    coverage_dependent_thermo = true,
)

Rg = 8.314459848e-2  # L bar mol^-1 K^-1
P_total = 1.01325     # bar
Cf_H2O = 0.208 * P_total / (Rg * consts.Te)
Cf_CO = 0.145 * P_total / (Rg * consts.Te)
Cf_H2 = 0.0
Cf_CO2 = 0.0
q = 101.4 / 1000.0 / 60.0  # L s^-1
V = q * 1.0                # L; one-second residence time
Nsites = consts.No * 1.0 / consts.Na  # mol; 1 m^2 catalyst area

p = (consts, Rg, V, q, Nsites, Cf_H2O, Cf_CO, Cf_H2, Cf_CO2, DESCRIPTORS)
x0 = zeros(12)
x0[1:4] .= (Cf_H2O, Cf_CO, Cf_H2, Cf_CO2)

startup_time = 0.001
startup = solve(
    ODEProblem(WGS_cstr!, x0, (0.0, startup_time), p),
    Rodas5();
    reltol = 1e-8,
    abstol = 1e-12,
    isoutofdomain = (u, p, t) -> any(u .< 0.0),
    maxiters = 10_000_000,
    dtmax = 1e-4,
)
startup.t[end] >= startup_time || error("Constrained startup stopped early: $(startup.retcode)")

x_init = copy(startup.u[end])
X_prev = clamp((Cf_CO - x_init[2]) / Cf_CO, 0.0, 1.0)
chunk_time = 1e3
max_chunks = 1000
ss_tol_theta = 1e-8
X_tol = 1e-10

function save_result(x_end, X_CO, chunks, max_dtheta_dt)
    tof_s_inv = q * (x_end[3] - Cf_H2) / Nsites
    tof_min_inv = 60.0 * tof_s_inv
    tof_reported = 600.0 * tof_s_inv
    output = joinpath(@__DIR__, "NEUTRAL_RESULT.CSV")
    open(output, "w") do io
        println(io, "CHARGE,MODEL,TOF_H2_S_INV,TOF_H2_MIN_INV,TOF_H2_REPORTED,X_CO,CHUNKS,MAX_DTHETA_DT_S_INV")
        @printf(io, "%.1f,ORIGINAL_JULIA_4D,%.12g,%.12g,%.12g,%.12g,%d,%.12g\n",
                CHARGE, tof_s_inv, tof_min_inv, tof_reported, X_CO, chunks, max_dtheta_dt)
    end
    @printf("Charge = %.1f; reported H2 TOF = %.9g; H2 TOF = %.9g s^-1; CO conversion = %.9g\n",
            CHARGE, tof_reported, tof_s_inv, X_CO)
    println("Saved: ", output)
end

for ichunk in 1:max_chunks
    global x_init, X_prev
    sol = solve(
        ODEProblem(WGS_cstr!, x_init, (0.0, chunk_time), p),
        Rodas5();
        reltol = 1e-8,
        abstol = 1e-12,
        isoutofdomain = (u, p, t) -> any(u .< 0.0),
        maxiters = 10_000_000,
    )
    sol.t[end] >= chunk_time || error("Solver stopped early in chunk $ichunk: $(sol.retcode)")
    x_end = copy(sol.u[end])
    dxdt_end = zeros(12)
    WGS_cstr!(dxdt_end, x_end, p, ichunk * chunk_time + startup_time)
    all(isfinite, x_end) && all(isfinite, dxdt_end) || error("Nonfinite state or derivative")

    X_curr = clamp((Cf_CO - x_end[2]) / Cf_CO, 0.0, 1.0)
    max_dtheta_dt = maximum(abs.(dxdt_end[5:12] ./ Nsites))
    if ichunk > 1 && max_dtheta_dt < ss_tol_theta && abs(X_curr - X_prev) < X_tol
        save_result(x_end, X_curr, ichunk, max_dtheta_dt)
        break
    end
    ichunk == max_chunks && error("Steady state not reached within $max_chunks chunks")
    x_init = x_end
    X_prev = X_curr
end
