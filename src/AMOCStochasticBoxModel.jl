"""
Atlantic meridional overturning circulation (AMOC) stochastic box model.

## Details

AMOC stochastic box model described in Soons, Grafke & Dijkstra (2024). This is a
stochastic extension of the five compartment AMOC box model described in Wood et al.
(2019).

$(EXPORTS)
    
## References:

1. Soons, J., Grafke, T., & Dijkstra, H. A. (2024). Optimal transition paths for AMOC
   collapse and recovery in a stochastic box model. Journal of Physical Oceanography,
   54(12), 2537-2552.
2. Wood, R. A., Rodríguez, J. M., Smith, R. S., Jackson, L. C., & Hawkins, E. (2019). 
   Observable, low-order dynamical controls on thresholds of the Atlantic meridional
   overturning circulation. Climate Dynamics, 53, 6A815-6834.
"""
module AMOCStochasticBoxModel

using StochasticDiffEq
using DiffEqNoiseProcess
using DocStringExtensions
using Plots
using Random
using SpecialFunctions

export AMOCStochasticBoxModelParameters, simulate_model, plot_solution

"""
$(TYPEDEF)

Parameters for AMOC stochastic box model.
    
$(TYPEDSIGNATURES)
    
## Details
    
Default values for parameters are taken from Soons, Grafke & Dijkstra (2024) with
exception of smoothed step function scale parameter `ϵ_θ` which defaults to 1e-2
(compared to 1e-10 in paper) as this was found to have negligible effect on simulated
paths while increasing smoothness of map from noise process to simulated paths, and
noise process variance parameter `ν` which is set to `0.2^2` (various different values
are used in paper) which was found to give sufficient variation in stochastic forcing to
cause high probability of switch between bistable `φ_ON` and `φ_OFF` states over
simulated time spans of order 500 years.

$(TYPEDFIELDS)
"""
@kwdef struct AMOCStochasticBoxModelParameters{T<:Real}
    "Volume scaling parameter / m³"
    V_0::T = 1e10
    "Non-dimensional volume of box 1 (northern, N)"
    V_1::T = 3.261
    "Non-dimensional volume of box 2 (Atlantic thermocline, T)"
    V_2::T = 7.777
    "Non-dimensional volume of box 3 (southern, S)"
    V_3::T = 8.897
    "Non-dimensional volume of box 4 (Indo-Pacific, IP)"
    V_4::T = 22.02
    "Non-dimensional volume of box 5 (bottom, B)"
    V_5::T = 86.49
    "Time scaling parameter / s"
    t_d::T = 3.1536e9
    "Salinity scaling parameter / psu"
    S_0::T = 35.0
    "Total salinity in basin / 10⁶⋅m³⋅psu"
    C::T = 4.446304026e13
    "Thermal coefficient / kg⋅m⁻³⋅K"
    α::T = 0.12
    "Saline coefficient / kg⋅m⁻³"
    β::T = 0.79
    "Proportion cold water path"
    γ::T = 0.39
    "S-B box mixing parameter / Sv"
    η::T = 74.492
    "Northern subtropical gyre coefficient / Sv"
    K_N::T = 5.456
    "Southern subtropical gyre coefficient / Sv"
    K_S::T = 5.447
    "Indo-Pacific gyre coefficient / Sv"
    K_IP::T = 96.817
    "Northern freshwater flux / Sv"
    F_N::T = 0.384
    "Atlantic thermocline freshwater flux / Sv"
    F_T::T = -0.723
    "Southern freshwater flux / Sv"
    F_S::T = 1.078
    "Indo-Pacific freshwater flux / Sv"
    F_IP::T = -0.739
    "Southern box temperature / K"
    T_S::T = 4.773
    "Base temperature / K"
    T_0::T = 2.65
    "Heat transport coefficient / 10⁻⁶⋅m⁻³sK"
    μ::T = 0.055
    "MOC-density difference coefficient / 10⁶⋅m⁶⋅(kg)⁻¹⋅s"
    λ::T = 27.9
    "Fraction of stochastic fresh water flux entering northern box"
    A_N::T = 0.070
    "Fraction of stochastic fresh water flux entering Atlantic thermocline box"
    A_T::T = 0.752
    "Fraction of stochastic fresh water flux entering southern box"
    A_S::T = -0.257
    "Fraction of stochastic fresh water flux entering Indo-Pacific box"
    A_IP::T = -0.565
    "Stochastic freshwater forcing non-dimensional process variance"
    ν::T = 0.2^2
    "Smoothed step function scale parameter"
    ϵ_θ::T = 1e-2
end

"""
$(SIGNATURES)

Smoothed step function in variable `q` with scale parameter `ϵ_θ`.
    
Converges to Heaviside step function in `q` in limit `ϵ_θ → 0`.
"""
θ(q, ϵ_θ) = (erf(q / ϵ_θ) + 1) / 2

"""
$(SIGNATURES)
    
Compute `κ` coefficient in definition of AMOC strength given model parameters `p`.
"""
κ(p) = p.λ / (1 + p.α * p.λ * p.μ)

function ϕ₅(φ_1, φ_2, φ_3, φ_4, p)
    (
        p.C / (p.V_0 * p.S_0) - p.V_1 * φ_1 - p.V_2 * φ_2 - p.V_3 * φ_3 - p.V_4 * φ_4
    ) / p.V_5
end

"""
$(SIGNATURES)

Compute non-dimensional state component φ₅ corresponding to salinity of bottom (B) box
from other state components `φ` and parameters `p`.
"""
ϕ₅(φ, p) = ϕ₅(φ..., p)

function amoc_strength(φ_1, φ_2, φ_3, φ_4, p)
    κ(p) * (p.α * (p.T_S - p.T_0) + p.β * p.S_0 * (φ_1 - φ_3))
end

"""
$(SIGNATURES)

AMOC strength (in Sv) as a function of non-dimensional state `φ` and parameters `p`.

Corresponds to downwelling in North Atlantic that transports salt between Atlantic
thermocline (T), northern (N) and bottom (B) boxes in model. 
"""
amoc_strength(φ, p) = amoc_strength(φ..., p)

"""
$(SIGNATURES)
    
Solve for a non-dimensional state `φ` which gives a target AMOC strength `a` given
parameters `p` and values for a subset of the state components.
"""
function solve_for_state_with_target_amoc_strength(a, p; φ_1=1, φ_2=1, φ_4=1)
    φ_3 = p.α * (p.T_S - p.T_0) / (p.β * p.S_0) + φ_1 - a / (κ(p) * p.β * p.S_0)
    [φ_1, φ_2, φ_3, φ_4]
end

"""
$(SIGNATURES)

Compute drift term in AMOC stochastic box model stochastic differential equation and
write in place to vector `dφ` given current non-dimensional state `φ`, parameters `p`
and time `t`.
"""
function drift!(dφ, φ, p, t)
    q = amoc_strength(φ, p)
    φ_5 = ϕ₅(φ, p)
    dφ[1] = (p.t_d / (p.V_0 * p.V_1)) * (
        q * (θ(q, p.ϵ_θ) * (φ[2] - φ[1]) - θ(-q, p.ϵ_θ) * (φ_5 - φ[1]))
        +
        p.K_N * (φ[2] - φ[1])
        -
        p.F_N
    )
    dφ[2] = (p.t_d / (p.V_0 * p.V_2)) * (
        q * (
            θ(q, p.ϵ_θ) * (p.γ * φ[3] + (1 - p.γ) * φ[4] - φ[2])
            -
            θ(-q, p.ϵ_θ) * (φ[1] - φ[2])
        )
        + p.K_S * (φ[3] - φ[2])
        + p.K_N * (φ[1] - φ[2])
        -
        p.F_T
    )
    dφ[3] = (p.t_d / (p.V_0 * p.V_3)) * (
        p.γ * q * (θ(q, p.ϵ_θ) * (φ_5 - φ[3]) - θ(-q, p.ϵ_θ) * (φ[2] - φ[3]))
        + p.K_IP * (φ[4] - φ[3])
        + p.K_S * (φ[2] - φ[3])
        + p.η * (φ_5 - φ[3])
        -
        p.F_S
    )
    dφ[4] = (p.t_d / (p.V_0 * p.V_4)) * (
        (1 - p.γ) * q * (θ(q, p.ϵ_θ) * (φ_5 - φ[4]) - θ(-q, p.ϵ_θ) * (φ[2] - φ[4]))
        +
        p.K_IP * (φ[3] - φ[4])
        -
        p.F_IP
    )
end

"""
$(SIGNATURES)

Compute diffusion-coefficient term in AMOC stochastic box model stochastic differential
equation and write in place to vector `dw` given current non-dimensional state `φ`,
parameters `p` and time `t`.
"""
function diffusion_coefficient!(dw, φ, p, t)
    dw[1] = p.A_N / p.V_1
    dw[2] = p.A_T / p.V_2
    dw[3] = p.A_S / p.V_3
    dw[4] = p.A_IP / p.V_4
    dw .*= (p.t_d / p.V_0) * √p.ν
end

"""
$(SIGNATURES)

Simulate AMOC stochastic box model with parameters `parameters` over time interval
`time_interval` with state initialised to achieve a target AMOC strength
`target_initial_amoc_strength` using stochastic differential equation solver algorithm
`solver_algorithm` and random seed `seed`.
"""
function simulate_model(;
    parameters=AMOCStochasticBoxModelParameters(),
    time_interval=(0., 500.),
    target_initial_amoc_strength=0.,
    solver_algorithm=SOSRI(),
    seed=1234,
)
    φ_0 = solve_for_state_with_target_amoc_strength(
        target_initial_amoc_strength, parameters
    )
    noise = WienerProcess(time_interval[1], 0., 0.)
    problem = SDEProblem(
        drift!, diffusion_coefficient!, φ_0, time_interval, parameters; noise, seed
    )
    solve(problem, solver_algorithm)
end

"""
$(SIGNATURES)

Plot solution `solution` of AMOC stochastic box model with figure size `fig_size`.

Plots three-dimensional projections of non-dimensional state trajectories in first
panel, simulated box salinities over time in second panel, simulated AMOC strength over
time in third panel and simulated stochastic freshwater forcing over time in last panel.
"""
function plot_solution(solution; fig_size=(800, 900))
    parameters = solution.prob.p
    projections = [
        plot(
            solution,
            xlabel="\$\\phi_$(idxs[1])\$",
            ylabel="\$\\phi_$(idxs[2])\$",
            zlabel="\$\\phi_$(idxs[3])\$",
            idxs=idxs,
            legend=false,
            ticks=nothing,
        )
        for idxs in [(1, 2, 3), (1, 2, 4), (1, 3, 4), (2, 3, 4)]
    ]
    salinity_vs_time = plot(
        solution,
        idxs=[((t, φ) -> parameters.S_0 * φ, 0, i) for i in 1:4],
        labels=["\$S_N\$" "\$S_T\$" "\$S_S\$" "\$S_{IP}\$"]
    )
    plot!(
        salinity_vs_time,
        solution,
        idxs=[((t, φ...) -> parameters.S_0 * ϕ₅(φ..., parameters), 0, 1, 2, 3, 4)],
        label="\$S_B\$",
        xlabel="\$t\$ / years",
        ylabel="Salinity / psu"
    )
    amoc_vs_time = plot(
        solution,
        idxs=[((t, φ...) -> amoc_strength(φ..., parameters), 0, 1, 2, 3, 4)],
        xlabel="\$t\$ / years",
        ylabel="AMOC strength / Sv",
        legend=false
    )
    dw = zeros(4)
    diffusion_coefficient!(dw, nothing, parameters, nothing)
    scale_factor = parameters.S_0 * sum(dw .^ 2) .^ 0.5
    forcing_vs_time = plot(
        solution.t,
        solution.W.u .* scale_factor,
        xlabel="\$t\$ / years",
        ylabel="Forcing / Sv",
        label=nothing
    )
    plot(
        plot(projections..., layout=(1, 4)),
        salinity_vs_time,
        amoc_vs_time,
        forcing_vs_time,
        layout=grid(4, 1, heights=[0.3, 0.3, 0.2, 0.2]),
        size=fig_size
    )
end

end
