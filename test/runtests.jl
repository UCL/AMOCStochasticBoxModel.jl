using AMOCStochasticBoxModel
using Test

@testset "AMOCStochasticBoxModel.jl" begin
    parameters = AMOCStochasticBoxModelParameters()
    solution = simulate_model(; parameters)
    @test !any(u -> any(isnan.(u)), solution.u)
    plot_output = plot_solution(solution)
    @test !isnothing(plot_output)
    for a in [-0.5, 0., 1.]
        ϕ_0 = AMOCStochasticBoxModel.solve_for_state_with_target_amoc_strength(
            a, parameters
        )
        @test AMOCStochasticBoxModel.amoc_strength(ϕ_0, parameters) ≈ a atol=1e-12
    end
    @test AMOCStochasticBoxModel.θ(-1., 1e-10) ≈ 0
    @test AMOCStochasticBoxModel.θ(+1., 1e-10) ≈ 1
end
