using AMOCStochasticBoxModel
using Test

@testset "AMOCStochasticBoxModel.jl" begin
    parameters = AMOCStochasticBoxModelParameters()
    solution = simulate_model(; parameters)
    @test !any(u -> any(isnan.(u)), solution.u)
end
