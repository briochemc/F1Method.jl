
using Test, F1Method

using LinearAlgebra
using DiffEqBase
using ForwardDiff

using AIBECS
using Unitful: m, d, s, yr, Myr, mol, mmol, μmol, μM
using Distributions
using WorldOceanAtlasTools
using FiniteDiff

# Set up:
# - overload `SteadyStateProblem` constructor
# - overload `solve` function
# - define solver algorithm (basic Newton here)
# - define type for that algorithm (here `MyAlg`)
include("simple_setup.jl")

@testset "quasi-Rosenbrock derivative" begin
    include("rosenbrock.jl")
end

@testset "AIBECS test" begin
    include("AIBECS_test.jl")
end

