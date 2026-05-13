module SubModuleLinearSolveDense

using Test
using F1Method
using LinearAlgebra
using SciMLBase
using ForwardDiff
using LinearSolve

# Reuse the quasi-Rosenbrock setup from test/rosenbrock.jl: same statefun,
# same ODEFunction, same MyAlg Newton solver, same mismatch f / ∇ₓf. The
# only thing we vary is the `linsolve` kwarg on F1Cache.

function newton_solve(F, ∇ₓF, x; Ftol = 1e-10)
    while norm(F(x)) ≥ Ftol
        x .-= ∇ₓF(x) \ F(x)
    end
    return x
end

struct MyAlg <: SciMLBase.AbstractSteadyStateAlgorithm end

function SciMLBase.solve(prob::SciMLBase.AbstractSteadyStateProblem,
                          alg::MyAlg;
                          Ftol = 1e-10)
    p = prob.p
    x0 = copy(prob.u0)
    F(x) = prob.f(x, p)
    ∇ₓF(x) = prob.f.jac(x, p)
    x_steady = newton_solve(F, ∇ₓF, x0; Ftol = Ftol)
    resid = F(x_steady)
    return SciMLBase.build_solution(prob, alg, x_steady, resid; retcode = ReturnCode.Success)
end

statefun(x, p, t = 0) = [
    -2 * (p[1] - x[1]) - 4 * p[2] * (x[2] - x[1]^2) * x[1]
    p[2] * (x[2] - x[1]^2)
]
F = ODEFunction{false}(statefun; jac = (x, p, t = 0) -> ForwardDiff.jacobian(x -> statefun(x, p), x))

function state_mismatch(x)
    δ(x) = x .- 1
    return 0.5δ(x)'δ(x)
end
function parameter_mismatch(p)
    δ(p) = log.(p)
    return 0.5δ(p)'δ(p)
end
f(x, p) = state_mismatch(x) + parameter_mismatch(p)
∇ₓf(x, p) = ForwardDiff.jacobian(x -> [f(x, p)], x)

x₀ = rand(2)
p₀ = rand(2)

exact_solution(p)  = [p[1], p[1]^2]
exact_objective(p) = f(exact_solution(p), p)

# `linsolve = nothing` skips LinearSolve entirely and uses Julia's
# `factorize`, which on a dense `Matrix{Float64}` returns a pivoted
# LAPACK `LU` (same family as LUFactorization, but invoked directly
# through LinearAlgebra rather than via a LinearSolve.LinearCache).
for linsolve in (nothing, LUFactorization(), GenericLUFactorization())
    @testset "linsolve = $(linsolve)" begin
        cache = F1Method.F1Cache(F, ∇ₓf, x₀, p₀, MyAlg(); linsolve = linsolve)
        F1_objective(p) = F1Method.objective(f, F, cache, p, MyAlg())
        F1_gradient(p)  = F1Method.gradient(f, F, ∇ₓf, cache, p, MyAlg())
        F1_hessian(p)   = F1Method.hessian(f, F, ∇ₓf, cache, p, MyAlg())

        @test exact_objective(p₀)                       ≈ F1_objective(p₀)
        @test exact_objective(2p₀)                      ≈ F1_objective(2p₀)
        @test ForwardDiff.gradient(exact_objective, p₀)  ≈ F1_gradient(p₀)
        @test ForwardDiff.gradient(exact_objective, 2p₀) ≈ F1_gradient(2p₀)
        @test ForwardDiff.hessian(exact_objective, p₀)   ≈ F1_hessian(p₀)
        @test ForwardDiff.hessian(exact_objective, 2p₀)  ≈ F1_hessian(2p₀)
    end
end

end # submodule
