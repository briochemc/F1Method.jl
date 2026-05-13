module F1Method

#======================================================================
This package implements the F-1 method as described in the work of
Pasquier et al. (2019). The numbers in parentheses refer to the
equation numbers in the manuscript. A BibTeX citation file is
available in the GitHub repository.
======================================================================#

using LinearAlgebra
using SciMLBase
using ADTypes
using DifferentiationInterface
const DI = DifferentiationInterface

"""
    F1Cache

Memory cache holding the steady-state solution and the linear-algebra
artifacts that the F-1 method reuses across `objective`, `gradient`, and
`hessian` calls at the same parameter point.

# Fields
- `linear_cache`: a linear-solve cache holding the factors of
    `𝐀 = ∇ₓ𝑭(𝒔, 𝒑)` — either a `Factorization` (default)
    or a `LinearSolve.LinearCache` when a `linsolve`
    algorithm is configured
- `s`: the steady-state solution, `𝒔(𝒑)`
- `∇s`: the parameter-Jacobian of the steady state, `∇𝒔(𝒑)`
- `∇ₓf`: the state-Jacobian of the objective at the steady state, `∇ₓ𝑓(𝒔, 𝒑)`
- `p`: the parameter vector that matches `linear_cache`, `∇s`, `∇ₓf` (and `s`)
- `psol`: the parameter vector that matches `s` (may be ahead of `p`
           when only the steady-state solve has been refreshed)
- `ad`: the AD backend used for parameter-side derivatives

Construct with `F1Cache`.
"""
mutable struct F1Cache{C, S, JACS, JACF, P, AD <: AbstractADType}
    linear_cache::C
    s::S
    ∇s::JACS
    ∇ₓf::JACF
    p::P
    psol::P
    ad::AD
end


# Solve-operation helpers on the mathematical operator 𝐀 = ∇ₓ𝑭. The
# `Factorization` branch lives here; the `LinearCache` branch lives in
# F1MethodLinearSolveExt.
set_A!(cache::F1Cache, J)   = cache.linear_cache = factorize(J)
solveA!(cache::F1Cache, b)  = cache.linear_cache \ b
solveAᵀ!(cache::F1Cache, B) = cache.linear_cache' \ B


function update_cache!(F, ∇ₓf, cache, p, alg; options...)
    if p ≠ cache.p                            # only update if 𝒑 has changed
        update_solution!(F, cache, p, alg; options...)
        ∇ₚF = DI.jacobian(p -> F(cache.s, p), cache.ad, p)
        set_A!(cache, F.jac(cache.s, p))      # update factors of ∇ₓ𝑭(𝒔, 𝒑)
        cache.∇s .= solveA!(cache, -∇ₚF)      # update ∇𝒔
        cache.∇ₓf .= ∇ₓf(cache.s, p)          # update ∇ₓ𝑓(𝒔, 𝒑)
        cache.p .= p                          # update 𝒑 for the cached values
    end
    return cache
end

function update_solution!(F, cache, p, alg; options...)
    if p ≠ cache.psol
        prob = SteadyStateProblem(F, cache.s, p)
        cache.s .= solve(prob, alg; options...).u
        cache.psol .= p
    end
    return cache
end

"""
    objective(f, F, cache, p, alg; options...)

Evaluate `𝑓̂(𝒑) = 𝑓(𝒔(𝒑), 𝒑)`, where `𝒔(𝒑)` is the steady-state solution
computed by the iterative solver `alg` (so that `𝑭(𝒔, 𝒑) = 0`). The
`cache` (built with `F1Cache`) is updated in place if `𝒑` has
changed since the last call.
"""
function objective(f, F, cache, p, alg; options...)
    update_solution!(F, cache, p, alg; options...)
    return f(cache.s, p)
end

"""
    gradient(f, F, ∇ₓf, cache, p, alg; options...)

Return the gradient `∇𝑓̂(𝒑)` as a `Vector` using the F-1 method. (Prior
to F1Method 0.6 this returned a `1 × m` row matrix; the new shape is
`Vector{T}` of length `m = length(p)`, which is what Optim.jl /
Optimization.jl expect.)
"""
function gradient(f, F, ∇ₓf, cache, p, alg; options...)
    update_cache!(F, ∇ₓf, cache, p, alg; options...)
    s, ∇s = cache.s, cache.∇s
    ∇ₚf = DI.gradient(p -> f(s, p), cache.ad, p)
    return vec(cache.∇ₓf * ∇s) + ∇ₚf
end

"""
    hessian(f, F, ∇ₓf, cache, p, alg; options...)

Return the Hessian `∇²𝑓̂(𝒑)` as an `m × m` `Matrix` using the F-1 method.
"""
function hessian(f, F, ∇ₓf, cache, p, alg; options...)
    update_cache!(F, ∇ₓf, cache, p, alg; options...)
    s, ∇s, m = cache.s, cache.∇s, length(p)
    A⁻ᵀ∇ₓfᵀ = vec(solveAᵀ!(cache, cache.∇ₓf'))   # independent of (𝑗, 𝑘)
    H(λ) = f(s + ∇s * λ, p + λ) - F(s + ∇s * λ, p + λ)' * A⁻ᵀ∇ₓfᵀ
    return DI.hessian(H, cache.ad, zeros(m))
end

"""
    F1Cache(F, ∇ₓf, x, p, alg; ad=AutoForwardDiff(), linsolve=nothing, options...)

Construct the F-1 memory cache. `ad` selects the AD backend used for
parameter-side derivatives (any `ADTypes.AbstractADType`); it defaults
to `AutoForwardDiff()`, matching pre-0.7 behaviour. Passing a
`linsolve <: SciMLBase.AbstractLinearAlgorithm` (from
[LinearSolve.jl](https://github.com/SciML/LinearSolve.jl)) routes the
internal linear solves through a `LinearSolve.LinearCache`, enabling
symbolic-factorization reuse for sparse direct solvers and buffer
reuse for dense; the default `linsolve = nothing` keeps the legacy
`factorize` / `\\` path. Remaining keyword arguments are forwarded to
`solve(::SteadyStateProblem, alg; ...)`.
"""
function F1Cache(F, ∇ₓf, x, p, alg;
                 ad::AbstractADType = AutoForwardDiff(),
                 linsolve = nothing,
                 options...)
    return _init_F1Cache(linsolve, F, ∇ₓf, x, p, alg, ad; options...)
end

function _init_F1Cache(::Nothing, F, ∇ₓf, x, p, alg, ad; options...)
    x = copy(x)
    p = copy(p)
    psol = copy(p)
    prob = SteadyStateProblem(F, x, p)
    s = solve(prob, alg; options...).u
    A = factorize(F.jac(s, p))
    ∇ₚF = DI.jacobian(p -> F(s, p), ad, p)
    return F1Cache(A, s, A \ -∇ₚF, ∇ₓf(s, p), p, psol, ad)
end

"""
    initialize_mem(args...; kwargs...)

Deprecated; use `F1Cache` instead. Retained for 0.7.x; will be
removed in 0.8.
"""
function initialize_mem(args...; kwargs...)
    Base.depwarn(
        "initialize_mem is deprecated; use F1Cache(...) instead",
        :initialize_mem,
    )
    return F1Cache(args...; kwargs...)
end

"""
    optimization_function(f, F, ∇ₓf, cache, alg; options...)

Wrap the F-1-method `objective` / `gradient` / `hessian` triple at the
F-1 `cache` into a `SciMLBase.OptimizationFunction`, ready to feed into
`OptimizationProblem` and any solver from the [Optimization.jl][1]
ecosystem.

This method lives in `F1MethodOptimizationExt`; it is only available
when `Optimization.jl` has been loaded.

[1]: https://docs.sciml.ai/Optimization/stable/
"""
function optimization_function end

# TODO(F1Method 0.8?): consider whether `optimization_function` should
# accept its own `linsolve` and rebuild the F1Cache internally, vs.
# always inheriting from the F1Cache the caller passes in. Today the
# latter is simpler — F1Cache is the carrier of every cached artifact
# including the linear solver — but if F1Cache construction grows more
# args this could become awkward.

end # module
