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
    Mem

Memory cache holding the steady-state solution and the linear-algebra
artifacts that the F-1 method reuses across `objective`, `gradient`, and
`hessian` calls at the same parameter point.

# Fields
- `s`    : the steady-state solution, `𝒔(𝒑)`
- `A`    : the factors of `𝐀 = ∇ₓ𝑭(𝒔, 𝒑)`
- `∇s`   : the parameter-Jacobian of the steady state, `∇𝒔(𝒑)`
- `∇ₓf`  : the state-Jacobian of the objective at the steady state,
           `∇ₓ𝑓(𝒔, 𝒑)`
- `p`    : the parameter vector that matches `A`, `∇s`, `∇ₓf` (and `s`)
- `psol` : the parameter vector that matches `s` (may be ahead of `p`
           when only the steady-state solve has been refreshed)
- `ad`   : the AD backend used for parameter-side derivatives

Initialise with [`initialize_mem`](@ref).
"""
mutable struct Mem{Ts, TA, T∇s, T∇ₓf, Tp, AD <: AbstractADType}
    s::Ts
    A::TA
    ∇s::T∇s
    ∇ₓf::T∇ₓf
    p::Tp
    psol::Tp
    ad::AD
end


function update_mem!(F, ∇ₓf, mem, p, alg; options...)
    if p ≠ mem.p                              # only update if 𝒑 has changed
        update_solution!(F, mem, p, alg; options...)
        ∇ₚF = DI.jacobian(p -> F(mem.s, p), mem.ad, p)
        mem.A = factorize(F.jac(mem.s, p))    # update factors of ∇ₓ𝑭(𝒔, 𝒑)
        mem.∇s .= mem.A \ -∇ₚF                # update ∇𝒔
        mem.∇ₓf .= ∇ₓf(mem.s, p)              # update ∇ₓ𝑓(𝒔, 𝒑)
        mem.p .= p                            # update 𝒑 for the cached values
    end
end

function update_solution!(F, mem, p, alg; options...)
    if p ≠ mem.psol
        prob = SteadyStateProblem(F, mem.s, p)
        mem.s .= solve(prob, alg; options...).u
        mem.psol .= p
    end
end

"""
    objective(f, F, mem, p, alg; options...)

Evaluate `𝑓̂(𝒑) = 𝑓(𝒔(𝒑), 𝒑)`, where `𝒔(𝒑)` is the steady-state solution
computed by the iterative solver `alg` (so that `𝑭(𝒔, 𝒑) = 0`). The
memory cache `mem` (built with [`initialize_mem`](@ref)) is updated
in place if `𝒑` has changed since the last call.
"""
function objective(f, F, mem, p, alg; options...)
    update_solution!(F, mem, p, alg; options...)
    return f(mem.s, p)
end

"""
    gradient(f, F, ∇ₓf, mem, p, alg; options...)

Return the gradient `∇𝑓̂(𝒑)` as a `Vector` using the F-1 method. (Prior
to F1Method 0.6 this returned a `1 × m` row matrix; the new shape is
`Vector{T}` of length `m = length(p)`, which is what Optim.jl /
Optimization.jl expect.)
"""
function gradient(f, F, ∇ₓf, mem, p, alg; options...)
    update_mem!(F, ∇ₓf, mem, p, alg; options...)
    s, ∇s = mem.s, mem.∇s
    ∇ₚf = DI.gradient(p -> f(s, p), mem.ad, p)
    return vec(mem.∇ₓf * ∇s) + ∇ₚf
end

"""
    hessian(f, F, ∇ₓf, mem, p, alg; options...)

Return the Hessian `∇²𝑓̂(𝒑)` as an `m × m` `Matrix` using the F-1 method.
"""
function hessian(f, F, ∇ₓf, mem, p, alg; options...)
    update_mem!(F, ∇ₓf, mem, p, alg; options...)
    s, A, ∇s, m = mem.s, mem.A, mem.∇s, length(p)
    A⁻ᵀ∇ₓfᵀ = vec(A' \ mem.∇ₓf')              # independent of (𝑗, 𝑘)
    H(λ) = f(s + ∇s * λ, p + λ) - F(s + ∇s * λ, p + λ)' * A⁻ᵀ∇ₓfᵀ
    return DI.hessian(H, mem.ad, zeros(m))
end

"""
    initialize_mem(F, ∇ₓf, x, p, alg; ad=AutoForwardDiff(), options...)

Initialise the F-1 memory cache. `ad` selects the AD backend used for
parameter-side derivatives (any `ADTypes.AbstractADType`); it defaults
to `AutoForwardDiff()`, matching pre-0.7 behaviour. Remaining keyword
arguments are forwarded to `solve(::SteadyStateProblem, alg; ...)`.
"""
function initialize_mem(F, ∇ₓf, x, p, alg;
                        ad::AbstractADType = AutoForwardDiff(), options...)
    x = copy(x)
    p = copy(p)
    psol = copy(p)
    prob = SteadyStateProblem(F, x, p)
    s = solve(prob, alg; options...).u
    A = factorize(F.jac(s, p))
    ∇ₚF = DI.jacobian(p -> F(s, p), ad, p)
    return Mem(s, A, A \ -∇ₚF, ∇ₓf(s, p), p, psol, ad)
end

"""
    optimization_function(f, F, ∇ₓf, mem, alg; options...)

Wrap the F-1-method `objective` / `gradient` / `hessian` triple at memory
cache `mem` into a `SciMLBase.OptimizationFunction`, ready to feed into
`OptimizationProblem` and any solver from the [Optimization.jl][1]
ecosystem.

This method lives in `F1MethodOptimizationExt`; it is only available
when `Optimization.jl` has been loaded.

[1]: https://docs.sciml.ai/Optimization/stable/
"""
function optimization_function end

end # module
