module F1Method

#======================================================================
This package (the code below) implements the F-1 method as described
in the work of Pasquier et al. (2019). The numbers in parentheses
refer to the Equation numbers in the above manuscript. A bibtex
citation file is available in the GitHub repository.
======================================================================#

using LinearAlgebra, ForwardDiff, DiffEqBase

"""
    Mem

Memory cache to store reusable objects.
Contains
- `s`   the steady-state solution, 𝒔(𝒑)
- `A`   the factors of 𝐀 = ∇ₓ𝑭(𝒔,𝒑)
- `∇s`  the derivative ∇𝒔(𝒑)
- `∇ₓf` the derivative ∇ₓ𝑓(𝒔,𝒑)
- `p`   the parameters 𝒑
The `Mem`-type object should be initialized with `initialize_mem`.
"""
mutable struct Mem
    s     # 𝒔(𝒑)
    A     # factors of 𝐀 = ∇ₓ𝑭(𝒔,𝒑)
    ∇s    # ∇𝒔(𝒑)
    ∇ₓf   # ∇ₓ𝑓(𝒔,𝒑)
    p     # 𝒑
end


function update_mem!(f, F, ∇ₓf, ∇ₓF, mem, p, alg; options...)
    if p ≠ mem.p                      # only update mem if 𝒑 has changed
        update_solution!(F, ∇ₓF, mem, p, alg; options...)
        s, m = mem.s.u, length(p)
        ∇ₚF = ForwardDiff.jacobian(λ -> F(s,p+λ), zeros(m))
        mem.A = factorize(∇ₓF(s,p))   # update factors of ∇ₓ𝑭(𝒔,𝒑)
        mem.∇s .= mem.A \ -∇ₚF        # update ∇𝒔
        mem.∇ₓf .= ∇ₓf(s,p)           # update ∇ₓ𝑓(𝒔,𝒑)
        mem.p = p                     # update 𝒑
    end
end

function update_solution!(F, ∇ₓF, mem, p, alg; options...)
    if ~(mem.s isa SteadyStateSolution) || p ≠ mem.s.prob.p
        mem.s isa SteadyStateSolution ? x = mem.s.u : x = mem.s
        prob = SteadyStateProblem(F, ∇ₓF, x, p)       # define problem
        mem.s = solve(prob, alg; options...)          # update 𝒔
    end
end

"""
    objective(f, F, ∇ₓF, mem, p, alg; options...)

Returns `f(x,p)` such that `F(x,p)=0` using the F-1 method.

Specifically, `objective(f, F, ∇ₓF, mem, p, alg; options...)`
evaluates the objective function defined by `f̂(p) = f(s(p),p)`, where
`s(p)`, which is the steady-state solution (i.e., such that `F(s(p),p)=0`)
is computed by the iterative Newton-type solver `alg`.
The Jacobian, `∇ₓF`, and the memory cache `mem` must be supplied.
"""
function objective(f, F, ∇ₓF, mem, p, alg; options...)
    update_solution!(F, ∇ₓF, mem, p, alg; options...)
    return f(mem.s,p)
end

"""
    gradient(f, F, ∇ₓf, ∇ₓF, mem, p, alg; options...)

Returns the gradient of the `objective` function using the F-1 method.
"""
function gradient(f, F, ∇ₓf, ∇ₓF, mem, p, alg; options...)
    update_mem!(f, F, ∇ₓf, ∇ₓF, mem, p, alg; options...)
    s, ∇s, m = mem.s, mem.∇s, length(p)
    ∇ₚf = ForwardDiff.jacobian(λ -> [f(s,p+λ)], zeros(m))
    return mem.∇ₓf * ∇s + ∇ₚf
end

"""
    hessian(f, F, ∇ₓf, ∇ₓF, mem, p, alg; options...)

Returns the Hessian of the `objective` function using the F-1 method.
"""
function hessian(f, F, ∇ₓf, ∇ₓF, mem, p, alg; options...)
    update_mem!(f, F, ∇ₓf, ∇ₓF, mem, p, alg; options...)
    s, A, ∇s, m = mem.s, mem.A, mem.∇s, length(p)
    A⁻ᵀ∇ₓfᵀ = vec(A' \ mem.∇ₓf') # independent of (𝑗,𝑘)
    H(λ) = f(s+∇s*λ, p+λ) - F(s+∇s*λ, p+λ)' * A⁻ᵀ∇ₓfᵀ
    return ForwardDiff.hessian(H, zeros(m))
end

"""
    initialize_mem(x, p)

Initializes the memory cache for the F-1 method.
"""
function initialize_mem(x, p)
    n, m = length(x), length(p)
    return Mem(copy(x), nothing, zeros(n,m), zeros(1,n), nothing)
end

end
