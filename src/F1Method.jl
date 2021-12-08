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
mutable struct Mem{Ts, TA, T∇s, T∇ₓf, Tp}
    s::Ts     # 𝒔(𝒑)
    A::TA     # factors of 𝐀 = ∇ₓ𝑭(𝒔,𝒑)
    ∇s::T∇s   # ∇𝒔(𝒑)
    ∇ₓf::T∇ₓf # ∇ₓ𝑓(𝒔,𝒑)
    p::Tp     # 𝒑 that matches all but 𝒔(𝒑)
    psol::Tp  # 𝒑 that matches 𝒔(𝒑)
end


function update_mem!(f, F, ∇ₓf, mem, p, alg; options...)
    if p ≠ mem.p                      # only update mem if 𝒑 has changed
        update_solution!(F, mem, p, alg; options...)
        ∇ₚF = ForwardDiff.jacobian(p -> F(mem.s, p), p)
        mem.A = factorize(F.jac(mem.s, p)) # update factors of ∇ₓ𝑭(𝒔,𝒑)
        mem.∇s .= mem.A \ -∇ₚF           # update ∇𝒔
        mem.∇ₓf .= ∇ₓf(mem.s, p)         # update ∇ₓ𝑓(𝒔,𝒑)
        mem.p .= p                  # update 𝒑 for the variables above
    end
end

function update_solution!(F, mem, p, alg; options...)
    if p ≠ mem.psol
        prob = SteadyStateProblem(F, mem.s, p) # define problem
        mem.s .= solve(prob, alg; options...).u      # update 𝒔
        mem.psol .= p                          # update 𝒑 for 𝒔
    end
end

"""
    objective(f, F, mem, p, alg; options...)

Returns `f(x,p)` such that `F(x,p)=0` using the F-1 method.

Specifically, `objective(f, F, mem, p, alg; options...)`
evaluates the objective function defined by `f̂(p) = f(s(p),p)`, where
`s(p)`, which is the steady-state solution (i.e., such that `F(s(p),p)=0`)
is computed by the iterative Newton-type solver `alg`.
The memory cache `mem` must be supplied.
"""
function objective(f, F, mem, p, alg; options...)
    update_solution!(F, mem, p, alg; options...)
    return f(mem.s, p)
end

"""
    gradient(f, F, ∇ₓf, mem, p, alg; options...)

Returns the gradient of the `objective` function using the F-1 method.
"""
function gradient(f, F, ∇ₓf, mem, p, alg; options...)
    update_mem!(f, F, ∇ₓf, mem, p, alg; options...)
    s, ∇s, m = mem.s, mem.∇s, length(p)
    ∇ₚf = ForwardDiff.jacobian(p -> [f(s,p)], p)
    return mem.∇ₓf * ∇s + ∇ₚf
end

"""
    hessian(f, F, ∇ₓf, mem, p, alg; options...)

Returns the Hessian of the `objective` function using the F-1 method.
"""
function hessian(f, F, ∇ₓf, mem, p, alg; options...)
    update_mem!(f, F, ∇ₓf, mem, p, alg; options...)
    s, A, ∇s, m = mem.s, mem.A, mem.∇s, length(p)
    A⁻ᵀ∇ₓfᵀ = vec(A' \ mem.∇ₓf') # independent of (𝑗,𝑘)
    H(λ) = f(s+∇s*λ, p+λ) - F(s+∇s*λ, p+λ)' * A⁻ᵀ∇ₓfᵀ
    return ForwardDiff.hessian(H, zeros(m))
end

"""
    initialize_mem(x, p)

Initializes the memory cache for the F-1 method.
"""
function initialize_mem(F, ∇ₓf, x, p, alg; options...)
    x = copy(x)
    p = copy(p)
    psol = copy(p)
    prob = SteadyStateProblem(F, x, p)
    s = solve(prob, alg; options...).u
    A = factorize(F.jac(s, p))
    ∇ₚF = ForwardDiff.jacobian(p -> F(s, p), p)
    return Mem(s, A, A \ -∇ₚF, ∇ₓf(s, p), p, psol)
end

end
