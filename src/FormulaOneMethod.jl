module FormulaOneMethod

using LinearAlgebra
using DualNumbers, HyperDualNumbers, DiffEqBase

mutable struct Buffer
    p   # p
    s   # s(p)
    A   # factors of ∇ₓF(s, p)
    ∇s  # ∇s(p)
    ∇ₓf # ∇ₓf(s, p)
end

function update_buffer!(f, F, ∇ₓf, ∇ₓF, buffer, p, alg; options...)
    if p ≠ buffer.p       # only update if p has changed
        s, m = buffer.s, length(p)
        prob = SteadyStateProblem(F, ∇ₓF, s, p) # define problem
        buffer.s .= solve(prob, alg, buffer.A; options...) # update s (inner solver)
        ∇ₚF = hcat([𝔇(F(s, p + ε * e(j,m))) for j in 1:m]...) # Eq.(?)
        buffer.A = factorize(∇ₓF(s,p))  # update factors of ∇ₓF(s, p)
        buffer.∇s .= buffer.A \ -∇ₚF    # update ∇s via Eq.(?)
        buffer.∇ₓf .= ∇ₓf(s,p)          # update ∇ₓf(s, p)
        buffer.p = p      # update p
    end
end

function f̂(f, F, ∇ₓf, ∇ₓF, buffer, p, alg; options...) # objective
    update_buffer!(f, F, ∇ₓf, ∇ₓF, buffer, p, alg; options...)
    s = buffer.s
    return f(s,p)
end

function ∇f̂(f, F, ∇ₓf, ∇ₓF, buffer, p, alg; options...) # gradient
    update_buffer!(f, F, ∇ₓf, ∇ₓF, buffer, p, alg; options...)
    s, ∇s, m = buffer.s, buffer.∇s, length(p)
    ∇ₚf = [𝔇(f(s,p + ε * e(j,m))) for j in 1:m]' # Eq. (?)
    return buffer.∇ₓf * ∇s + ∇ₚf               # Eq. (?)
end

function ∇²f̂(f, F, ∇ₓf, ∇ₓF, buffer, p, alg; options...) # Hessian
    update_buffer!(f, F, ∇ₓf, ∇ₓF, buffer, p, alg; options...)
    s, A, ∇s, m = buffer.s, buffer.A, buffer.∇s, length(p)
    A⁻ᵀ∇ₓfᵀ = vec(A' \ buffer.∇ₓf') # independent of (j,k)
    out = zeros(m,m)       # preallocate
    for j in 1:m, k in j:m # Loop for Eq.(?)
        pⱼₖ = p + ε₁ * e(j,m) + ε₂ * e(k,m)           # Hyperdual p
        xⱼₖ = s + ε₁ * ∇s * e(j,m) + ε₂ * ∇s * e(k,m) # Hyperdual x
        out[j,k] = ℌ(f(xⱼₖ,pⱼₖ)) - ℌ(F(xⱼₖ,pⱼₖ))' * A⁻ᵀ∇ₓfᵀ # Eq.(?)
        j ≠ k ? out[k,j] = out[j,k] : nothing   # symmetry
    end
    return out
end

function initialize_buffer(f, F, ∇ₓf, ∇ₓF, x₀, p, alg; options...)
    m = length(p)
    prob = SteadyStateProblem(F, ∇ₓF, x₀, p)
    s = solve(prob, alg; options...).u
    ∇ₚF = hcat([𝔇(F(s, p + ε * e(j,m))) for j in 1:m]...)
    A = factorize(∇ₓF(s,p))
    ∇s = A \ -∇ₚF
    return Buffer(p, s, A, ∇s, ∇ₓf(s,p))
end

# Helper functions
e(j, m) = [i == j for i in 1:m]      # j-th basis vector
𝔇(x) = DualNumbers.dualpart.(x)      # dual part
ℌ(x) = HyperDualNumbers.ε₁ε₂part.(x) # hyperdual part

export f̂, ∇f̂, ∇²f̂, initialize_buffer

end # module
