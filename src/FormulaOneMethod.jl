module FormulaOneMethod

using DualNumbers, HyperDualNumbers

mutable struct Buffer
    p   # p
    s   # s(p)
    A   # factors of ∇ₓF(s, p)
    ∇s  # ∇s(p)
    ∇ₓf # ∇ₓf(s, p)
end

function update_buffer!(f, F, ∇ₓf, ∇ₓF, buffer, p; options...)
    if p ≠ buffer.p       # only update if p has changed
        buffer.p = p      # update p
        s, A, ∇s = buffer.s, buffer.A, buffer.∇s  # unpack buffer
        prob = steady_state_problem(F, ∇ₓF, s, p) # define problem
        s .= solve(prob, options...) # update s (inner solver)
        ∇ₚF = hcat([𝔇(F(s, p + ε * e(j))) for j in 1:m]) # Eq.(?)
        A .= factorize(∇ₓF(s, p))      # update factors of ∇ₓF(s, p)
        ∇s .= A \ -∇ₚF                 # update ∇s via Eq.(?)
        buffer.∇ₓf .= ∇ₓf(s, p)        # update ∇ₓf(s, p)
    end
end

function f̂(f, F, ∇ₓf, ∇ₓF, buffer, p; options...) # objective
    update_buffer!(f, F, ∇ₓf, ∇ₓF, buffer, p; options...)
    s = buffer.s
    return f(s, p)
end

function ∇f̂(f, F, ∇ₓf, ∇ₓF, buffer, p; options...) # gradient
    update_buffer!(f, F, ∇ₓf, ∇ₓF, buffer, p; options...)
    s, ∇s = buffer.s, buffer.∇s
    ∇ₚf = [𝔇(f(s, p + ε * e(j))) for j in 1:m] # Eq. (?)
    return buffer.∇ₓf * ∇s + ∇ₚf               # Eq. (?)
end

function ∇²f̂(f, F, ∇ₓf, ∇ₓF, buffer, p; options...) # Hessian
    update_buffer!(f, F, ∇ₓf, ∇ₓF, buffer, p; options...)
    s, A, ∇s = buffer.s, buffer.A, buffer.∇s
    A⁻ᵀ∇ₓfᵀ = vec(A' \ buffer.∇ₓf') # independent of (j,k)
    out = zeros(m, m)      # preallocate
    for j in 1:m, k in j:m # Loop for Eq.(?)
        pⱼₖ = p + ε₁ * e(j) + ε₂ * e(k)           # Hyperdual p
        xⱼₖ = s + ε₁ * ∇s * e(j) + ε₂ * ∇s * e(k) # Hyperdual x
        out[j, k] = ℌ(f(xⱼₖ, pⱼₖ)) - ℌ(F(xⱼₖ, pⱼₖ))' * A⁻ᵀ∇ₓfᵀ # Eq.(?)
        j ≠ k ? out[k, j] = out[j, k] : nothing   # symmetry
    end
    return out
end

# Helper functions
e(j) = [i == j for i in 1:m]         # j-th basis vector
𝔇(x) = DualNumbers.dualpart.(x)      # dual part
ℌ(x) = HyperDualNumbers.ε₁ε₂part.(x) # hyperdual part

export f̂, ∇f̂, ∇²f̂

end # module
