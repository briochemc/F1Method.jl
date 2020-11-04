
# This AIBECS test is derived from the P-model tutorial in AIBECS.jl
# But it uses one of the samll circulations to test derivatives of the cost function


# AIBECS model
grd, T_Circ = Primeau_2x2x2.load()
T_DIP(p) = T_Circ
T_POP(p) = transportoperator(grd, z -> w(z,p))
function w(z,p)
    @unpack w₀, w′ = p
    return @. w₀ + w′ * z
end
const z_top = topdepthvec(grd) # uptake only in top layer
function U(x,p)
    @unpack σ, τ_DIP, k = p
    return @. σ * x/τ_DIP * x/(x+k) * (z_top==0) * (x≥0)
end
function R(x,p)
    @unpack τ_POP = p
    return x / τ_POP
end
function G_DIP(DIP, POP, p)
    @unpack DIP_geo, τ_geo = p
    return @. -$U(DIP,p) + $R(POP,p) + (DIP_geo - DIP) / τ_geo
end
function G_POP(DIP, POP, p)
    @unpack τ_geo = p
    return @. $U(DIP,p) - $R(POP,p) - POP / τ_geo
end
import AIBECS: @units, units
import AIBECS: @limits, limits
import AIBECS: @initial_value, initial_value
import AIBECS: @flattenable, flattenable, flatten
const ∞ = Inf
@initial_value @units @flattenable @limits struct PmodelParameters{U} <: AbstractParameters{U}
    w₀::U       |  0.64 | m/d      | true | (0,∞)
    w′::U       |  0.13 | m/d/m    | true | (0,∞)
    τ_DIP::U    | 230.0 | d        | true | (0,∞)
    k::U        |  6.62 | μmol/m^3 | true | (0,∞)
    τ_POP::U    |   5.0 | d        | true | (0,∞)
    τ_geo::U    |   1.0 | Myr      | true | (0,∞)
    DIP_geo::U  |  2.12 | mmol/m^3 | true | (-∞,∞)
    σ::U        |  0.3  | NoUnits  | true | (0,1)
end
import AIBECS: @prior, prior
function prior(::Type{T}, s::Symbol) where {T<:AbstractParameters}
    if flattenable(T, s)
        U = units(T, s)
        if limits(T, s) == (0,∞)
            μ = log(ustrip(upreferred(initial_value(T, s) * U)))
            return LogNormal(μ, 1.0)
        elseif limits(T, s) == (-∞,∞)
            μ = ustrip(upreferred(initial_value(T, s) * U))
            σ = ustrip(upreferred(10.0U)) # Assumes that a sensible unit is chosen!
            return Normal(μ, σ)
        elseif limits(T, s) == (0,1)
            return Uniform(0,1)
        end
    else
        return nothing
    end
end
prior(::T, s::Symbol) where {T<:AbstractParameters} = prior(T,s)
prior(::Type{T}) where {T<:AbstractParameters} = Tuple(prior(T,s) for s in AIBECS.symbols(T))
prior(::T) where {T<:AbstractParameters} = prior(T)
p = PmodelParameters()
nb = sum(iswet(grd))
𝐹, ∇ₓ𝐹 = state_function_and_Jacobian((T_DIP, T_POP), (G_DIP, G_POP), nb)
@unpack DIP_geo = p
x = DIP_geo * ones(2nb) # initial guess
prob = SteadyStateProblem(𝐹, ∇ₓ𝐹, x, p)
sol = solve(prob, CTKAlg()).u

# Observations from World Ocean Atlas used to evaluate the
# AIBECS model mismatch with observations
# and generate the objective function
ρSW = 1.035u"kg/L" # approximate mean sea water density
DIPobs = ustrip(upreferred(WorldOceanAtlasTools.observations("PO₄") * ρSW))
modify(DIP, POP) = (DIP,)
ωs = (1.0,) # the weight for the mismatch (weight of POP = 0)
ωp = 1e-4       # the weight for the parameters prior estimates
obs = (DIPobs,)
𝑓, ∇ₓ𝑓, ∇ₚ𝑓 = generate_objective_and_derivatives(ωs, ωp, grd, modify, obs)

# Now we apply the F1 method
mem = F1Method.initialize_mem(x, p)
objective(p) = F1Method.objective(𝑓, 𝐹, ∇ₓ𝐹, mem, p, CTKAlg(), preprint="obj ", τstop=ustrip(u"s", 1e3u"Myr"))
gradient(p) = F1Method.gradient(𝑓, 𝐹, ∇ₓ𝑓, ∇ₓ𝐹, mem, p, CTKAlg(), preprint="grad", τstop=ustrip(u"s", 1e3u"Myr"))
hessian(p) = F1Method.hessian(𝑓, 𝐹, ∇ₓ𝑓, ∇ₓ𝐹, mem, p, CTKAlg(), preprint="hess ", τstop=ustrip(u"s", 1e3u"Myr"))

# and convert p::PmodelParameters to λ::Vector according to AIBECS change of variables
λ2p = subfun(typeof(p))
∇λ2p = ∇subfun(typeof(p))
∇²λ2p = ∇²subfun(typeof(p))
p2λ = invsubfun(typeof(p))
λ = p2λ(p)
function obj(λ)
    show(λ2p(λ))
    return objective(λ2p(λ))
end
function grad(λ)
    return gradient(λ2p(λ)) * Diagonal(∇λ2p(λ))
end
function hess(λ)
    ∇p = Diagonal(∇λ2p(λ)) # for variable change
    ∇²p = Diagonal(∇²λ2p(λ)) # for variable change
    G = vec(gradient(λ2p(λ)))
    H = hessian(λ2p(λ))
    return ∇p * H * ∇p + Diagonal(G) * ∇²p
end

# Finally we test the result with the "reliable" FiniteDiff :)
λ = p2λ(p)
@test FiniteDiff.finite_difference_gradient(obj, λ)' ≈ grad(λ) rtol=1e-3
@test FiniteDiff.finite_difference_hessian(obj, λ) ≈ hess(λ) rtol=1e-3






