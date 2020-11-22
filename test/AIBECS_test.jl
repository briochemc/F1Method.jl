
# This AIBECS test is derived from the P-model tutorial in AIBECS.jl
# But it uses one of the samll circulations to test derivatives of the cost function
module SubModuleAIBECS

using Test
using F1Method
using AIBECS
import AIBECS: @units, units
import AIBECS: @limits, limits
import AIBECS: @initial_value, initial_value
import AIBECS: @flattenable, flattenable, flatten
import AIBECS: @prior, prior
using Unitful: m, d, s, yr, Myr, mol, mmol, μmol, μM
using Distributions
using WorldOceanAtlasTools
using FiniteDiff

# AIBECS model
grd, T_Circ = Primeau_2x2x2.load()
T_DIP(p) = T_Circ
T_POP(p) = transportoperator(grd, z -> w(z,p))
function w(z,p)
    @unpack w₀, w′ = p
    w₀ + w′ * z
end
z_top = topdepthvec(grd) # uptake only in top layer
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

∞ = Inf
@initial_value @units @flattenable @limits struct PmodelParameters{U} <: AbstractParameters{U}
    w₀::U       |  0.64 | m/d      | true  | (0,∞)
    w′::U       |  0.13 | m/d/m    | true  | (0,∞)
    τ_DIP::U    | 230.0 | d        | true  | (0,∞)
    k::U        |  6.62 | μmol/m^3 | true  | (0,∞)
    τ_POP::U    |   5.0 | d        | true  | (0,∞)
    τ_geo::U    |   1.0 | Myr      | false | (0,∞)
    DIP_geo::U  |  2.12 | mmol/m^3 | true  | (-∞,∞)
    σ::U        |  0.3  | NoUnits  | true  | (0,1)
end
function prior(::Type{T}, s::Symbol) where {T<:AbstractParameters}
    if flattenable(T, s)
        lb, ub = limits(T, s)
        if (lb, ub) == (0,∞)
            μ = log(initial_value(T, s))
            return LogNormal(μ, 1.0)
        elseif (lb, ub) == (-∞,∞)
            μ = initial_value(T, s)
            σ = 10.0 # Assumes that a sensible unit is chosen (i.e., that within 10.0 * U)
            return Normal(μ, σ)
        else
            return LocationScale(lb, ub-lb, LogitNormal()) # <- The LogitNormal works well for Optim?
        end
    else
        return nothing
    end
end
prior(::T, s::Symbol) where {T<:AbstractParameters} = prior(T,s)
prior(::Type{T}) where {T<:AbstractParameters} = Tuple(prior(T,s) for s in AIBECS.symbols(T))
prior(::T) where {T<:AbstractParameters} = prior(T)
p = PmodelParameters()
λ = p2λ(p)
nb = sum(iswet(grd))
F, ∇ₓF = F_and_∇ₓF((T_DIP, T_POP), (G_DIP, G_POP), nb, PmodelParameters)
@unpack DIP_geo = p
x = DIP_geo * ones(2nb) # initial guess
prob = SteadyStateProblem(F, ∇ₓF, x, p)
sol = solve(prob, CTKAlg())

# Observations from World Ocean Atlas used to evaluate the
# AIBECS model mismatch with observations
# and generate the objective function
ρSW = 1.035u"kg/L" # approximate mean sea water density
DIPobs = ustrip(upreferred(WorldOceanAtlasTools.observations("PO₄") * ρSW))
modify(DIP, POP) = (DIP,)
ωs = (1.0,) # the weight for the mismatch (weight of POP = 0)
ωp = 1.0       # the weight for the parameters prior estimates
obs = (DIPobs,)
f, ∇ₓf = f_and_∇ₓf(ωs, ωp, grd, modify, obs, PmodelParameters)

# Now we apply the F1 method
τ = ustrip(u"s", 1e3u"Myr")
mem = F1Method.initialize_mem(F, ∇ₓf, ∇ₓF, x, λ, CTKAlg(), τstop=τ)
objective(λ) = F1Method.objective(f, F, ∇ₓF, mem, λ, CTKAlg(), τstop=τ)
gradient(λ) = F1Method.gradient(f, F, ∇ₓf, ∇ₓF, mem, λ, CTKAlg(), τstop=τ)
hessian(λ) = F1Method.hessian(f, F, ∇ₓf, ∇ₓF, mem, λ, CTKAlg(), τstop=τ)

# Finally we test the result with the "reliable" FiniteDiff :)
@test FiniteDiff.finite_difference_gradient(objective, 2λ) ≈ gradient(2λ)' rtol=1e-3
@test FiniteDiff.finite_difference_hessian(objective, 2λ) ≈ hessian(2λ) rtol=1e-3
@test FiniteDiff.finite_difference_gradient(objective, λ) ≈ gradient(λ)' rtol=1e-3
@test FiniteDiff.finite_difference_hessian(objective, λ) ≈ hessian(λ) rtol=1e-3




end #submodule

