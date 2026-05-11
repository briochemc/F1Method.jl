module F1MethodOptimizationExt

using F1Method
using Optimization: OptimizationFunction

function F1Method.optimization_function(f, F, ∇ₓf, mem, alg; options...)
    obj  = (λ, _) -> F1Method.objective(f, F, mem, λ, alg; options...)
    grad! = (g, λ, _) -> (g .= F1Method.gradient(f, F, ∇ₓf, mem, λ, alg; options...); nothing)
    hess! = (H, λ, _) -> (H .= F1Method.hessian(f, F, ∇ₓf, mem, λ, alg; options...); nothing)
    return OptimizationFunction(obj; grad = grad!, hess = hess!)
end

end # module
