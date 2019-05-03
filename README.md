
<img src="/img/Logo_v2.png" alt="logo" title="F1method" align="right" height="200"/>

F-1 differentiation method
==========================

This package implements the F-1 method descibed in Pasquier et al. (2019).
It allows for efficient quasi-auto-differentiation of an objective function define implicitly by a the solution of a steady-state problem represented by discretized nonlinear partial differential equations taking the form

<img src="https://latex.codecogs.com/svg.latex?&space;\frac{\partial&space;\boldsymbol{x}}{\partial&space;t}&space;=&space;\boldsymbol{F}(\boldsymbol{x},&space;\boldsymbol{p})" title="Eq1"/>

where <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{x}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\boldsymbol{x}" title="\boldsymbol{x}" /></a> is a column vector of the model state variables and <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{p}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\boldsymbol{p}" title="\boldsymbol{p}" /></a> is a vector of parameters.
Specifically, the F-1 method allows for an efficient computation of both the gradient and the Hessian matrix of an objective function <a href="https://www.codecogs.com/eqnedit.php?latex=\hat{f}(\boldsymbol{p})&space;\equiv&space;f(\boldsymbol{s},\boldsymbol{p})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\hat{f}(\boldsymbol{p})&space;\equiv&space;f(\boldsymbol{s},\boldsymbol{p})" title="\hat{f}(\boldsymbol{p}) \equiv f(\boldsymbol{s},\boldsymbol{p})" /></a> where <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{s}(\boldsymbol{p})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\boldsymbol{s}(\boldsymbol{p})" title="\boldsymbol{s}(\boldsymbol{p})" /></a> is the steady-state solution of the system, i.e., such that  <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{F}(\boldsymbol{s},&space;\boldsymbol{p})&space;=&space;0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\boldsymbol{F}(\boldsymbol{s},&space;\boldsymbol{p})&space;=&space;0" title="\boldsymbol{F}(\boldsymbol{s}, \boldsymbol{p}) = 0" /></a>.

A requirement of the F-1 method is that the Jacobian matrix <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{A}&space;=&space;\nabla_{\boldsymbol{x}}\boldsymbol{F}(\boldsymbol{x},\boldsymbol{p})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{A}&space;=&space;\nabla_{\boldsymbol{x}}\boldsymbol{F}(\boldsymbol{x},\boldsymbol{p})" title="\mathbf{A} = \nabla_{\boldsymbol{x}}\boldsymbol{F}(\boldsymbol{x},\boldsymbol{p})" /></a> can be created, stored, and factorized.

To use the F-1 method, the user must:

- Make sure that there is a suitable algorithm `alg` to solve the steady-state equation, that is implemented with the `solve` function of [DiffEqBase](https://github.com/JuliaDiffEq/DiffEqBase.jl). (An example is given in the CI tests — see the `test/` directory.)
- Provide the derivatives of `f` and `F` with respect to the state, `x`.

### Simple usage

Make sure you have olverloaded `solve` from DiffEqBase.
Once an initial state, `x₀`, and some parameters, `p₀`, are chosen, simply evaluate the derivatives with

```julia
# Initialize the cache for storing reusable objects
mem = initialize_mem(x₀, p₀)

# Compute the objective function, 𝑓̂(𝒑)
f̂(p) = F1.f̂(f, F, ∇ₓF, mem, p, myAlg(); my_options...)
f̂(p₀)

# Compute the gradient, ∇𝑓̂(𝒑)
∇f̂(p) = F1.∇f̂(f, F, ∇ₓf, ∇ₓF, mem, p₀, myAlg(); my_options...)
∇f̂(p₀)

# Compute the Hessian matrix, ∇𝑓̂(𝒑)
∇²f̂(p) = F1.∇²f̂(f, F, ∇ₓf, ∇ₓF, mem, p₀, myAlg(); my_options...)
∇²f̂(p₀)
```
