
<img src="/img/Logo_v2.png" alt="logo" title="F1method" align="right" height="200"/>

F-1 method
==========

This package implements the efficient computation derived in (Citation).

This method applies to models that can be represented by a system of discrete nonlinear partial differential equations taking the form

<img src="https://latex.codecogs.com/svg.latex?&space;\frac{\partial&space;\boldsymbol{x}}{\partial&space;t}&space;=&space;\boldsymbol{F}(\boldsymbol{x},&space;\boldsymbol{p})" title="Eq1"/>

where <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{x}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\boldsymbol{x}" title="\boldsymbol{x}" /></a> is a column vector of the model state variables and <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{p}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\boldsymbol{p}" title="\boldsymbol{p}" /></a> is a vector of parameters.
Specifically, if one aims to minimize an objective function <a href="https://www.codecogs.com/eqnedit.php?latex=\hat{f}(\boldsymbol{p})&space;\equiv&space;f(\boldsymbol{s},\boldsymbol{p})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\hat{f}(\boldsymbol{p})&space;\equiv&space;f(\boldsymbol{s},\boldsymbol{p})" title="\hat{f}(\boldsymbol{p}) \equiv f(\boldsymbol{s},\boldsymbol{p})" /></a> where <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{s}(\boldsymbol{p})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\boldsymbol{s}(\boldsymbol{p})" title="\boldsymbol{s}(\boldsymbol{p})" /></a> is the steady-state solution of the system, i.e., such that  <a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{F}(\boldsymbol{s},&space;\boldsymbol{p})&space;=&space;0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\boldsymbol{F}(\boldsymbol{s},&space;\boldsymbol{p})&space;=&space;0" title="\boldsymbol{F}(\boldsymbol{s}, \boldsymbol{p}) = 0" /></a>, the F-1 method allows for an efficient computation of the Hessian matrix.
(As long as the Jacobian matrix <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{A}&space;=&space;\nabla_{\boldsymbol{x}}\boldsymbol{F}(\boldsymbol{x},\boldsymbol{p})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{A}&space;=&space;\nabla_{\boldsymbol{x}}\boldsymbol{F}(\boldsymbol{x},\boldsymbol{p})" title="\mathbf{A} = \nabla_{\boldsymbol{x}}\boldsymbol{F}(\boldsymbol{x},\boldsymbol{p})" /></a> can be created, factorized, and stored.)
