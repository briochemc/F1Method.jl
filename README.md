
<img src="/img/Logo_v2.pdf" alt="logo" title="F1method" align="right" height="200"/>

F-1 method
==========

This package implements the efficient computation derived in (Citation).

This method applies to models that can be represented by a system of discrete nonlinear partial differential equations taking the form

<img src="https://latex.codecogs.com/svg.latex?&space;\frac{\partial&space;\boldsymbol{x}}{\partial&space;t}&space;=&space;\boldsymbol{F}(\boldsymbol{x},&space;\boldsymbol{p})" title="Eq1"/>

where <img src="https://latex.codecogs.com/svg.latex?\inline&space;\large&space;\boldsymbol{x}" title="\boldsymbol{x}" /> is a column vector of the model state variables and <img src="https://latex.codecogs.com/svg.latex?\inline&space;\large&space;\boldsymbol{p}" title="p"/> is a vector of parameters.
Specifically, if one aims to minimize an objective function <img src="https://latex.codecogs.com/svg.latex?\inline&space;\large&space;\hat{f}(\boldsymbol{p})&space;=&space;f(\boldsymbol{s},&space;\boldsymbol{p})" title="\boldsymbol{f}" /> where <img src="https://latex.codecogs.com/svg.latex?\inline&space;\large&space;\boldsymbol{s}" title="s"/> is the steady-state solution of the system, i.e., such that  <img src="https://latex.codecogs.com/svg.latex?\inline&space;\large&space;F(\boldsymbol{s},&space;\boldsymbol{p})&space;=&space;0" title="\boldsymbol{sdef}" />, the F-1 method allows for an efficient computation of the Hessian matrix.
(As long as the Jacobian matrix <img src="https://latex.codecogs.com/svg.latex?\inline&space;\large&space;\mathbf{A} = \nabla_{\boldsymbol{x}}F(\boldsymbol{x},&space;\boldsymbol{p})&space;=&space;0" title="\boldsymbol{A}" /> can be created, factorized, and stored.)
