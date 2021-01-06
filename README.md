# EntropyStableEuler

[![Build Status](https://travis-ci.com/jlchan/EntropyStableEuler.jl.svg?branch=master)](https://travis-ci.com/jlchan/EntropyStableEuler.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/jlchan/EntropyStableEuler.jl?svg=true)](https://ci.appveyor.com/project/jlchan/EntropyStableEuler-jl)
[![Build status](https://github.com/jlchan/EntropyStableEuler.jl/workflows/CI/badge.svg)](https://github.com/jlchan/EntropyStableEuler.jl/actions)
[![Codecov](https://codecov.io/gh/jlchan/EntropyStableEuler.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jlchan/EntropyStableEuler.jl)

Entropy stable finite volume fluxes and formulas for compressible Euler and Navier-Stokes. Code based off of formulas in [Chandrashekar 2012](https://doi.org/10.4208/cicp.170712.010313a) and [Winters et al. 2019](https://link.springer.com/article/10.1007/s10543-019-00789-w). Formulas for entropy variables are from [Hughes, Mallet, Franca 1986](https://doi.org/10.1016/0045-7825(86)90127-1) and from [Parsani, Carpenter, Nielsen 2015](https://doi.org/10.1016/j.jcp.2015.03.026)

# Usage

The package exports
- The `Euler{d}` type, which dispatches to 1D, 2D, and 3D Euler
- `fS(Euler{d}(),UL,UR)`, which evaluates entropy conservative fluxes in d-dimensions using conservative variables
- `fS_prim(QL,QR)`, which evaluates entropy conservative fluxes using primitive variables. `fS_prim(QL,QR,QlogL,QlogR)` uses precomputed logs too.
- `u_vfun(V), v_ufun(U)` to convert between conservative variables `U` and entropy variables `V`
- `cons_to_prim_beta` to convert between conservative and "primitive" variables (involving inverse temperature β) used to evaluate fluxes.

# Example 
```
using EntropyStableEuler

# construct solution at two states
UL = map(x->x.*ones(4),(1,.1,.2,2))
UR = map(x->x.+.1*randn(4),(rhoL,rhouL,rhovL,EL)) # small perturbation

# evaluate fluxes
Fx,Fy = fS(Euler{2}(),UL,UR)

# pass in primitive vars/precomputed logs for efficiency
QL,QR = cons_to_prim_beta.((UL,UR))
Fx,Fy = fS_prim(Euler{2}(),QL,QR)

QlogL = map(x->log.(x),(first(QL),last(QL)))
QlogR = map(x->log.(x),(first(QR),last(QR)))
Fx,Fy = fS_prim(Euler{2}(),QL,QR,QlogL,QlogR)
```

# To-do
- add Lax-Friedrichs penalty and matrix dissipation from [Winters et al. 2017](https://doi.org/10.1016/j.jcp.2016.12.006)
- Jacobians for transforms between conservative and entropy variables
- viscous entropy variable matrices for compressible Navier-Stokes
