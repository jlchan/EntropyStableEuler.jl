# EntropyStableEuler

[![Build Status](https://travis-ci.com/jlchan/EntropyStableEuler.jl.svg?branch=master)](https://travis-ci.com/jlchan/EntropyStableEuler.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/jlchan/EntropyStableEuler.jl?svg=true)](https://ci.appveyor.com/project/jlchan/EntropyStableEuler-jl)
[![Build status](https://github.com/jlchan/EntropyStableEuler.jl/workflows/CI/badge.svg)](https://github.com/jlchan/EntropyStableEuler.jl/actions)
[![Codecov](https://codecov.io/gh/jlchan/EntropyStableEuler.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jlchan/EntropyStableEuler.jl)

Provides entropy stable finite volume fluxes for the compressible Euler equations, as well as formulas for transformations between conservative and entropy variables. 

Code based off of formulas in [Chandrashekar 2012](https://doi.org/10.4208/cicp.170712.010313a) and [Winters et al. 2019](https://link.springer.com/article/10.1007/s10543-019-00789-w). Formulas for entropy variables are from [Hughes, Mallet, Franca 1986](https://doi.org/10.1016/0045-7825(86)90127-1) and from [Parsani, Carpenter, Nielsen 2015](https://doi.org/10.1016/j.jcp.2015.03.026)

# Usage

The package exports
- The `Euler{d}` type, which is used to dispatch on 1D, 2D, and 3D formulas.
- `fS(Euler{d}(),UL,UR)`, which evaluates d-dimensional EC fluxes using conservative variables.
- `fS_prim(Euler{d}(),QL,QR)`, which evaluates entropy conservative fluxes using primitive variables. `fS_prim(Euler{d}(),QL,QR,QlogL,QlogR)` uses precomputed logs too.
- `u_vfun(Euler{d}(),V), v_ufun(Euler{d}(),U)` to convert between conservative variables `U` and entropy variables `V`
- `cons_to_prim_beta(Euler{d}(),U)` to convert between conservative and "primitive" variables (involving inverse temperature β) used to evaluate fluxes.

# Example

```julia
using EntropyStableEuler

# construct solution at two states
UL = map(x->x.*ones(4),(1,.1,.2,2))
UR = map(x->x.*ones(4),(1.1,.2,.3,2.5))

# evaluate fluxes
Fx,Fy = fS(Euler{2}(),UL,UR)

# pass in primitive vars/precomputed logs for efficiency
QL,QR = cons_to_prim_beta.(Euler{2}(),(UL,UR))
Fx,Fy = fS_prim(Euler{2}(),QL,QR)

QlogL = map(x->log.(x),(first(QL),last(QL)))
QlogR = map(x->log.(x),(first(QR),last(QR)))
Fx,Fy = fS_prim(Euler{2}(),QL,QR,QlogL,QlogR)
```

# To-do
- add Lax-Friedrichs penalty and matrix dissipation from [Winters et al. 2017](https://doi.org/10.1016/j.jcp.2016.12.006)
- Jacobians for transforms between conservative and entropy variables
- viscous entropy variable matrices for compressible Navier-Stokes
