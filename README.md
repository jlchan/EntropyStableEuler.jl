# EntropyStableEuler

[![Build Status](https://travis-ci.com/jlchan/EntropyStableEuler.jl.svg?branch=master)](https://travis-ci.com/jlchan/EntropyStableEuler.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/jlchan/EntropyStableEuler.jl?svg=true)](https://ci.appveyor.com/project/jlchan/EntropyStableEuler-jl)
[![Codecov](https://codecov.io/gh/jlchan/EntropyStableEuler.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jlchan/EntropyStableEuler.jl)

Entropy stable finite volume fluxes and formulas for compressible Euler and Navier-Stokes. Code based off of formulas in [Chandrashekar 2012](https://doi.org/10.4208/cicp.170712.010313a) and [Winters et al. 2019](https://link.springer.com/article/10.1007/s10543-019-00789-w). Formulas for entropy variables are from [Hughes, Mallet, Franca 1986](https://doi.org/10.1016/0045-7825(86)90127-1) and from [Parsani, Carpenter, Nielsen 2015](https://doi.org/10.1016/j.jcp.2015.03.026)

# To-do
- add Lax-Friedrichs penalty and matrix dissipation from [Winters et al. 2017](https://doi.org/10.1016/j.jcp.2016.12.006)
- add Jacobians for transforms between conservative and entropy variables
