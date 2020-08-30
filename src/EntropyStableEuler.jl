"""
    Module EntropyStableEuler

Includes general math tools
"""

module EntropyStableEuler
using StaticArrays

const γ = 1.4

export logmean
include("./logmean.jl")

export primitive_to_conservative
export u_vfun, v_ufun, Sfun
include("./entropy_variables_2D.jl")

# export u_vfun, v_ufun, betafun, pfun, rhoe_ufun
# export dVdU_explicit, dUdV_explicit
# export wavespeed # c
# export Sfun, sfun # math/physical entropies
# export u_vfun1D, v_ufun1D, betafun1D # specialization to 1D
# export primitive_to_conservative
# include("./euler_variables.jl")

export euler_fluxes
export euler_flux_x, euler_flux_y # separate x,y fluxes for faster implicit assembly using ForwardDiff
include("./euler_fluxes.jl")

export vortex
include("./analytic_solutions.jl")

end
