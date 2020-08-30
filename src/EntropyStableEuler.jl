"""
    Module EntropyStableEuler

Includes general math tools
"""

module EntropyStableEuler
using StaticArrays

const γ = 1.4

export logmean
include("./logmean.jl")

#####
##### two-dimensional fluxes
#####
module Fluxes2D
import ..γ
using ..EntropyStableEuler # for logmean
export primitive_to_conservative,conservative_to_primitive_beta
export u_vfun, v_ufun, Sfun
export euler_fluxes_2D, euler_fluxes_2D_x, euler_fluxes_2D_y
include("./entropy_variables_2D.jl")
include("./euler_fluxes_2D.jl")
end

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
