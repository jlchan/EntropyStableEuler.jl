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
export u_vfun, v_ufun

# dispatch to n-dimensional constitutive routines, with optional entropy scaling
function primitive_to_conservative(rho,u,v,p)
   rho,rhoU,E = primitive_to_conservative_nd(rho,(u,v),p)
   return rho,rhoU...,E
end
function v_ufun(rho,rhou,rhov,E)
    v1,vU,vE = v_ufun_nd(rho,(rhou,rhov),E)
    return scale_entropy_vars_output(v1,vU...,vE)
end
function u_vfun(v1,vU1,vU2,vE)
    v1,vU1,vU2,vE = scale_entropy_vars_input(v1,vU1,vU2,vE)
    rho,rhoU,E = u_vfun_nd(v1,(vU1,vU2),vE)
    return rho,rhoU...,E
end
function conservative_to_primitive_beta(rho,rhou,rhov,E)
    rho,U,beta = conservative_to_primitive_beta_nd(rho,(rhou,rhov),E)
    return rho,U...,beta
end

export Sfun,pfun
Sfun(rho,rhou,rhov,E) = Sfun_nd(rho,(rhou,rhov),E)
pfun(rho,rhou,rhov,E) = pfun_nd(rho,(rhou,rhov),E)
include("./entropy_variables.jl")

export euler_fluxes_2D, euler_fluxes_2D_x, euler_fluxes_2D_y
include("./euler_fluxes_2D.jl")
end

# export u_vfun, v_ufun, betafun, pfun, rhoe_ufun
# export dVdU_explicit, dUdV_explicit
# export wavespeed # c
# export Sfun, sfun # math/physical entropies
# export u_vfun1D, v_ufun1D, betafun1D # specialization to 1D
# export primitive_to_conservative
# include("./euler_variables.jl")

# export euler_fluxes
# export euler_flux_x, euler_flux_y # separate x,y fluxes for faster implicit assembly using ForwardDiff
# include("./euler_fluxes.jl")

export vortex
include("./analytic_solutions.jl")

end
