"""
    Module EntropyStableEuler

Includes general math tools
"""

module EntropyStableEuler
using StaticArrays

const γ = 1.4

export logmean
include("./logmean.jl")

# submodules
export Fluxes1D,Fluxes2D,Fluxes3D

# for lax-friedrichs flux penalization - make available to all submodules
export wavespeed_1D

# export entropy_scaling,scale_entropy_input,scale_entropy_output
# changes definition of entropy variables by a constant scaling
const entropy_scaling = 1/γ # implies vE = -1 /(ι * γ) and γι = non-dim temperature if R* = cp
scale_entropy_output(V...) = (x->@. x*entropy_scaling).(V)
scale_entropy_input(V...) = (x->@. x/entropy_scaling).(V)

#####
##### one-dimensional fluxes
#####
module Fluxes1D
import ..γ
import ..entropy_scaling
import ..EntropyStableEuler: logmean
import ..EntropyStableEuler: scale_entropy_output, scale_entropy_input

export primitive_to_conservative,conservative_to_primitive_beta
export u_vfun, v_ufun

export wavespeed_1D

# dispatch to n-dimensional constitutive routines, with optional entropy scaling
function primitive_to_conservative(rho,u,p)
   rho,rhou,E = primitive_to_conservative_nd(rho,tuple(u),p)
   return rho,rhou[1],E
end
function v_ufun(rho,rhou,E)
    v1,vU,vE = v_ufun_nd(rho,tuple(rhou),E)
    return scale_entropy_output(v1,vU[1],vE)
end
function u_vfun(v1,vU,vE)
    v1,vU,vE = scale_entropy_input(v1,vU,vE)
    rho,rhoU,E = u_vfun_nd(v1,tuple(vU),vE)
    return rho,rhoU[1],E
end
function conservative_to_primitive_beta(rho,rhou,E)
    rho,U,beta = conservative_to_primitive_beta_nd(rho,tuple(rhou),E)
    return rho,U[1],beta
end

export Sfun,pfun
Sfun(rho,rhou,E) = Sfun_nd(rho,tuple(rhou),E)
pfun(rho,rhou,E) = pfun_nd(rho,tuple(rhou),E)
include("./entropy_variables.jl")

export euler_fluxes_1D
include("./euler_fluxes.jl")
end

import .Fluxes1D: wavespeed_1D

#####
##### two-dimensional fluxes
#####
module Fluxes2D
import ..γ
import ..entropy_scaling
import ..EntropyStableEuler: logmean, entropy_scaling
import ..EntropyStableEuler: scale_entropy_output, scale_entropy_input

export primitive_to_conservative,conservative_to_primitive_beta
export u_vfun, v_ufun

# dispatch to n-dimensional constitutive routines, with optional entropy scaling
function primitive_to_conservative(rho,u,v,p)
   rho,rhoU,E = primitive_to_conservative_nd(rho,(u,v),p)
   return rho,rhoU...,E
end
function v_ufun(rho,rhou,rhov,E)
    v1,vU,vE = v_ufun_nd(rho,(rhou,rhov),E)
    return scale_entropy_output(v1,vU...,vE)
end
function u_vfun(v1,vU1,vU2,vE)
    v1,vU1,vU2,vE = scale_entropy_input(v1,vU1,vU2,vE)
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
include("./euler_fluxes.jl")
end

#####
##### three-dimensional fluxes
#####
module Fluxes3D
import ..γ
import ..entropy_scaling
import ..EntropyStableEuler: logmean, entropy_scaling
import ..EntropyStableEuler: scale_entropy_output, scale_entropy_input

export primitive_to_conservative,conservative_to_primitive_beta
export u_vfun, v_ufun

# dispatch to n-dimensional constitutive routines, with optional entropy scaling
function primitive_to_conservative(rho,u,v,w,p)
   rho,rhoU,E = primitive_to_conservative_nd(rho,(u,v,w),p)
   return rho,rhoU...,E
end
function v_ufun(rho,rhou,rhov,rhow,E)
    v1,vU,vE = v_ufun_nd(rho,(rhou,rhov,rhow),E)
    return scale_entropy_output(v1,vU...,vE)
end
function u_vfun(v1,vU1,vU2,vU3,vE)
    v1,vU1,vU2,vU3,vE = scale_entropy_input(v1,vU1,vU2,vU3,vE)
    rho,rhoU,E = u_vfun_nd(v1,(vU1,vU2,vU3),vE)
    return rho,rhoU...,E
end
function conservative_to_primitive_beta(rho,rhou,rhov,rhow,E)
    rho,U,beta = conservative_to_primitive_beta_nd(rho,(rhou,rhov,rhow),E)
    return rho,U...,beta
end

export Sfun,pfun
Sfun(rho,rhou,rhov,rhow,E) = Sfun_nd(rho,(rhou,rhov,rhow),E)
pfun(rho,rhou,rhov,rhow,E) = pfun_nd(rho,(rhou,rhov,rhow),E)
include("./entropy_variables.jl")

export euler_fluxes_3D
include("./euler_fluxes.jl")
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

# export vortex
# include("./analytic_solutions.jl")

end
