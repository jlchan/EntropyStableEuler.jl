"""
    Module EntropyStableEuler

Includes general math tools
"""

module EntropyStableEuler

using StaticArrays
using UnPack

abstract type AbstractEqn{d} end # d = dimension
Base.@kwdef struct Euler{d} <: AbstractEqn{d}
    γ::Float64 = 1.4 # default value of γ
end

export AbstractEqn,Euler

export logmean
include("logmean.jl")

export betafun, pfun, prim_to_cons, cons_to_prim_beta
export Sfun, u_vfun, v_ufun
include("variables.jl")

export fS, fS_prim
include("fluxes.jl")

# #####
# ##### one-dimensional fluxes
# #####
#
# module Fluxes1D
# using StaticArrays
# import ..γ
# import ..EntropyStableEuler: logmean
#
# export wavespeed_1D
#
# export primitive_to_conservative, conservative_to_primitive_beta
# export u_vfun, v_ufun
# export euler_flux_prim
# export Sfun,pfun,betafun
# include("./euler_fluxes_1D.jl")
# include("./entropy_variables.jl")
# end
#
# #####
# ##### two-dimensional fluxes
# #####
# module Fluxes2D
# using StaticArrays
# import ..γ
# import ..EntropyStableEuler: logmean
#
# export primitive_to_conservative,conservative_to_primitive_beta
# export u_vfun, v_ufun
# export Sfun,pfun,betafun
# export euler_flux_prim
#
# include("./entropy_variables.jl")
# include("./euler_fluxes_2D.jl")
# end
#
# #####
# ##### three-dimensional fluxes
# #####
# module Fluxes3D
# using StaticArrays
# import ..γ
# import ..EntropyStableEuler: logmean
#
# export primitive_to_conservative,conservative_to_primitive_beta
# export u_vfun, v_ufun
# export Sfun,pfun,betafun
# export euler_flux_prim
#
# include("./entropy_variables.jl")
# include("./euler_fluxes_3D.jl")
# end

end
