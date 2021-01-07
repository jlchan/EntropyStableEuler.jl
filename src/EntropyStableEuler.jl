module EntropyStableEuler

using StaticArrays
using UnPack

export Euler
Base.@kwdef struct Euler{d}
    γ::Float64 = 1.4 # default value of γ
end

export logmean
include("logmean.jl")

export betafun, pfun, prim_to_cons, cons_to_prim_beta
export Sfun, u_vfun, v_ufun
include("variables.jl")

export fS, fS_prim
include("fluxes.jl")

end
