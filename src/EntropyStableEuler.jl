module EntropyStableEuler

using StaticArrays
using UnPack

export Euler
Base.@kwdef struct Euler{d}
    γ::Float64 = 1.4 # default value of γ
end

export nfields
nfields(equation::Euler{d}) where {d} = d+2

export logmean
include("logmean.jl")

export betafun, pfun, prim_to_cons, cons_to_prim_beta, cons_to_prim_beta_log
export Sfun, cons_to_entropy, entropy_to_cons
include("variables.jl")

export fS, fS_prim, fS_prim_log
export normal_wavespeed, LxF_dissipation
include("fluxes.jl")

end
