unorm(U) = sum(map((x->x.^2),U))

function unpackfields(eqn::Euler{d},U) where {d}
    return first(U),ntuple(i->U[i+1],d),last(U)
end

"""
    function prim_to_cons(eqn::Euler{d},rho,U,p) where {d}

convert primitive variables (ρ,U...,p) to conservative vars (ρ,ρU,E).
n-dimensional version where U = tuple(u1,...,u_d)
"""
function prim_to_cons(eqn::Euler{d},Q) where {d}
    rho,U,p = unpackfields(eqn,Q)
    rhoU = map(x->rho.*x,U)
    Unorm = unorm(U)
    E = @. p./(eqn.γ-1) + .5*rho*Unorm
    return SVector{d+2}(rho,rhoU...,E)
end

"""
    function cons_to_prim_beta(rho,rhoU,E)

converts conservative variables to `primitive' variables which make evaluating EC
fluxes cheaper.
"""
function cons_to_prim_beta(eqn::Euler{d},U) where {d}
    rho,rhoU,E = unpackfields(eqn,U)
    return SVector{d+2}(rho, map(x->x./rho,rhoU)..., betafun(eqn,U))
end

function pfun(eqn::Euler{d},U) where {d}
    rho,rhoU,E = unpackfields(eqn,U)
    rhoUnorm = unorm(rhoU)
    return @. (eqn.γ-1)*(E-.5*rhoUnorm./rho)
end

function betafun(eqn::Euler{d},U) where {d}
    rho = first(U)
    p = pfun(eqn,U)
    return @. rho/(2*p)
end

function sfun(eqn::Euler{d},U) where {d}
    rho = first(U)
    p = pfun(eqn,U)
    return @. log(p/(rho^eqn.γ))
end

function Sfun(eqn::Euler{d},U) where {d}
    rho = first(U)
    return -rho.*sfun(eqn,U)
end

function v_ufun(eqn::Euler{d}, U) where {d}
    @unpack γ = eqn

    s = sfun(eqn,U)
    p = pfun(eqn,U)

    rho,rhoU,E = unpackfields(eqn,U)
    v1 = @. (γ+1-s) - (γ-1)*E/p
    vU = (x->@. x*(γ-1)/p).(rhoU)
    vE = (x->@. x*(γ-1)/p)(-rho)

    return SVector(v1,vU...,vE)
end

function s_vfun(eqn::Euler{d},V) where {d}
    @unpack γ = eqn
    v1,vU,vE = unpackfields(eqn,V)
    vUnorm = unorm(vU)
    return @. γ - v1 + vUnorm/(2*vE)
end

function rhoe_vfun(eqn::Euler{d},V) where {d}
    @unpack γ = eqn
    s = s_vfun(eqn,V)
    v1,vU,vE = unpackfields(eqn,V)
    return @. ((γ-1)/((-vE)^γ))^(1/(γ-1)) * exp(-s/(γ-1))
end

function u_vfun(eqn::Euler{d},V) where {d}
    v1,vU,vE  = unpackfields(eqn,V)
    rhoeV     = rhoe_vfun(eqn,V)
    vUnorm    = unorm(vU)
    rho       = (@. rhoeV*(-vE))
    rhoU      = (x->rhoeV.*x).(vU)
    E         = (@. rhoeV*(1-vUnorm/(2*vE)))
    return SVector{d+2}(rho,rhoU...,E)
end

"""
    wavespeed(eqn::Euler{1},U)

one-dimensional wavespeed (for DG penalization terms)
"""
function wavespeed(eqn::Euler{1},U)
    rho,rhou,_ = unpackfields(eqn,U)
    p = pfun(eqn,U)
    cvel = @. sqrt(γ*p/rho)
    return @. abs(rhou/rho) + cvel
end

# function dUdV_explicit(v1,vU1,vU2,vE)
#     rho,rhou,rhov,E = u_vfun(v1,vU1,vU2,vE)
#     u,v = (x->x./rho).((rhou,rhov))
#     p = pfun(rho,rhou,rhov,E)
#     a2 = γ*p/rho
#     H = a2/(γ-1) + (u^2+v^2)/2
#
#     dUdV = @SMatrix [rho  rhou        rhov        E;
#                      rhou rhou*u + p  rhou*v      rhou*H;
#                      rhov rhov*u      rhov*v + p  rhov*H;
#                      E    rhou*H      rhov*H      rho*H^2-a2*p/(γ-1)]
#
#     return dUdV*(1/(γ-1))
# end
#
# function dVdU_explicit(rho,rhou,rhov,E)
#     rhoe = rhoe_ufun(rho,rhou,rhov,E)
#     V = v_ufun(rho,rhou,rhov,E)
#     k = .5*(V[2]^2+V[3]^2)/V[4]
#
#     dVdU = @SMatrix [γ+k^2      k*V[2]          k*V[3]         V[4]*(k+1);
#                     k*V[2]      V[2]^2-V[4]     V[2]*V[3]      V[2]*V[4];
#                     k*V[3]      V[2]*V[3]       V[3]^2-V[4]    V[3]*V[4]
#                     V[4]*(k+1)  V[2]*V[4]       V[3]*V[4]      V[4]^2]
#     return -dVdU/(rhoe*V[4])
# end
