"""
    fS(eqn::Euler{d},UL,UR) where {d}

Entropy conservative fluxes for the compressible Euler equations, where UL,UR are
tuples of left and right states in terms of the conservative variables Q = (ρ,ρuvw...,E)
"""
function fS(eqn::Euler{d},UL,UR) where {d}
    QL = cons_to_prim_beta(eqn,UL)
    QR = cons_to_prim_beta(eqn,UR)
    return fS_prim(eqn,QL,QR)
end

"""
    fS_prim(eqn::Euler{d},UL,UR)
    #fS_prim(eqn::Euler{d},UL,UR,rhologL,betalogL,rhologR,betalogR)

Entropy conservative fluxes for the compressible Euler equations, where
UL,UR are tuples of left and right solutions in terms of β-primitive variables
Q = (rho,uvw...,β) (where uvw = velocity components and β is inverse temperature).
"""
function fS_prim(eqn::Euler{d},UL,UR) where {d}
    UlogL = map(x->log.(x),(first(UL),last(UL)))
    UlogR = map(x->log.(x),(first(UR),last(UR)))
    return fS_prim_log(eqn,(UL...,UlogL...),(UR...,UlogR...))
end

"""
    fS_prim_log(eqn::Euler{d},UL,UR)

Entropy conservative fluxes for the compressible Euler equations, where
UL,UR are tuples of left and right solutions (in β-primitive variables) and logs
of rho, β, e.g., U = (rho,uvw...,β,log.(rho),log.(β)).
"""
function fS_prim_log(eqn::Euler{1},UL,UR)

    @unpack γ = eqn

    rhoL,uL,betaL,rhologL,betalogL = UL
    rhoR,uR,betaR,rhologR,betalogR = UR
    # rhologL,betalogL = UlogL
    # rhologR,betalogR = UlogR

    rholog  = logmean.(rhoL,rhoR,rhologL,rhologR)
    betalog = logmean.(betaL,betaR,betalogL,betalogR)

    # arithmetic avgs
    rhoavg = (@. .5*(rhoL+rhoR))
    uavg   = (@. .5*(uL+uR))

    unorm = (@. uL*uR)
    pa    = (@. rhoavg/(betaL+betaR))
    f4aux = (@. rholog/(2*(γ-1)*betalog) + pa + .5*rholog*unorm)

    FxS1 = (@. rholog*uavg)
    FxS2 = (@. FxS1*uavg + pa)
    FxS3 = (@. f4aux*uavg)

    return SVector(FxS1,FxS2,FxS3)
end

function fS_prim_log(eqn::Euler{2},UL,UR)

    @unpack γ = eqn

    rhoL,uL,vL,betaL,rhologL,betalogL = UL
    rhoR,uR,vR,betaR,rhologR,betalogR = UR
    # rhologL,betalogL = UlogL
    # rhologR,betalogR = UlogR

    rholog = logmean.(rhoL,rhoR,rhologL,rhologR)
    betalog = logmean.(betaL,betaR,betalogL,betalogR)

    # arithmetic avgs
    rhoavg = (@. .5*(rhoL+rhoR))
    uavg   = (@. .5*(uL+uR))
    vavg   = (@. .5*(vL+vR))

    unorm = (@. uL*uR + vL*vR)
    pa    = (@. rhoavg/(betaL+betaR))
    f4aux = (@. rholog/(2*(γ-1)*betalog) + pa + .5*rholog*unorm)

    FxS1 = (@. rholog*uavg)
    FxS2 = (@. FxS1*uavg + pa)
    FxS3 = (@. FxS1*vavg)
    FxS4 = (@. f4aux*uavg)

    FyS1 = (@. rholog*vavg)
    FyS2 = (@. FxS3)
    FyS3 = (@. FyS1*vavg + pa)
    FyS4 = (@. f4aux*vavg)

    return SVector(FxS1,FxS2,FxS3,FxS4),SVector(FyS1,FyS2,FyS3,FyS4)
end

function fS_prim_log(eqn::Euler{3},UL,UR)

    @unpack γ = eqn

    rhoL,uL,vL,wL,betaL,rhologL,betalogL = UL
    rhoR,uR,vR,wR,betaR,rhologR,betalogR = UR
    # rhologL,betalogL = UlogL
    # rhologR,betalogR = UlogR

    rholog = logmean.(rhoL,rhoR,rhologL,rhologR)
    betalog = logmean.(betaL,betaR,betalogL,betalogR)

    # arithmetic avgs
    rhoavg = (@. .5*(rhoL+rhoR))
    uavg   = (@. .5*(uL+uR))
    vavg   = (@. .5*(vL+vR))
    wavg   = (@. .5*(wL+wR))

    unorm = (@. uL*uR + vL*vR + wL*wR)
    pa    = (@. rhoavg/(betaL+betaR))
    E_plus_p  = (@. rholog/(2*(γ-1)*betalog) + pa + .5*rholog*unorm)

    FxS1 = (@. rholog*uavg)
    FxS2 = (@. FxS1*uavg + pa)
    FxS3 = (@. FxS1*vavg) # rho * u * v
    FxS4 = (@. FxS1*wavg) # rho * u * w
    FxS5 = (@. E_plus_p*uavg)

    FyS1 = (@. rholog*vavg)
    FyS2 = (@. FxS3) # rho * u * v
    FyS3 = (@. FyS1*vavg + pa)
    FyS4 = (@. FyS1*wavg) # rho * v * w
    FyS5 = (@. E_plus_p*vavg)

    FzS1 = (@. rholog*wavg)
    FzS2 = (@. FxS4) # rho * u * w
    FzS3 = (@. FyS4) # rho * v * w
    FzS4 = (@. FzS1*wavg + pa) # rho * w^2 + p
    FzS5 = (@. E_plus_p*wavg)

    Fx = SVector(FxS1,FxS2,FxS3,FxS4,FxS5)
    Fy = SVector(FyS1,FyS2,FyS3,FyS4,FyS5)
    Fz = SVector(FzS1,FzS2,FzS3,FzS4,FzS5)
    return Fx,Fy,Fz
end
