"""
    fS(eqn::Euler{d},UL,UR) where {d}

Entropy conservative fluxes for the compressible Euler equations, where UL,UR are
tuples of left and right states in terms of the conservative variables Q = (ρ,ρuvw...,E)
"""
@inline function fS(eqn::Euler{d},UL,UR) where {d}
    QL = cons_to_prim_beta(eqn,UL)
    QR = cons_to_prim_beta(eqn,UR)
    return fS_prim(eqn,QL,QR)
end

"""
    fS_prim(eqn::Euler{d},UL,UR)

Entropy conservative fluxes for the compressible Euler equations, where
UL,UR are tuples of left and right solutions in terms of β-primitive variables
Q = (rho,uvw...,β) (where uvw = velocity components and β is inverse temperature).
"""
@inline function fS_prim(eqn::Euler{d},UL,UR) where {d}
    logL = map(x->log(x),(first(UL),last(UL)))
    logR = map(x->log(x),(first(UR),last(UR)))

    return fS_prim_log(eqn,(UL...,logL...),(UR...,logR...))
end

# manually unroll for d = 1,2,3 to avoid type instability
@inline function fS_prim(eqn::Euler{1},UL,UR)
    UlogL = (UL[1],UL[2],UL[3],log(UL[1]),log(UL[3]))
    UlogR = (UR[1],UR[2],UR[3],log(UR[1]),log(UR[3]))
    return fS_prim_log(eqn,UlogL,UlogR)
end

@inline function fS_prim(eqn::Euler{2},UL,UR)
    UlogL = (UL[1],UL[2],UL[3],UL[4],log(UL[1]),log(UL[4]))
    UlogR = (UR[1],UR[2],UR[3],UR[4],log(UR[1]),log(UR[4]))
    return fS_prim_log(eqn,UlogL,UlogR)
end

@inline function fS_prim(eqn::Euler{3},UL,UR)
    UlogL = (UL[1],UL[2],UL[3],UL[4],UL[5],log(UL[1]),log(UL[5]))
    UlogR = (UR[1],UR[2],UR[3],UR[4],UR[5],log(UR[1]),log(UR[5]))
    return fS_prim_log(eqn,UlogL,UlogR)
end

"""
    fS_prim_log(eqn::Euler{d},UL,UR)

Entropy conservative fluxes for the compressible Euler equations, where
UL,UR are tuples of left and right solutions (in β-primitive variables) and logs
of rho, β, e.g., U = (rho,uvw...,β,log(rho),log(β)).
"""
@inline function fS_prim_log(eqn::Euler{1},UL,UR)

    @unpack γ = eqn

    rhoL,uL,betaL,rhologL,betalogL = UL
    rhoR,uR,betaR,rhologR,betalogR = UR

    rholog  = logmean(rhoL,rhoR,rhologL,rhologR)
    betalog = logmean(betaL,betaR,betalogL,betalogR)

    # arithmetic avgs
    rhoavg = ( .5*(rhoL+rhoR))
    uavg   = ( .5*(uL+uR))

    unorm = ( uL*uR)
    pa    = ( rhoavg/(betaL+betaR))
    f4aux = ( rholog/(2*(γ-1)*betalog) + pa + .5*rholog*unorm)

    FxS1 = ( rholog*uavg)
    FxS2 = ( FxS1*uavg + pa)
    FxS3 = ( f4aux*uavg)

    return SVector(FxS1,FxS2,FxS3)
end

@inline function fS_prim_log(eqn::Euler{2},UL,UR)

    @unpack γ = eqn

    rhoL,uL,vL,betaL,rhologL,betalogL = UL
    rhoR,uR,vR,betaR,rhologR,betalogR = UR

    rholog = logmean(rhoL,rhoR,rhologL,rhologR)
    betalog = logmean(betaL,betaR,betalogL,betalogR)

    # arithmetic avgs
    rhoavg = ( .5*(rhoL+rhoR))
    uavg   = ( .5*(uL+uR))
    vavg   = ( .5*(vL+vR))

    unorm = ( uL*uR + vL*vR)
    pa    = ( rhoavg/(betaL+betaR))
    f4aux = ( rholog/(2*(γ-1)*betalog) + pa + .5*rholog*unorm)

    FxS1 = ( rholog*uavg)
    FxS2 = ( FxS1*uavg + pa)
    FxS3 = ( FxS1*vavg)
    FxS4 = ( f4aux*uavg)

    FyS1 = ( rholog*vavg)
    FyS2 = ( FxS3)
    FyS3 = ( FyS1*vavg + pa)
    FyS4 = ( f4aux*vavg)

    return SVector(FxS1,FxS2,FxS3,FxS4),SVector(FyS1,FyS2,FyS3,FyS4)
end

@inline function fS_prim_log(eqn::Euler{3},UL,UR)

    @unpack γ = eqn

    rhoL,uL,vL,wL,betaL,rhologL,betalogL = UL
    rhoR,uR,vR,wR,betaR,rhologR,betalogR = UR

    rholog = logmean(rhoL,rhoR,rhologL,rhologR)
    betalog = logmean(betaL,betaR,betalogL,betalogR)

    # arithmetic avgs
    rhoavg = ( .5*(rhoL+rhoR))
    uavg   = ( .5*(uL+uR))
    vavg   = ( .5*(vL+vR))
    wavg   = ( .5*(wL+wR))

    unorm = ( uL*uR + vL*vR + wL*wR)
    pa    = ( rhoavg/(betaL+betaR))
    E_plus_p  = ( rholog/(2*(γ-1)*betalog) + pa + .5*rholog*unorm)

    FxS1 = ( rholog*uavg)
    FxS2 = ( FxS1*uavg + pa)
    FxS3 = ( FxS1*vavg) # rho * u * v
    FxS4 = ( FxS1*wavg) # rho * u * w
    FxS5 = ( E_plus_p*uavg)

    FyS1 = ( rholog*vavg)
    FyS2 = ( FxS3) # rho * u * v
    FyS3 = ( FyS1*vavg + pa)
    FyS4 = ( FyS1*wavg) # rho * v * w
    FyS5 = ( E_plus_p*vavg)

    FzS1 = ( rholog*wavg)
    FzS2 = ( FxS4) # rho * u * w
    FzS3 = ( FyS4) # rho * v * w
    FzS4 = ( FzS1*wavg + pa) # rho * w^2 + p
    FzS5 = ( E_plus_p*wavg)

    Fx = SVector(FxS1,FxS2,FxS3,FxS4,FxS5)
    Fy = SVector(FyS1,FyS2,FyS3,FyS4,FyS5)
    Fz = SVector(FzS1,FzS2,FzS3,FzS4,FzS5)
    return Fx,Fy,Fz
end

### didn't add coordinates b/c timing doesn't seem to make a diff.
# julia> Qlog = (1.0, 0.1, 0.2, 2.5250000000000004, 0.0, 0.9262410627273233)
# julia> @btime fS_prim_log($Euler{2}(),$Qlog,$Qlog);
#   2.963 μs (0 allocations: 0 bytes)
#
# julia> @btime fS_prim_log($Euler{2}(),$Qlog,$Qlog,$1);
#   3.044 μs (0 allocations: 0 bytes)
