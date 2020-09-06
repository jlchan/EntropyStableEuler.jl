#####
##### helper functions
#####

"function wavespeed_1D(rho,rhou,E)
    one-dimensional wavespeed (for DG penalization terms)"
function wavespeed_1D(rho,rhou,E)
    cvel = @. sqrt(γ*pfun(rho,rhou,E)/rho)
    return @. abs(rhou/rho) + cvel
end

#####
##### one-dimensional fluxes
#####

"function euler_fluxes_1D(rhoL,uL,betaL,rhoR,uR,betaR)"
function euler_fluxes_1D(rhoL,uL,betaL,rhoR,uR,betaR)
    rhologL,betalogL,rhologR,betalogR = log.((rhoL,betaL,rhoR,betaR))
    return euler_fluxes_1D(rhoL,uL,betaL,rhologL,betalogL,
                           rhoR,uR,betaR,rhologR,betalogR)
end

"function euler_fluxes_1D(rhoL,uL,betaL,rhologL,betalogL,
                        rhoR,uR,betaR,rhologR,betalogR)"
function euler_fluxes_1D(rhoL,uL,betaL,rhologL,betalogL,
                        rhoR,uR,betaR,rhologR,betalogR)

    rholog = logmean.(rhoL,rhoR,rhologL,rhologR)
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

    return (FxS1,FxS2,FxS3)
end

#####
##### two-dimensional fluxes
#####
"function euler_fluxes_2D(rhoL,uL,vL,betaL,rhoR,uR,vR,betaR)
assumes primitive variables ordering: UL = (rhoL,uL,...,betaL),
                                      UR = (rhoR,uR,...,betaR)"
function euler_fluxes_2D(rhoL,uL,vL,betaL,rhoR,uR,vR,betaR)
    rhologL,betalogL,rhologR,betalogR = log.((rhoL,betaL,rhoR,betaR))
    return euler_fluxes_2D(rhoL,uL,vL,betaL,rhologL,betalogL,
                           rhoR,uR,vR,betaR,rhologR,betalogR)
end

"
function euler_fluxes_2D(rhoL,uL,vL,betaL,rhoR,uR,vR,betaR,
                         rhologL,betalogL,rhologR,betalogR)
assumes primitive variables ordering: UL = (rhoL,uL,...,betaL),
                                      UR = (rhoR,uR,...,betaR)
"
function euler_fluxes_2D(rhoL,uL,vL,betaL,rhologL,betalogL,
                         rhoR,uR,vR,betaR,rhologR,betalogR)

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
    return (FxS1,FxS2,FxS3,FxS4),(FyS1,FyS2,FyS3,FyS4)
end

function euler_fluxes_2D_x(rhoL,uL,vL,betaL,rhoR,uR,vR,betaR)
    rhologL,betalogL,rhologR,betalogR = log.((rhoL,betaL,rhoR,betaR))
    return euler_fluxes_2D_x(rhoL,uL,vL,betaL,rhologL,betalogL,
                             rhoR,uR,vR,betaR,rhologR,betalogR)
end
function euler_fluxes_2D_y(rhoL,uL,vL,betaL,rhoR,uR,vR,betaR)
    rhologL,betalogL,rhologR,betalogR = log.((rhoL,betaL,rhoR,betaR))
    return euler_fluxes_2D_y(rhoL,uL,vL,betaL,rhologL,betalogL,
                             rhoR,uR,vR,betaR,rhologR,betalogR)
end

function euler_fluxes_2D_x(rhoL,uL,vL,betaL,rhologL,betalogL,
                           rhoR,uR,vR,betaR,rhologR,betalogR)

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

    return (FxS1,FxS2,FxS3,FxS4)
end

function euler_fluxes_2D_y(rhoL,uL,vL,betaL,rhologL,betalogL,
                           rhoR,uR,vR,betaR,rhologR,betalogR)

    rholog = logmean.(rhoL,rhoR,rhologL,rhologR)
    betalog = logmean.(betaL,betaR,betalogL,betalogR)

    # arithmetic avgs
    rhoavg = (@. .5*(rhoL+rhoR))
    uavg   = (@. .5*(uL+uR))
    vavg   = (@. .5*(vL+vR))

    unorm = (@. uL*uR + vL*vR)
    pa    = (@. rhoavg/(betaL+betaR))
    f4aux = (@. rholog/(2*(γ-1)*betalog) + pa + .5*rholog*unorm)

    FyS1 = (@. rholog*vavg)
    FyS2 = (@. FyS1*uavg)
    FyS3 = (@. FyS1*vavg + pa)
    FyS4 = (@. f4aux*vavg)
    return (FyS1,FyS2,FyS3,FyS4)
end

#####
##### three-dimensional fluxes
#####
"function euler_fluxes_3D(rhoL,uL,vL,wL,betaL,rhoR,uR,vR,wR,betaR)"
function euler_fluxes_3D(rhoL,uL,vL,wL,betaL,rhoR,uR,vR,wR,betaR)
    rhologL,betalogL,rhologR,betalogR = log.((rhoL,betaL,rhoR,betaR))
    return euler_fluxes_3D(rhoL,uL,vL,wL,betaL,rhologL,betalogL,
                           rhoR,uR,vR,wR,betaR,rhologR,betalogR)
end

"function euler_fluxes_3D(rhoL,uL,vL,wL,betaL,rhologL,betalogL,
                         rhoR,uR,vR,wR,betaR,rhologR,betalogR)"
function euler_fluxes_3D(rhoL,uL,vL,wL,betaL,rhologL,betalogL,
                         rhoR,uR,vR,wR,betaR,rhologR,betalogR)

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

    Fx = (FxS1,FxS2,FxS3,FxS4,FxS5)
    Fy = (FyS1,FyS2,FyS3,FyS4,FyS5)
    Fz = (FzS1,FzS2,FzS3,FzS4,FzS5)
    return Fx,Fy,Fz
end
