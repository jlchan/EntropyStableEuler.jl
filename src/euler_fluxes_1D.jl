"function euler_fluxes_1D(rhoL,uL,betaL,rhoR,uR,betaR)"
function euler_fluxes_1D(rhoL,uL,betaL,rhoR,uR,betaR)
    rhologL,betalogL,rhologR,betalogR = log.((rhoL,betaL,rhoR,betaR))
    return euler_fluxes_1D(rhoL,uL,betaL,rhoR,uR,betaR,
                           rhologL,betalogL,rhologR,betalogR)
end

"function euler_fluxes_1D(rhoL,uL,betaL,rhoR,uR,betaR,
                         rhologL,betalogL,rhologR,betalogR)"
function euler_fluxes_1D(rhoL,uL,betaL,rhoR,uR,betaR,
                         rhologL,betalogL,rhologR,betalogR)

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
