unorm(U) = sum((x->x.^2).(U))

"function primitive_to_conservative(rho,u,v,p)

    convert primitive variables (ρ,u,v,p) to conservative vars (ρ,ρu,ρv,E)."
function primitive_to_conservative_nd(rho,U,p)
    rhoU = (x->rho.*x).(U)
    Unorm = unorm(U)
    E = @. p/(γ-1) + .5*rho*Unorm
    return (rho,rhoU,E)
end

#####
##### functions of conservative variables
#####
"function pfun(rho, rhou, rhov, E)
    pressure as a function of conservative variables"
function pfun_nd(rho, rhoU, E)
    rhoUnorm2 = unorm(rhoU)./rho
    return @. (γ-1)*(E - .5*rhoUnorm2)
end

"function betafun(rho,rhou,rhov,E)
    inverse temperature (used in entropy conservative fluxes)"
function betafun_nd(rho,rhoU,E)
    p = pfun_nd(rho,rhoU,E)
    return (@. rho/(2*p))
end

"function rhoe_ufun(rho, rhou, rhov, E)
    specific energy as a function of conservative variables"
function rhoe_ufun_nd(rho, rhoU, E)
    return pfun_nd(rho, rhoU, E) / (γ-1)
end

"function sfun(rho, rhou, rhov, E)
    Specific entropy as a function of conservative variables"
function sfun_nd(rho, rhoU, E)
    p = pfun_nd(rho, rhoU, E)
    return @. log(p/(rho^γ))
end

"function Sfun(rho,rhou,rhov,E)
    Mathematical entropy as a function of conservative variables"
function Sfun_nd(rho, rhoU, E)
    return -rho.*sfun_nd(rho, rhoU, E)*entropy_scale
end

"function v_ufun(rho, rhou, rhov, E)
    Entropy variables as functions of conservative vars"
function v_ufun_nd(rho, rhoU, E)
    s = sfun_nd(rho,rhoU,E)
    p = pfun_nd(rho,rhoU,E)

    v1 = (@. (γ + 1 - s) - (γ-1)*E/p)
    vU = (x->@. x*(γ-1)/p).(rhoU)
    vE = (x->@. x*(γ-1)/p)(-rho)

    # v1,vU,vE = scale_entropy_vars_output(v1, vU, vE)
    return v1,vU,vE
end

#####
##### functions of entropy variables
#####
"function s_vfun(v1,vU,vE)
    entropy as function of entropy variables"
function s_vfun_nd(v1,vU,vE)
    vUnorm = unorm(vU)
    return @. γ - v1 + vUnorm/(2*vE)
end

"function rhoe_vfun(v1,vU,vE)
    specific energy as function of entropy variables"
function rhoe_vfun_nd(v1,vU,vE)
    s = s_vfun_nd(v1,vU,vE)
    return (@. ((γ-1)/((-vE)^γ))^(1/(γ-1)) * exp(-s/(γ-1)))
end

"function u_vfun(v1,vU,vE)
    Conservative vars as functions of entropy variables"
function u_vfun_nd(v1,vU,vE)
    # v1,vU,vE = scale_entropy_vars_input(v1,vU,vE)
    rhoeV     = rhoe_vfun_nd(v1,vU,vE)
    vUnorm    = unorm(vU)
    rho       = (@. rhoeV*(-vE))
    rhoU      = (x->rhoeV.*x).(vU)
    E         = (@. rhoeV*(1-vUnorm/(2*vE)))
    return (rho,rhoU,E)
end

function conservative_to_primitive_beta_nd(rho,rhoU,E)
    return rho, (x->x./rho).(rhoU), betafun_nd(rho,rhoU,E)
end
