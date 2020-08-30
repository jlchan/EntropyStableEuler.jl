"function primitive_to_conservative(rho,u,v,p)

    convert primitive variables (ρ,u,v,p) to conservative vars (ρ,ρu,ρv,E)."
function primitive_to_conservative(rho,u,v,p)
    unorm = @. u^2 + v^2
    rhou,rhov = (x->rho.*x).((u,v))
    E = (@. p/(γ-1) + .5*rho*unorm)
    return (rho,rhou,rhov,E)
end

#####
##### functions of conservative variables
#####

"function pfun(rho, rhou, rhov, E)
    pressure as a function of conservative variables"
function pfun(rho, rhou, rhov, E)
    rhoUnorm2 = @. (rhou^2 + rhov^2)/rho
    return @. (γ-1)*(E - .5*rhoUnorm2)
end

"function rhoe_ufun(rho, rhou, rhov, E)
    specific energy as a function of conservative variables"
function rhoe_ufun(rho, rhou, rhov, E)
    return pfun(rho, rhou, rhov, E) / (γ-1)
end

"function sfun(rho, rhou, rhov, E)
    Specific entropy as a function of conservative variables"
function sfun(rho, rhou, rhov, E)
    p = pfun(rho, rhou, rhov, E)
    return @. log(p/(rho^γ))
end

"function Sfun(rho,rhou,rhov,E)
    Mathematical entropy as a function of conservative variables"
function Sfun(rho, rhou, rhov, E)
    # return -rho.*sfun(rho, rhou, rhov, E)
    return -rho.*sfun(rho, rhou, rhov, E)*entropy_scale
end

const entropy_scale = 1/(γ-1)
scale_entropy_vars_output(V...) = (x->x*entropy_scale).(V)
scale_entropy_vars_input(V...) = (x->x/entropy_scale).(V)


"function v_ufun(rho, rhou, rhov, E)
    Entropy variables as functions of conservative vars"
function v_ufun(rho, rhou, rhov, E)
    s = sfun(rho,rhou,rhov,E)
    p = pfun(rho,rhou,rhov,E)

    v1 = (@. (γ + 1 - s) - (γ-1)*E/p)
    vU1,vU2,vE = (x->@. x*(γ-1)/p).((rhou,rhov,-rho))

    # v1 = (@. (γ + 1 - s)/(γ-1) - E/p)
    # vU1,vU2,vE = (x->x./p).((rhou,rhov,-rho))

    return scale_entropy_vars_output(v1, vU1, vU2, vE)
end

#####
##### functions of entropy variables
#####

"function s_vfun(v1,vU1,vU2,vE)
    entropy as function of entropy variables"
function s_vfun(v1,vU1,vU2,vE)
    vUnorm = @. vU1.^2 + vU2.^2
    return @. γ - v1 + vUnorm/(2*vE)
end

"function rhoe_vfun(v1,vU1,vU2,vE)
    specific energy as function of entropy variables"
function rhoe_vfun(v1,vU1,vU2,vE)
    s = s_vfun(v1,vU1,vU2,vE)
    return (@. ((γ-1)/((-vE)^γ))^(1/(γ-1)) * exp(-s/(γ-1)))
end

"function u_vfun(v1,vU1,vU2,vE)
    Conservative vars as functions of entropy variables"
function u_vfun(v1,vU1,vU2,vE)
    v1,vU1,vU2,vE = scale_entropy_vars_input(v1,vU1,vU2,vE)
    rhoeV     = rhoe_vfun(v1,vU1,vU2,vE)
    vUnorm    = @. vU1.^2 + vU2.^2
    rho       = (@. rhoeV*(-vE))
    rhou,rhov = (x->rhoeV.*x).((vU1,vU2))
    E         = (@. rhoeV*(1-vUnorm/(2*vE)))
    return (rho,rhou,rhov,E)
end
