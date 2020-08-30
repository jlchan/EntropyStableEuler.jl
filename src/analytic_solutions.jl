"
function vortex(x, y, t; x0=5, y0=0)

2D isentropic vortex solution, centered at x0,y0 = (5,0) by default.
Assumes domain around [0,20] x [-5,5].
Returns variables (rho, u, v, p)
"
function vortex(x, y, t; x0=5, y0=0)
    beta = 5
    r2 = @. (x-x0-t)^2 + (y-y0)^2

    u = @. 1 - beta*exp(1-r2)*(y-y0)/(2*pi)
    v = @. beta*exp(1-r2)*(x-x0-t)/(2*pi)
    rho = @. 1 - (1/(8*γ*pi^2))*(γ-1)/2*(beta*exp(1-r2))^2
    rho = @. rho^(1/(γ-1))
    p = @. rho^γ

    return (rho, u, v, p)
end

"
function vortex(x, y, z, t)

3D isentropic vortex 'column' solution, centered at x0,y0,z0 = (7.5,7.5,0) by default.
Usually solved on the domain [0,15] x [0,20] x [-5,5].
Returns variables (rho, u, v, w, p)
"
function vortex(x, y, z, t)
    p0 = 1/γ
    Πmax = .4
    c = 7.5
    r1 = @. (y - c - t)
    r2 = @. (x - c)
    r3 = zeros(size(x))
    Π = @. Πmax*exp(.5*(1 - (r1^2 + r2^2 + r3^2)))

    rho = @. 1 - .5*(γ-1)*Π^2
    rho = @. rho^(1/(γ-1))
    u   = @. Π*r1
    v   = @. Π*r2 + 1
    w   = @. Π*r3
    Unorm = @. u^2 + v^2 + w^2
    E = @. (p0/(γ-1) * ρ^γ) + .5*ρ*Unorm

    return rho, u, v, w, p
end
