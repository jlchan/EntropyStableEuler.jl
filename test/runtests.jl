using EntropyStableEuler
using Test
using StaticArrays

@testset "Logmean tests" begin
    uL,uR = 1,2
    @test logmean(uL,uR) == logmean(uL,uR,log(uL),log(uR))
    @test logmean(uL,uR) == logmean(uR,uL)
    @test logmean(uL,uL) ≈ uL
end

function init_prim(d)
    if d==1
        rho,u,p = 1,.1,2
        Q = (rho,u,p)
    elseif d==2
        rho,u,v,p = 1,.1,.2,2
        Q = (rho,u,v,p)
    elseif d==3
        rho,u,v,w,p = 1,.1,.2,.3,2
        Q = (rho,u,v,w,p)
    end
    return Q
end

@testset "Entropy variable tests for d=$d" for d in (1:3)

    U = prim_to_cons(Euler{d}(),init_prim(d))
    V = v_ufun(Euler{d}(),U)

    h = 1e-7
    central_diff(f,x) = (f(x+h) - f(x-h))/(2*h)
    swapentry(x,y,i) = (x[begin:i-1]...,y,x[i+1:end]...)
    for j = 1:d+2
        @test abs(V[j] - central_diff(x->Sfun(Euler{d}(),swapentry(U,x,j)),U[j])) < h
    end

    UV = u_vfun(Euler{d}(),V)
    @test all(UV .≈ U)
end

@testset "Symmetry of entropy conservative fluxes for d=$d" for d in (1:3)
    # test symmetry
    UL = prim_to_cons(Euler{d}(),init_prim(d))
    UR = prim_to_cons(Euler{d}(),init_prim(d).*1.1)
    FLR = fS(Euler{d}(),UL,UR)
    FRL = fS(Euler{d}(),UR,UL)
    @test all(FLR .≈ FRL)
end

@testset "Consistency of entropy conservative fluxes for d=$d" for d in (1:3)
    Q = init_prim(d)
    U = prim_to_cons(Euler{d}(),Q)
    F = fS(Euler{d}(),U,U)

    p = pfun(Euler{d}(),U)
    if d==1
        rho,rhou,E = U
        u = rhou./rho
        Fexact = SVector{3}(rho*u, rho*u^2 + p, u*(E+p))
    elseif d==2
        rho,rhou,rhov,E = U
        u,v = rhou./rho, rhov./rho
        Fx = SVector{4}(rho*u, rho*u^2 + p, rho*u*v,     u*(E+p))
        Fy = SVector{4}(rho*v, rho*u*v,     rho*v^2 + p, v*(E+p))
        Fexact = (Fx,Fy)
    elseif d==3
        rho,rhou,rhov,rhow,E = U
        u,v,w = rhou./rho, rhov./rho, rhow./rho
        Fx = SVector{5}(rho*u, rho*u^2 + p, rho*u*v,     rho*u*w,     u*(E+p))
        Fy = SVector{5}(rho*v, rho*u*v,     rho*v^2 + p, rho*v*w,     v*(E+p))
        Fz = SVector{5}(rho*w, rho*u*w,     rho*v*w,     rho*w^2 + p, w*(E+p))
        Fexact = (Fx,Fy,Fz)
    end

    @test all(F .≈ Fexact)
end

@testset "Entropy conservation property for d=$d" for d in (1:3)
    γ = Euler{d}().γ

    UL = prim_to_cons(Euler{d}(),init_prim(d))
    UR = prim_to_cons(Euler{d}(),init_prim(d).*1.1)
    VL = v_ufun(Euler{d}(),UL)
    VR = v_ufun(Euler{d}(),UR)

    ψ(U) = (γ-1).*U[2:d+1]
    F = fS(Euler{d}(),UL,UR)
    vTF(VL,VR,F) = sum((VL .- VR).*F)
    for j = 1:d
        if d==1
            @test vTF(VL,VR,F) ≈ ψ(UL)[j]-ψ(UR)[j]
        else
            @test vTF(VL,VR,F[j]) ≈ ψ(UL)[j]-ψ(UR)[j]
        end
    end
end

@testset "Check type stability for d=$d" for d in (1:3)
    Q = init_prim(d)
    U = prim_to_cons(Euler{d}(),Q)
    V = v_ufun(Euler{d}(),U)

    @inferred prim_to_cons(Euler{d}(),Q)
    @inferred cons_to_prim_beta(Euler{d}(),U)
    @inferred v_ufun(Euler{d}(),U)
    @inferred u_vfun(Euler{d}(),V)
    @inferred fS(Euler{d}(),U,U)
end
