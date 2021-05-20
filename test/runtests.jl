using EntropyStableEuler
using Test
using StaticArrays

norm(u) = sqrt(sum(u.^2))

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

@testset "Tests for d = $d" for d in (1:3)

    @testset "Entropy variable tests" begin

        U = prim_to_cons(Euler{d}(),init_prim(d))
        V = cons_to_entropy(Euler{d}(),U)

        h = 1e-7
        central_diff(f,x) = (f(x+h) - f(x-h))/(2*h)
        swapentry(x,y,i) = (x[begin:i-1]...,y,x[i+1:end]...)
        for j = 1:d+2
            @test abs(V[j] - central_diff(x->Sfun(Euler{d}(),swapentry(U,x,j)),U[j])) < h
        end

        UV = entropy_to_cons(Euler{d}(),V)
        @test all(UV .≈ U)
    end

    @testset "Symmetry" begin
        # test symmetry
        UL = prim_to_cons(Euler{d}(),init_prim(d))
        UR = prim_to_cons(Euler{d}(),init_prim(d).*1.1)
        FLR = fS(Euler{d}(),UL,UR)
        FRL = fS(Euler{d}(),UR,UL)
        @test all(FLR .≈ FRL)
    end

    @testset "Consistency" begin
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

    @testset "Entropy conservation property" begin
        γ = Euler{d}().γ

        UL = prim_to_cons(Euler{d}(),init_prim(d))
        UR = prim_to_cons(Euler{d}(),init_prim(d).*1.1)
        VL = cons_to_entropy(Euler{d}(),UL)
        VR = cons_to_entropy(Euler{d}(),UR)

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

    @testset "Dissipation" begin
        UL = prim_to_cons(Euler{d}(),init_prim(d))
        UR = prim_to_cons(Euler{d}(),init_prim(d).*1.1)
        VL = cons_to_entropy(Euler{d}(),UL)
        VR = cons_to_entropy(Euler{d}(),UR)

        normal = randn(d)
        normal = normal/norm(normal)
        dissipation = LxF_dissipation(Euler{d}(),normal,UL,UL)
        @test norm(dissipation) < 50*eps()

        dissipation = LxF_dissipation(Euler{d}(),normal,UL,UR)
        vTF(VL,VR,F) = sum((VL .- VR).*F)
        @test vTF(VL,VR,dissipation) > 0.0
    end

    @testset "Type stability tests" begin
        Q = init_prim(d)
        U = prim_to_cons(Euler{d}(),Q)
        V = cons_to_entropy(Euler{d}(),U)

        @inferred prim_to_cons(Euler{d}(),Q)
        @inferred cons_to_prim_beta(Euler{d}(),U)
        @inferred cons_to_entropy(Euler{d}(),U)
        @inferred entropy_to_cons(Euler{d}(),V)
        @inferred fS(Euler{d}(),U,U)
    end
end
