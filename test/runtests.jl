using EntropyStableEuler
using Test

@testset "Logmean tests" begin
    uL,uR = 1,2
    @test logmean(uL,uR) == logmean(uL,uR,log(uL),log(uR))
    @test logmean(uL,uR) == logmean(uR,uL)
    @test logmean(uL,uL) ≈ uL
end

@testset "Entropy variable tests" begin
    @testset "1D tests" begin
        # @test true
    end

    @testset "2D tests" begin
        rho,u,v,p = 1,.1,.2,2
        rho,rhou,rhov,E = primitive_to_conservative(rho,u,v,p)
        v1,v2,v3,v4 = v_ufun(rho,rhou,rhov,E)

        h = 1e-7
        central_diff(f,x) = (f(x+h) - f(x-h))/(2*h)

        # @show [abs(v1 - central_diff(rho->Sfun(rho,rhou,rhov,E),rho))
        # abs(v2 - central_diff(rhou->Sfun(rho,rhou,rhov,E),rhou))
        # abs(v3 - central_diff(rhov->Sfun(rho,rhou,rhov,E),rhov))
        # abs(v4 - central_diff(E->Sfun(rho,rhou,rhov,E),E))]

        @test abs(v1 - central_diff(rho->Sfun(rho,rhou,rhov,E),rho)) < h
        @test abs(v2 - central_diff(rhou->Sfun(rho,rhou,rhov,E),rhou)) < h
        @test abs(v3 - central_diff(rhov->Sfun(rho,rhou,rhov,E),rhov)) < h
        @test abs(v4 - central_diff(E->Sfun(rho,rhou,rhov,E),E)) < h

        u1,u2,u3,u4 = u_vfun(v1,v2,v3,v4)
        @test u1 ≈ rho
        @test u2 ≈ rhou
        @test u3 ≈ rhov
        @test u4 ≈ E
    end
end
