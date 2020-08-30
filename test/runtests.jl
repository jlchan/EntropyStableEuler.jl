using EntropyStableEuler
using Test

@testset "Logmean tests" begin
    uL,uR = 1,2
    @test logmean(uL,uR) == logmean(uL,uR,log(uL),log(uR))
    @test logmean(uL,uR) == logmean(uR,uL)
    @test logmean(uL,uL) ≈ uL
end

@testset "2D entropy variable tests" begin
    using EntropyStableEuler.Fluxes2D
    import EntropyStableEuler.Fluxes2D: γ
    
    rho,u,v,p = 1,.1,.2,2
    rho,rhou,rhov,E = primitive_to_conservative(rho,u,v,p)
    v1,v2,v3,v4 = v_ufun(rho,rhou,rhov,E)

    h = 1e-7
    central_diff(f,x) = (f(x+h) - f(x-h))/(2*h)

    @test abs(v1 - central_diff(rho->Sfun(rho,rhou,rhov,E),rho)) < h
    @test abs(v2 - central_diff(rhou->Sfun(rho,rhou,rhov,E),rhou)) < h
    @test abs(v3 - central_diff(rhov->Sfun(rho,rhou,rhov,E),rhov)) < h
    @test abs(v4 - central_diff(E->Sfun(rho,rhou,rhov,E),E)) < h

    u1,u2,u3,u4 = u_vfun(v1,v2,v3,v4)
    @test u1 ≈ rho
    @test u2 ≈ rhou
    @test u3 ≈ rhov
    @test u4 ≈ E

    # test symmetry
    UL = copy.((rho,rhou,rhov,E))
    VL = copy.((v1,v2,v3,v4))
    UR = primitive_to_conservative(1.1,.2,.3,2.1)
    VR = v_ufun(UR...)
    QL = conservative_to_primitive_beta(UL...)
    QR = conservative_to_primitive_beta(UR...)
    Fx,Fy = euler_fluxes_2D(QL...,QR...)
    Fx2,Fy2 = euler_fluxes_2D(QR...,QL...)
    @test all(isapprox.(Fx,Fx2))
    @test all(isapprox.(Fy,Fy2))

    # test consistency
    p = Fluxes2D.pfun(rho,rhou,rhov,E)
    exact_flux_x = (rho*u, rho*u^2 + p, rhou*v, u*(E+p))
    exact_flux_y = (rho*v, rhou*v, rho*v^2 + p, v*(E+p))
    FxL,FyL = euler_fluxes_2D(QL...,QL...)
    @test all(isapprox.(FxL,exact_flux_x))
    @test all(isapprox.(FyL,exact_flux_y))

    # test entropy conservation property
    vTFx = sum(((x,y,z)->((x-y)*z)).(VL,VR,Fx))
    vTFy = sum(((x,y,z)->((x-y)*z)).(VL,VR,Fy))
    ψx(U) = (γ-1)*U[2]*Fluxes2D.entropy_scale # entropy potentials
    ψy(U) = (γ-1)*U[3]*Fluxes2D.entropy_scale
    @test vTFx ≈ ψx(UL)-ψx(UR)
    @test vTFy ≈ ψy(UL)-ψy(UR)
end

# module bar
# export f
# f(x) = f(x...) # dispatch
# module foo1
# import ..bar: f
# f(x1,x2) = x1+x2
# end
# module foo2
# import ..bar: f
# f(x1,x2,x3) = x1+x2+x3
# end
# end
