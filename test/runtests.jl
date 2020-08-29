using EntropyStableEuler
using Test

@testset "Logmean tests" begin
    uL,uR = 1,2
    @test logmean(uL,uR) == logmean(uL,uR,log(uL),log(uR))
    @test logmean(uL,uL) ≈ uL
end

@testset "Entropy variable tests" begin
    @test true    
end
