using ABCDBeamTrace
@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

@testset "beamtrace" begin
    f = 100e-3
    L = 1000e-3
    expander_2x = [ThinLens(f), FreeSpace(3f), ThinLens(2f)]
    system = [expander_2x; FreeSpace(L); reverse(expander_2x)]
    beam = Beam(1000e-9, 1.0e-3)
    @test beamtrace(system, beam) isa Vector{Beam}
end

@testset "discretize" begin
    @test_broken discretize(FreeSpace(100), 2) ==
        [FreeSpace(50), FreeSpace(50)]
    @test_broken discretize([FreeSpace(100)], 2) ==
        [FreeSpace(50), FreeSpace(50)]
end
