@testset "projector_function C4v E irrep" begin
    V = SU2Space(1//2=>1, 0=>1)
    P = SU2Space(1//2=>1)
    T = rand(ComplexF64, P, V^4)

    fproj = SpatiallySymmetricTensors.projector_function(C4v(), :E)
    Tp = fproj(T)
    @test norm(fproj(Tp) - Tp) < 1e-10
end
