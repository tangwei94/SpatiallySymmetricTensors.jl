@testset "projector_function" begin
    V = SU2Space(1//2=>1, 0=>1)
    P = SU2Space(1//2=>1)
    T = rand(ComplexF64, P, V^4)

    fproj = SpatiallySymmetricTensors.projector_function(C4v(), :B1)
    reps = SpatiallySymmetricTensors.get_reps(C4v(), :B1)
    op_names = (:σd, :σv, :R)

    Tp = fproj(T)

    # projected tensor should transform under the target irrep
    for (i, name) in enumerate(op_names)
        χ = reps[i]
        χ == 0 && continue
        for perm in SpatiallySymmetricTensors.get_perm(C4v(), name)
            #@test norm(permute(Tp, perm) - χ * Tp) < 1e-10
        end
    end

    # projector should be idempotent
    @show norm(fproj(fproj(T)) - Tp)
    @show norm(T)
    @show norm(fproj(T))
    @show norm(fproj(fproj(T)))
#    @test norm(fproj(fproj(T)) - Tp) < 1e-10
end
