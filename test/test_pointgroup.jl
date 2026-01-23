@testset "projector_function C4 C4v, SU2Space" begin
    V = SU2Space(1//2=>1, 0=>1)
    P = SU2Space(1//2=>1)

    for (group, irrep_name) in [
        (C4v(), :A1), (C4v(), :A2), (C4v(), :B1), (C4v(), :B2),
        (C4(), :A), (C4(), :B),
    ]
        fproj = SpatiallySymmetricTensors.projector_function(group, irrep_name)

        T = rand(ComplexF64, P, V^4)
        Tp = fproj(T)
    
        # projector should be idempotent
        @test norm(fproj(Tp) - Tp) < 1e-10

        # projected tensor should transform under the target irrep
        ops = SpatiallySymmetricTensors.group_elements(group)
        reps = SpatiallySymmetricTensors.irrep_chars(group, irrep_name)
        for (name, perm) in ops
            χ = get(reps, name, nothing)
            @test χ !== nothing # check that the representation value exists
            @test norm(permute(Tp, perm) - χ * Tp) < 1e-10
        end
    end
end

@testset "projector_function C3v, SU2Space" begin
    V = SU2Space(1//2=>1, 0=>1)
    P = SU2Space(1//2=>1)

    for (group, irrep_name) in [
        (C3v(), :A1), (C3v(), :A2),
    ]
        fproj = SpatiallySymmetricTensors.projector_function(group, irrep_name)

        T = rand(ComplexF64, P, V^3)
        Tp = fproj(T)

        # projector should be idempotent
        @test norm(fproj(Tp) - Tp) < 1e-10

        # projected tensor should transform under the target irrep
        ops = SpatiallySymmetricTensors.group_elements(group)
        reps = SpatiallySymmetricTensors.irrep_chars(group, irrep_name)
        for (name, perm) in ops
            χ = get(reps, name, nothing)
            @test χ !== nothing
            @test norm(permute(Tp, perm) - χ * Tp) < 1e-10
        end
    end
end

@testset "projector_function C6v, SU2Space" begin
    V = SU2Space(1//2=>1, 0=>1)
    P = SU2Space(1//2=>1)

    for (group, irrep_name) in [
        (C6v(), :A1), (C6v(), :A2), (C6v(), :B1), (C6v(), :B2),
    ]
        fproj = SpatiallySymmetricTensors.projector_function(group, irrep_name)

        T = rand(ComplexF64, P, V^6)
        Tp = fproj(T)

        # projector should be idempotent
        @test norm(fproj(Tp) - Tp) < 1e-10

        # projected tensor should transform under the target irrep
        ops = SpatiallySymmetricTensors.group_elements(group)
        reps = SpatiallySymmetricTensors.irrep_chars(group, irrep_name)
        for (name, perm) in ops
            χ = get(reps, name, nothing)
            @test χ !== nothing
            @test norm(permute(Tp, perm) - χ * Tp) < 1e-10
        end
    end
end

@testset "projector_function D2, SU2Space" begin
    V = SU2Space(1//2=>1, 0=>1)
    P = SU2Space(1//2=>1)

    for (group, irrep_name) in [
        (D2(), :A), (D2(), :B1), (D2(), :B2), (D2(), :B3),
    ]
        fproj = SpatiallySymmetricTensors.projector_function(group, irrep_name)

        T = rand(ComplexF64, P, V^4)
        Tp = fproj(T)

        # projector should be idempotent
        @test norm(fproj(Tp) - Tp) < 1e-10

        # projected tensor should transform under the target irrep
        ops = SpatiallySymmetricTensors.group_elements(group)
        reps = SpatiallySymmetricTensors.irrep_chars(group, irrep_name)
        for (name, perm) in ops
            χ = get(reps, name, nothing)
            @test χ !== nothing
            @test norm(permute(Tp, perm) - χ * Tp) < 1e-10
        end
    end
end
