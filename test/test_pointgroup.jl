function _check_projector(group, irrep_name, T)
    fproj = SpatiallySymmetricTensors.projector_function(group, irrep_name)
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

@testset "projector_function C4 C4v, SU2Space" begin
    V = SU2Space(1//2=>1, 0=>1)
    P = SU2Space(1//2=>1)

    for (group, irrep_name) in [
        (C4v(), :A1), (C4v(), :A2), (C4v(), :B1), (C4v(), :B2),
        (C4(), :A), (C4(), :B),
    ]
        _check_projector(group, irrep_name, rand(ComplexF64, P, V^4))
    end
end

@testset "projector_function C4 C4v, U1Space" begin
    V = U1Space(-1=>1, 0=>1, 1=>1)
    P = U1Space(0=>1, 1=>1, 2=>1)

    for (group, irrep_name) in [
        (C4v(), :A1), (C4v(), :A2), (C4v(), :B1), (C4v(), :B2),
        (C4(), :A), (C4(), :B),
    ]
        _check_projector(group, irrep_name, rand(ComplexF64, P, V^4))
    end
end

@testset "projector_function C3v, SU2Space" begin
    V = SU2Space(1//2=>1, 0=>1)
    P = SU2Space(1//2=>1)

    for (group, irrep_name) in [
        (C3v(), :A1), (C3v(), :A2),
    ]
        _check_projector(group, irrep_name, rand(ComplexF64, P, V^3))
    end
end

@testset "projector_function C3v, U1Space" begin
    V = U1Space(-1=>1, 0=>1, 1=>1)
    P = U1Space(0=>1, 1=>1, 2=>1)

    for (group, irrep_name) in [
        (C3v(), :A1), (C3v(), :A2),
    ]
        _check_projector(group, irrep_name, rand(ComplexF64, P, V^3))
    end
end

@testset "projector_function C6v, SU2Space" begin
    V = SU2Space(1//2=>1, 0=>1)
    P = SU2Space(1//2=>1)

    for (group, irrep_name) in [
        (C6v(), :A1), (C6v(), :A2), (C6v(), :B1), (C6v(), :B2),
    ]
        _check_projector(group, irrep_name, rand(ComplexF64, P, V^6))
    end
end

@testset "projector_function C6v, U1Space" begin
    V = U1Space(-1=>1, 0=>1, 1=>1)
    P = U1Space(0=>1, 1=>1, 2=>1)

    for (group, irrep_name) in [
        (C6v(), :A1), (C6v(), :A2), (C6v(), :B1), (C6v(), :B2),
    ]
        _check_projector(group, irrep_name, rand(ComplexF64, P, V^6))
    end
end

@testset "projector_function D2, SU2Space" begin
    V = SU2Space(1//2=>1, 0=>1)
    P = SU2Space(1//2=>1)

    for (group, irrep_name) in [
        (D2(), :A), (D2(), :B1), (D2(), :B2), (D2(), :B3),
    ]
        _check_projector(group, irrep_name, rand(ComplexF64, P, V^4))
    end
end

@testset "projector_function D2, U1Space" begin
    V = U1Space(-1=>1, 0=>1, 1=>1)
    P = U1Space(0=>1, 1=>1, 2=>1)

    for (group, irrep_name) in [
        (D2(), :A), (D2(), :B1), (D2(), :B2), (D2(), :B3),
    ]
        _check_projector(group, irrep_name, rand(ComplexF64, P, V^4))
    end
end

@testset "find_subspace and find_solution C4v, SU2Space" begin
    V = SU2Space(1//2=>1, 0=>1)
    P = SU2Space(1//2=>1)
    T = zeros(ComplexF64, P, V^4)

    group = C4v()
    for irrep_name in (:A1, :B1)
        P_sol = find_subspace(group, T, irrep_name)
        if size(P_sol, 2) > 0
            @test norm(P_sol' * P_sol - Matrix{ComplexF64}(I, size(P_sol, 2), size(P_sol, 2))) < 1e-10
        end

        sols = find_solution(group, T, irrep_name)
        fproj = SpatiallySymmetricTensors.projector_function(group, irrep_name)
        ops = SpatiallySymmetricTensors.group_elements(group)
        irrep_chars = SpatiallySymmetricTensors.irrep_chars(group, irrep_name)
        for sol in sols
            @test abs(norm(sol) - 1) < 1e-10
            @test norm(fproj(sol) - sol) < 1e-10
            for (name, perm) in ops
                χ = get(irrep_chars, name, nothing)
                @test χ !== nothing
                @test norm(permute(sol, perm) - χ * sol) < 1e-10
            end
        end
    end
end
