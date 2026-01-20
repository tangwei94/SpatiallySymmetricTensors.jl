@testset "test set_data_by_vector and vec" begin 
    # round-trip between vectorized parameters and symmetric tensor storage
    V = SU2Space(1//2=>1, 0=>1);
    P = SU2Space(1//2=>1);
    T = zeros(ComplexF64, P, V^4);

    mt = mapping_table(T)
    num_paras = num_free_parameters(T; _mapping_table=mt)

    v = rand(num_paras)

    T1 = set_data_by_vector(T, v; _mapping_table=mt)
    @test norm(vec(T1) - v) < 1e-12

    # in-place setter returns the same tensor and mutates as expected
    T2 = similar(T)
    T2_ref = set_data_by_vector!(T2, v; _mapping_table=mt)
    @test T2_ref === T2
    @test norm(vec(T2) - v) < 1e-12
end

@testset "spatial_operations.jl" begin #for _ in 1:3
    V = SU2Space(1//2=>1, 0=>1)
    P = SU2Space(1//2=>1)
    T = zeros(ComplexF64, P, V^6)

    permutations = [((1, ), (3, 4, 5, 6, 7, 2)), 
                    ((1, ), (4, 5, 6, 7, 2, 3)),
                    ((1, ), (5, 6, 7, 2, 3, 4)),
                    ((1, ), (6, 7, 2, 3, 4, 5)),
                    ((1, ), (7, 2, 3, 4, 5, 6)),
                    ((1, ), (4, 3, 6, 2, 7, 5))]

    for perm in permutations
        R_mat = matrix_for_spatial_operation(T, perm)
        num_paras = num_free_parameters(T)
        for ix0 in 1:num_paras
            T0 = begin
                paras = zeros(num_paras)
                paras[ix0] = 1
                set_data_by_vector(T, paras)
            end
        
            T1 = set_data_by_vector(T, R_mat[:, ix0])
            @test norm(permute(T0, perm) - T1) < 1e-12
        end
    end
end

@testset "empty eigenspace" begin
    # no eigenvalue match should return an empty subspace (hermitian case)
    V = SU2Space(1//2=>1, 0=>1)
    P = SU2Space(1//2=>1)
    T = zeros(ComplexF64, P, V^4)

    mt = mapping_table(T)
    num_paras = num_free_parameters(T; _mapping_table=mt)
    P_init = Matrix{ComplexF64}(I, num_paras, num_paras)

    f_op(T1) = T1
    P_sol = find_subspace(T, P_init, f_op; λ=2.0, is_hermitian=true, _mapping_table=mt)
    @test size(P_sol, 1) == num_paras
    @test size(P_sol, 2) == 0
end

@testset "empty eigenspace (non-hermitian)" begin
    # same as above, but for non-hermitian branch
    V = SU2Space(1//2=>1, 0=>1)
    P = SU2Space(1//2=>1)
    T = zeros(ComplexF64, P, V^4)

    mt = mapping_table(T)
    num_paras = num_free_parameters(T; _mapping_table=mt)
    P_init = Matrix{ComplexF64}(I, num_paras, num_paras)

    f_op(T1) = T1
    P_sol = find_subspace(T, P_init, f_op; λ=2.0, is_hermitian=false, _mapping_table=mt)
    @test size(P_sol, 1) == num_paras
    @test size(P_sol, 2) == 0
end

@testset "invalid point-group reps" begin
    # unknown representation names should throw for each point group
    V = SU2Space(1//2=>1, 0=>1)
    P = SU2Space(1//2=>1)
    T = zeros(ComplexF64, P, V^4)

    @test_throws ArgumentError find_subspace(C3v(), T, :BadRep)
    @test_throws ArgumentError find_subspace(C4v(), T, :BadRep)
    @test_throws ArgumentError find_subspace(C6v(), T, :BadRep)
    @test_throws ArgumentError find_subspace(C4(), T, :BadRep)
    @test_throws ArgumentError find_subspace(D2(), T, :BadRep)
end

@testset "invalid point-group perm ops" begin
    # unknown operation names should throw for each point group
    for g in (C3v(), C4v(), C6v(), C4(), D2())
        @test_throws ArgumentError SpatiallySymmetricTensors.get_perm(g, :BadOp)
    end
end

@testset "pointgroup verbosity" begin
    # default verbosity is silent and should return a valid subspace
    V = SU2Space(1//2=>1, 0=>1)
    P = SU2Space(1//2=>1)
    T = zeros(ComplexF64, P, V^4)

    P_sol = find_subspace(C4v(), T, :A1; verbose=false)
    @test size(P_sol, 1) == num_free_parameters(T)
end

@testset "non-orthonormal tolerance" begin
    # orthonormality check should respect the tolerance parameter
    V = SU2Space(1//2=>1, 0=>1)
    P = SU2Space(1//2=>1)
    T = zeros(ComplexF64, P, V^4)

    mt = mapping_table(T)
    num_paras = num_free_parameters(T; _mapping_table=mt)
    P_init = Matrix{ComplexF64}(I, num_paras, num_paras)
    P_init[1, 1] += 1e-3

    f_op(T1) = T1
    @test_throws ErrorException find_subspace(T, P_init, f_op; _mapping_table=mt)
    P_sol = find_subspace(T, P_init, f_op; tol=1e-2, _mapping_table=mt)
    @test size(P_sol, 1) == num_paras
end

@testset "selector identity" begin
    # selector with a true condition should be the identity projector
    V = SU2Space(1//2=>1, 0=>1)
    P = SU2Space(1//2=>1)
    T = zeros(ComplexF64, P, V^4)

    mt = mapping_table(T)
    num_paras = num_free_parameters(T; _mapping_table=mt)
    Psel = selector(T, (f1, f2) -> true; _mapping_table=mt)
    @test Psel == Matrix{ComplexF64}(I, num_paras, num_paras)
end

@testset "selector columns" begin
    # selector should match the explicit column subset for a condition
    V = SU2Space(1//2=>1, 0=>1)
    P = SU2Space(1//2=>1)
    T = zeros(ComplexF64, P, V^4)

    mt = mapping_table(T)
    num_paras = num_free_parameters(T; _mapping_table=mt)
    condition(f1, f2) = length(findall(rep -> rep == SU2Irrep(1//2), f2.uncoupled)) == 1

    indices = Int[]
    for (f1, f2, a, n) in mt
        if condition(f1, f2)
            append!(indices, a:a+n-1)
        end
    end

    Psel = selector(T, condition; _mapping_table=mt)
    @test size(Psel, 1) == num_paras
    @test size(Psel, 2) == length(indices)
    @test Psel == Matrix{ComplexF64}(I, num_paras, num_paras)[:, indices]
end

@testset "pointgroup get_perm" begin
    # basic sanity: counts and permutation structure per point group/operator
    groups = [
        (C3v(), Dict(:σd=>0, :σv=>3, :R=>2), 4),
        (C4v(), Dict(:σd=>2, :σv=>2, :R=>2), 5),
        (C6v(), Dict(:σd=>3, :σv=>3, :R=>2), 7),
        (C4(), Dict(:σd=>0, :σv=>0, :R=>2), 5),
        (D2(), Dict(:σd=>0, :σv=>2, :R=>0), 5),
    ]
    for (g, counts, nlegs) in groups
        for (op, expected) in counts
            perms = SpatiallySymmetricTensors.get_perm(g, op)
            @test length(perms) == expected
            for perm in perms
                @test perm[1] == (1,)
                @test sort(collect(perm[2])) == collect(2:nlegs)
            end
        end
    end
end

@testset "pointgroup perm algebra" begin
    # check involution for reflections and group order/inverse for rotations
    function compose_perm(p1, p2)
        a = p1[2]
        b = p2[2]
        composed = ntuple(i -> a[b[i] - 1], length(a))
        return (p1[1], composed)
    end

    function is_identity_perm(p, nlegs)
        return p[1] == (1,) && p[2] == ntuple(i -> i + 1, nlegs - 1)
    end

    function perm_power(p, k)
        res = p
        for _ in 2:k
            res = compose_perm(res, p)
        end
        return res
    end

    for g in (C3v(), C4v(), C6v(), D2(), C4())
        nlegs = g isa C6v ? 7 : (g isa C3v ? 4 : 5)
        for op in (:σd, :σv)
            for p in SpatiallySymmetricTensors.get_perm(g, op)
                @test is_identity_perm(compose_perm(p, p), nlegs)
            end
        end
    end

    for p in SpatiallySymmetricTensors.get_perm(C4v(), :R)
        @test is_identity_perm(perm_power(p, 4), 5)
    end
    for p in SpatiallySymmetricTensors.get_perm(C6v(), :R)
        @test is_identity_perm(perm_power(p, 6), 7)
    end
    for p in SpatiallySymmetricTensors.get_perm(C3v(), :R)
        @test is_identity_perm(perm_power(p, 3), 4)
    end

    c4v_r = SpatiallySymmetricTensors.get_perm(C4v(), :R)
    @test is_identity_perm(compose_perm(c4v_r[1], c4v_r[2]), 5)
    @test is_identity_perm(compose_perm(c4v_r[2], c4v_r[1]), 5)

    c6v_r = SpatiallySymmetricTensors.get_perm(C6v(), :R)
    @test is_identity_perm(compose_perm(c6v_r[1], c6v_r[2]), 7)
    @test is_identity_perm(compose_perm(c6v_r[2], c6v_r[1]), 7)

    c3v_r = SpatiallySymmetricTensors.get_perm(C3v(), :R)
    @test is_identity_perm(compose_perm(c3v_r[1], c3v_r[2]), 4)
    @test is_identity_perm(compose_perm(c3v_r[2], c3v_r[1]), 4)
end

@testset "C4v E irrep" begin
    P = ℂ^2
    V = ℂ^2
    T = zeros(ComplexF64, P, V^4)

    mt = mapping_table(T)
    P_E = find_subspace(C4v(), T, :E)
    @test size(P_E, 2) % 2 == 0
    nblocks = size(P_E, 2) ÷ 2
    @test nblocks > 0

    reps = SpatiallySymmetricTensors.get_reps(C4v(), :E)
    ops = [(:σd, reps[1]), (:σv, reps[2]), (:R, reps[3])]
    for (op, rep_list) in ops
        perms = SpatiallySymmetricTensors.get_perm(C4v(), op)
        @test length(perms) == length(rep_list)
        for (perm, D) in zip(perms, rep_list)
            M = matrix_for_spatial_operation(T, perm; _mapping_table=mt)
            for b in 1:nblocks
                cols = P_E[:, (b - 1) * 2 + 1:b * 2]
                @test norm(M * cols - cols * D) < 1e-10
            end
        end
    end
end

@testset "C4v E irrep with U1GradedSpace" begin
    P = U1Space(0=>1, 1=>1, -1=>1)
    V = U1Space(0=>3, 1=>1, -1=>1)
    T = zeros(ComplexF64, P, V^4)

    mt = mapping_table(T)
    P_E = find_subspace(C4v(), T, :E)
    @test size(P_E, 2) % 2 == 0
    @test size(P_E, 2) > 0
    reps = SpatiallySymmetricTensors.get_reps(C4v(), :E)
    ops = [(:σd, reps[1]), (:σv, reps[2]), (:R, reps[3])]
    for (op, rep_list) in ops
        perms = SpatiallySymmetricTensors.get_perm(C4v(), op)
        @test length(perms) == length(rep_list)
        for perm in perms
            M = matrix_for_spatial_operation(T, perm; _mapping_table=mt)
            MC = M * P_E
            MC_proj = P_E * (P_E \ MC)
            @test norm(MC - MC_proj) / max(norm(MC), eps()) < 1e-2
        end
    end
end
