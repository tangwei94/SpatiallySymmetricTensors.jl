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

@testset "test for matrix_for_spatial_operation" begin
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
        for _ in 1:6
            T0 = rand(ComplexF64, P, V^6)

            T1 = permute(T0, perm)
            @test norm(R_mat * vec(T0) - vec(T1)) < 1e-12
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

    f_identity(T1) = T1
    P_sol = find_subspace(T, P_init, f_identity; λ=2.0, is_hermitian=true, _mapping_table=mt)
    @test size(P_sol, 1) == num_paras
    @test size(P_sol, 2) == 0
   
    f_rotation = SpatiallySymmetricTensors.linear_function_for_spatial_operation(((1, ), (3, 4, 5, 2)))
    P_sol = find_subspace(T, P_init, f_rotation; λ=2.0, is_hermitian=false, _mapping_table=mt)
    @test size(P_sol, 1) == num_paras
    @test size(P_sol, 2) == 0
end

@testset "non-orthonormal tolerance for P_init" begin
    # orthonormality check; should respect the tolerance parameter
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

@testset "test selector: identity" begin
    # selector with a true condition should be the identity projector
    V = SU2Space(1//2=>1, 0=>1)
    P = SU2Space(1//2=>1)
    T = zeros(ComplexF64, P, V^4)

    mt = mapping_table(T)
    num_paras = num_free_parameters(T; _mapping_table=mt)
    Psel = selector(T, (f1, f2) -> true; _mapping_table=mt)
    @test Psel == Matrix{ComplexF64}(I, num_paras, num_paras)
end

@testset "test selector: fusion tree condition" begin
    V = SU2Space(1//2=>1, 0=>1)
    P = SU2Space(1//2=>1)
    T = zeros(ComplexF64, P, V^4)

    mt = mapping_table(T)
    num_paras = num_free_parameters(T; _mapping_table=mt)

    # for the given P and V, there are only two possiblities
    # looking for the subspace with (1//2 ⊕ 0 ⊕ 0 ⊕ 0 -> 1//2) 
    condition_S(f1, f2) = length(findall(rep -> rep == SU2Irrep(1//2), f2.uncoupled)) == 1
    # looking for the subspace with (1//2 ⊕ 1//2 ⊕ 1//2 ⊕ 0 -> 1//2) 
    condition_L(f1, f2) = length(findall(rep -> rep == SU2Irrep(1//2), f2.uncoupled)) == 3

    Psel_S = selector(T, condition_S; _mapping_table=mt)
    Psel_L = selector(T, condition_L; _mapping_table=mt)

    @test size(Psel_S, 1) == num_paras
    @test size(Psel_S, 2) + size(Psel_L, 2) == num_paras
    @test norm(Psel_S' * Psel_L) < 1e-12
    @test rank(hcat(Psel_S, Psel_L)) == num_paras
end
