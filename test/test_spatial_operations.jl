@testset "test set_data_by_vector and vec" begin 
    # TODO. more test cases
    V = SU2Space(1//2=>1, 0=>1);
    P = SU2Space(1//2=>1);
    T = zeros(ComplexF64, P, V^4);

    mt = mapping_table(T)
    num_paras = num_free_parameters(T; _mapping_table=mt)

    v = rand(num_paras)

    T1 = set_data_by_vector(T, v; _mapping_table=mt)
    @test norm(vec(T1) - v) < 1e-12

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
    V = SU2Space(1//2=>1, 0=>1)
    P = SU2Space(1//2=>1)
    T = zeros(ComplexF64, P, V^4)

    @test_throws ArgumentError find_subspace(C4v(), T, :BadRep)
    @test_throws ArgumentError find_subspace(C6v(), T, :BadRep)
    @test_throws ArgumentError find_subspace(C4(), T, :BadRep)
    @test_throws ArgumentError find_subspace(D2(), T, :BadRep)
end

@testset "pointgroup verbosity" begin
    V = SU2Space(1//2=>1, 0=>1)
    P = SU2Space(1//2=>1)
    T = zeros(ComplexF64, P, V^4)

    P_sol = find_subspace(C4v(), T, :A1; verbose=false)
    @test size(P_sol, 1) == num_free_parameters(T)
end

@testset "non-orthonormal tolerance" begin
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
    V = SU2Space(1//2=>1, 0=>1)
    P = SU2Space(1//2=>1)
    T = zeros(ComplexF64, P, V^4)

    mt = mapping_table(T)
    num_paras = num_free_parameters(T; _mapping_table=mt)
    Psel = selector(T, (f1, f2) -> true; _mapping_table=mt)
    @test Psel == Matrix{ComplexF64}(I, num_paras, num_paras)
end

@testset "selector columns" begin
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
