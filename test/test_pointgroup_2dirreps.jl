function _check_projector_function_2dirreps(group, irrep_name, num_virtual_space::Int)
    V = SU2Space(1//2=>1, 0=>1)
    P = SU2Space(1//2=>1)
    T = rand(ComplexF64, P, V^num_virtual_space)

    fproj = SpatiallySymmetricTensors.projector_function(group, irrep_name)
    Tp = fproj(T)
    @test norm(fproj(Tp) - Tp) < 1e-10
end

@testset "projector_function point group E irreps" begin
    _check_projector_function_2dirreps(C4v(), :E, 4)
    _check_projector_function_2dirreps(C3v(), :E, 3)
    _check_projector_function_2dirreps(C6v(), :E1, 6)
    _check_projector_function_2dirreps(C6v(), :E2, 6)
end

function _check_find_solution_2dirreps(group, irrep_name, num_virtual_space::Int)
    V = SU2Space(1//2=>1, 0=>1)
    P = SU2Space(1//2=>1)
    T = rand(ComplexF64, P, V^num_virtual_space)

    T_sol = find_solution(group, T, irrep_name)
    @test length(T_sol) % 2 == 0
    nsol = length(T_sol) ÷ 2

    ops = SpatiallySymmetricTensors.group_elements(group)
    rep = SpatiallySymmetricTensors.irrep_rep(group, irrep_name)

    for ix in 1:nsol
        T1, T2 = T_sol[ix], T_sol[ix+nsol]
        for (gname, gperm) in ops
            rep_g = rep[gname]
            # [g(T1), g(T2)] =[T1 T2] * rep_g
            @test norm(T1 * rep_g[1, 1] + T2 * rep_g[2, 1] - permute(T1, gperm)) < 1e-10
            @test norm(T1 * rep_g[1, 2] + T2 * rep_g[2, 2] - permute(T2, gperm)) < 1e-10
        end
    end
end

@testset "find_solution point group  E irreps" begin
    _check_find_solution_2dirreps(C4v(), :E, 4)
    _check_find_solution_2dirreps(C3v(), :E, 3)
    _check_find_solution_2dirreps(C6v(), :E1, 6)
    _check_find_solution_2dirreps(C6v(), :E2, 6)
end