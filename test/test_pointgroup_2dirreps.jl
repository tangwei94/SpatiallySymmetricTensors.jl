@testset "projector_function C4v E irrep" begin
    V = SU2Space(1//2=>1, 0=>1)
    P = SU2Space(1//2=>1)
    T = rand(ComplexF64, P, V^4)

    fproj = SpatiallySymmetricTensors.projector_function(C4v(), :E)
    Tp = fproj(T)
    @test norm(fproj(Tp) - Tp) < 1e-10
end

@testset "find_solution C4v E irrep" begin

    V = SU2Space(1//2=>1, 0=>1)
    P = SU2Space(1//2=>1)
    T = rand(ComplexF64, P, V^4)

    T_sol = find_solution(C4v(), T, :E)
    @test length(T_sol) % 2 == 0
    nsol = length(T_sol) ÷ 2

    C4v_ops = SpatiallySymmetricTensors.group_elements(C4v())
    C4v_E_rep = SpatiallySymmetricTensors.irrep_rep(C4v(), :E)

    #T1, T2 = T_sol[1], T_sol[nsol+1]
    #for (gname, gperm) in C4v_ops
    #    @show gname
    #    rep_g = C4v_E_rep[gname]
    #    @show norm(rep_g[1, 1] * T1 + rep_g[2, 1] * T2 - permute(T1, gperm))
    #    #@test norm(rep_g[1, 1] * T1 + rep_g[2, 1] * T2 - permute(T1, gperm)) < 1e-10
    #end
end