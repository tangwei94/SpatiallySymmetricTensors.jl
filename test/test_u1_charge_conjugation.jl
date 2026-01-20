@testset "test u1_charge_conjugation" begin
    # physical and virtual spaces
    P = U1Space(0=>1, 1=>1, 2=>1);
    Pa = U1Space(-1=>1);
    V = U1Space(-1=>1, 0=>1, 1=>1);

    # A1 PEPS tensors
    T = zeros(ComplexF64, fuse(P, Pa), V^4);

    # zeroth order, only one solution: (0, 0, 0, 0) -> (1, -1) 
    P0 = selector(T, 
        (f1, f2) -> (length(findall(rep -> rep == U1Irrep(0), f2.uncoupled)) == 4)
        )
    T0 = find_solution(C4v(), T, :A1; P_filter=P0)[1];

    # first order, two solutions: (0, 0, 0, 1) -> (2, -1) and (0, 0, 0, -1) -> (0, -1)
    P1 = selector(T, 
        (f1, f2) -> (length(findall(rep -> rep == U1Irrep(0), f2.uncoupled)) == 3)
        )
    T1s = find_solution(C4v(), T, :A1; P_filter=P1);
    # make sure the phase of both first-order solutions are fixed
    function fix_phase(v) 
        for ix in eachindex(v.data)
            if norm(v.data[ix]) > 1e-12
                α = v.data[ix]
                return v / α
            end
        end
        return v
    end
    T1_even = (fix_phase(T1s[1]) + fix_phase(T1s[2])) 
    T1_odd = (fix_phase(T1s[1]) - fix_phase(T1s[2])) 

    @test norm(T0 - u1_charge_conjugation(T0)) < 1e-12
    @test norm(T1_even - u1_charge_conjugation(T1_even)) < 1e-12
    @test norm(T1_odd + u1_charge_conjugation(T1_odd)) < 1e-12

    P1_A1 = find_subspace(C4v(), T, :A1; P_filter=P1)
    P1_A1_conj_even = find_subspace_for_u1_charge_conjugation(T, P1_A1)
    T1_even_auto = find_solution(T, P1_A1_conj_even)[1]

    @test norm(T1_even_auto - T1_even / norm(T1_even)) < 1e-12
end

@testset "u1_charge_conjugation missing blocks" begin
    P = U1Space(1=>1)
    V = U1Space(0=>1, 1=>1)
    T = zeros(ComplexF64, P, V^2)

    @test_throws ArgumentError u1_charge_conjugation(T)
end
