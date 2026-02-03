function T_1_3_A1()
    V = SU2Space(1//2=>1, 0=>1)
    P = SU2Space(1//2=>1)
    T = zeros(ComplexF64, P, V^4)

    # a projector to the subspace with (1//2 ⊕ 0 ⊕ 0 ⊕ 0 -> 1//2) (short-range RVB)
    P_nocc_1_3 = begin
        _condition(f1, f2) = length(findall(rep-> rep == SU2Irrep(1//2), f2.uncoupled)) == 1
        selector(T, _condition)
    end

    T_1_3_A1 = find_solution(C4v(), T, :A1; P_filter=P_nocc_1_3)[1]
    return T_1_3_A1
end

function T_1_3_A1_from_plain()
    A = zeros(2, 3, 3, 3, 3)
    A[1, 2, 1, 1, 1] = A[1, 1, 2, 1, 1] = A[1, 1, 1, 2, 1] = A[1, 1, 1, 1, 2] = 1/2
    A[2, 3, 1, 1, 1] = A[2, 1, 3, 1, 1] = A[2, 1, 1, 3, 1] = A[2, 1, 1, 1, 3] = 1/2

    V = SU2Space(0=>1, 1//2=>1) # 1 corresponds to 0=>1, 2~3 corresponds to 1//2=>1
    P = SU2Space(1//2=>1)
    T = TensorMap(A, P, V^4)
    return T
end

function T_3_1_A1()
    V = SU2Space(1//2=>1, 0=>1)
    P = SU2Space(1//2=>1)
    T = zeros(ComplexF64, P, V^4)

    # a projector to the subspace with (1//2 ⊕ 1//2 ⊕ 1//2 ⊕ 0 -> 1//2) (long-range RVB)
    P_nocc_3_1 = begin
        _condition(f1, f2) = length(findall(rep-> rep == SU2Irrep(1//2), f2.uncoupled)) == 3
        selector(T, _condition)
    end

    T_3_1_A1 = find_solution(C4v(), T, :A1; P_filter=P_nocc_3_1)[1]
    return T_3_1_A1
end

function T_3_1_A1_from_plain()
    A = zeros(2, 3, 3, 3, 3)
    A[1, 2, 2, 3, 1] = A[1, 2, 2, 1, 3] = A[1, 2, 3, 1, 2] = A[1, 2, 1, 3, 2] = A[1, 3, 2, 2, 1] = A[1, 3, 1, 2, 2] = A[1, 1, 2, 2, 3] = A[1, 1, 3, 2, 2] = -1/2/sqrt(6)
    A[1, 2, 3, 2, 1] = A[1, 2, 1, 2, 3] = A[1, 3, 2, 1, 2] = A[1, 1, 2, 3, 2] = 1/sqrt(6)

    A[2, 2, 3, 3, 1] = A[2, 2, 1, 3, 3] = A[2, 3, 2, 1, 3] = A[2, 3, 3, 2, 1] = A[2, 3, 3, 1, 2] = A[2, 3, 1, 2, 3]=  A[2, 1, 2, 3, 3] = A[2, 1, 3, 3, 2] = 1/2/sqrt(6)
    A[2, 2, 3, 1, 3] = A[2, 3, 2, 3, 1] = A[2, 3, 1, 3, 2] = A[2, 1, 3, 2, 3] = -1/sqrt(6)

    V = SU2Space(0=>1, 1//2=>1) # 1 corresponds to 0=>1, 2~3 corresponds to 1//2=>1
    P = SU2Space(1//2=>1)
    T = TensorMap(A, P, V^4)
    return T
end

@testset "test T_1_3_A1" begin
    @show "test T_1_3_A1, short-range RVB state"
    T1 = T_1_3_A1()
    T2 = T_1_3_A1_from_plain()

    @test norm(T1 - T2) < 1e-12
end

@testset "test T_3_1_A1" begin
    @show "test T_3_1_A1, long-range RVB state"
    T1 = T_3_1_A1()
    T2 = T_3_1_A1_from_plain()

    @test norm(T1 - T2) < 1e-12
end
