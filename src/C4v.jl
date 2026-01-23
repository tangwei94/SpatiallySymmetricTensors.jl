"""C4v point group (square lattice with reflections and rotations)."""
struct C4v <: AbstractPointGroup end

const C4v_perm_type = Tuple{Tuple{Int}, NTuple{4, Int}}
const C4v_ops = Dict{Symbol, C4v_perm_type}(
    :Id => ((1, ), (2, 3, 4, 5)),
    :σd1 => ((1, ), (3, 2, 5, 4)),
    :σd2 => ((1, ), (5, 4, 3, 2)),
    :σv1 => ((1, ), (4, 3, 2, 5)),
    :σv2 => ((1, ), (2, 5, 4, 3)),
    :R1 => ((1, ), (3, 4, 5, 2)),
    :R3 => ((1, ), (5, 2, 3, 4)),
    :C2 => ((1, ), (4, 5, 2, 3)),
)
const C4v_A1_reps = Dict{Symbol, Int}(
    :Id => 1,
    :σd1 => 1, :σd2 => 1,
    :σv1 => 1, :σv2 => 1,
    :R1 => 1, :R3 => 1,
    :C2 => 1,
)
const C4v_A2_reps = Dict{Symbol, Int}(
    :Id => 1,
    :σd1 => -1, :σd2 => -1,
    :σv1 => -1, :σv2 => -1,
    :R1 => 1, :R3 => 1,
    :C2 => 1,
)
const C4v_B1_reps = Dict{Symbol, Int}(
    :Id => 1,
    :σd1 => -1, :σd2 => -1,
    :σv1 => 1, :σv2 => 1,
    :R1 => -1, :R3 => -1,
    :C2 => 1,
)
const C4v_B2_reps = Dict{Symbol, Int}(
    :Id => 1,
    :σd1 => 1, :σd2 => 1,
    :σv1 => -1, :σv2 => -1,
    :R1 => -1, :R3 => -1,
    :C2 => 1,
)
const C4v_E_reps = Dict{Symbol, Int}(
    :Id => 2,
    :σd1 => 0, :σd2 => 0,
    :σv1 => 0, :σv2 => 0,
    :R1 => 0, :R3 => 0,
    :C2 => -2,
)
const C4v_E_rep = Dict{Symbol, Matrix{ComplexF64}}(
    :Id => ComplexF64[1 0; 0 1],
    :R1 => ComplexF64[0 -1; 1 0],
    :R3 => ComplexF64[0 1; -1 0],
    :C2 => ComplexF64[-1 0; 0 -1],
    :σv1 => ComplexF64[1 0; 0 -1],
    :σv2 => ComplexF64[-1 0; 0 1],
    :σd1 => ComplexF64[0 1; 1 0],
    :σd2 => ComplexF64[0 -1; -1 0],
)

group_elements(::C4v) = C4v_ops
irrep_chars(::C4v, ::Val{:A1}) = C4v_A1_reps
irrep_chars(::C4v, ::Val{:A2}) = C4v_A2_reps
irrep_chars(::C4v, ::Val{:B1}) = C4v_B1_reps
irrep_chars(::C4v, ::Val{:B2}) = C4v_B2_reps
irrep_chars(::C4v, ::Val{:E}) = Dict(name => tr(mat) for (name, mat) in C4v_E_rep)

irrep_rep(::C4v, ::Val{:A1}) = Dict(name => [χ;;] for (name, χ) in C4v_A1_reps)
irrep_rep(::C4v, ::Val{:A2}) = Dict(name => [χ;;] for (name, χ) in C4v_A2_reps)
irrep_rep(::C4v, ::Val{:B1}) = Dict(name => [χ;;] for (name, χ) in C4v_B1_reps)
irrep_rep(::C4v, ::Val{:B2}) = Dict(name => [χ;;] for (name, χ) in C4v_B2_reps)
irrep_rep(::C4v, ::Val{:E}) = C4v_E_rep

"""
    split_multiplets(::C4v, T::AbstractTensorMap, ::Val{:E}, P_sol::Matrix{<:Number}; _mapping_table=mapping_table(T))

Split the E-irrep subspace `P_sol` into two multiplet components using the
eigenspaces of `-im * R1` and the C4v E-representation basis.
Returns a pair `(P_sol_a, P_sol_b)` with `P_sol_b = R1 * P_sol_a`.
"""
function split_multiplets(::C4v, T::AbstractTensorMap, ::Val{:E}, P_sol::Matrix{<:Number}; _mapping_table::MappingTable=mapping_table(T))

    mt = _mapping_table

    f_R = linear_function_for_spatial_operation(C4v_ops[:R1])
    mat_R = matrix_for_linear_function(T, f_R; _mapping_table=mt)

    rep_R = irrep_rep(C4v(), :E)[:R1]
    Λ, Q = eigen( Hermitian(-im * rep_R) )

    @show Λ
    P_1 = find_subspace(T, P_sol, (x -> -im * f_R(x)); λ= Λ[1], _mapping_table=mt)
    P_2 = find_subspace(T, P_sol, (x -> -im * f_R(x)); λ= Λ[2], _mapping_table=mt)

    if size(P_1, 2) != size(P_2, 2)
        throw(ArgumentError("split_multiplets expects equal multiplicities for ± eigen-subspaces, got $(size(P_1, 2)) and $(size(P_2, 2))"))
    end

    P_sol_a = Matrix(qr(P_1 * Q[1, 1] + P_2 * Q[2, 1]).Q)

    f_σv1 = linear_function_for_spatial_operation(C4v_ops[:σv1])
    mat_σv1 = matrix_for_linear_function(T, f_σv1; _mapping_table=mt)

    rep_σv1 = irrep_rep(C4v(), :E)[:σv1]
    Λ, Q = eigen( Hermitian(rep_σv1) )
    @show Λ
    P_1 = find_subspace(T, P_sol_a, (x -> f_σv1(x)); λ= Λ[1], _mapping_table=mt)
    P_2 = find_subspace(T, P_sol_a, (x -> f_σv1(x)); λ= Λ[2], _mapping_table=mt)

    @show size(P_sol_a)
    @show size(P_1), size(P_2)

    P_sol_b = mat_R * P_sol_a

    return P_sol_a, P_sol_b
end
