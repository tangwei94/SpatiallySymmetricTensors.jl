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

    Split the E-irrep subspace `P_sol` into two multiplet components.
"""
function split_multiplets(::C4v, T::AbstractTensorMap, ::Val{:E}, P_sol::Matrix{<:Number}; _mapping_table::MappingTable=mapping_table(T))

    mt = _mapping_table
    
    f_σv1 = linear_function_for_spatial_operation(C4v_ops[:σv1])
    rep_σv1 = irrep_rep(C4v(), :E)[:σv1] # rep_σv1 is diagonal
    P_1 = find_subspace(T, P_sol, f_σv1; λ= real(rep_σv1[1,1]), _mapping_table=mt)
    P_2_rotated = find_subspace(T, P_sol, f_σv1; λ= real(rep_σv1[2,2]), _mapping_table=mt)

    # `P_2_rotated` is not uniquely determined: it is defined only up to a gauge transformation `Q` within the irrep subspace, i.e. `P_2_rotated = P_2 * Q`.
    # We fix this gauge so that, in the chosen basis, the representation matches the E-irrep matrices in `C4v_E_rep`.
    # In practice this may be unnecessary, since the matrices in `C4v_E_rep` are not unique themselves and can be rotated by unitary transformations.
    f_σd1 = linear_function_for_spatial_operation(C4v_ops[:σd1])
    rep_σd1 = irrep_rep(C4v(), :E)[:σd1] # rep_σd1 not diagonal
    mat_σd1 = matrix_for_linear_function(T, f_σd1; _mapping_table=mt)
    Qdag = inv(rep_σd1[2, 1]) * P_2_rotated' * mat_σd1 * P_1 
    P_2 = P_2_rotated * Qdag

    return P_1, P_2
end
