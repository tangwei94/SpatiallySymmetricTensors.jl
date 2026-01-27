"""C3v point group (triangular lattice with reflections and 3-fold rotations)."""
struct C3v <: AbstractPointGroup end

const C3v_perm_type = Tuple{Tuple{Int}, NTuple{3, Int}}
const C3v_ops = Dict{Symbol, C3v_perm_type}(
    :Id => ((1, ), (2, 3, 4)),
    :R1 => ((1, ), (3, 4, 2)),
    :R2 => ((1, ), (4, 2, 3)),
    :σv1 => ((1, ), (3, 2, 4)),
    :σv2 => ((1, ), (2, 4, 3)),
    :σv3 => ((1, ), (4, 3, 2)),
)

const C3v_A1_reps = Dict{Symbol, Int}(
    :Id => 1,
    :R1 => 1, :R2 => 1,
    :σv1 => 1, :σv2 => 1, :σv3 => 1,
)
const C3v_A2_reps = Dict{Symbol, Int}(
    :Id => 1,
    :R1 => 1, :R2 => 1,
    :σv1 => -1, :σv2 => -1, :σv3 => -1,
)
const C3v_E_reps = Dict{Symbol, Int}(
    :Id => 2,
    :R1 => -1, :R2 => -1,
    :σv1 => 0, :σv2 => 0, :σv3 => 0,
)
const C3v_E_rep = Dict{Symbol, Matrix{ComplexF64}}(
    :Id => ComplexF64[1 0; 0 1],
    :R1 => ComplexF64[-1/2 -sqrt(3)/2; sqrt(3)/2 -1/2],
    :R2 => ComplexF64[-1/2 sqrt(3)/2; -sqrt(3)/2 -1/2],
    :σv1 => ComplexF64[1 0; 0 -1],
    :σv2 => ComplexF64[-1/2 -sqrt(3)/2; -sqrt(3)/2 1/2],
    :σv3 => ComplexF64[-1/2 sqrt(3)/2; sqrt(3)/2 1/2],
)

group_elements(::C3v) = C3v_ops
irrep_chars(::C3v, ::Val{:A1}) = C3v_A1_reps
irrep_chars(::C3v, ::Val{:A2}) = C3v_A2_reps
irrep_chars(::C3v, ::Val{:E}) = Dict(name => tr(mat) for (name, mat) in C3v_E_rep)
irrep_rep(::C3v, ::Val{:A1}) = Dict(name => [χ;;] for (name, χ) in C3v_A1_reps)
irrep_rep(::C3v, ::Val{:A2}) = Dict(name => [χ;;] for (name, χ) in C3v_A2_reps)
irrep_rep(::C3v, ::Val{:E}) = C3v_E_rep

"""
    split_multiplets(::C3v, T::AbstractTensorMap, ::Val{:E}, P_sol::Matrix{<:Number}; _mapping_table=mapping_table(T))

Split the E-irrep subspace `P_sol` into two multiplet components.
"""
function split_multiplets(::C3v, T::AbstractTensorMap, ::Val{:E}, P_sol::Matrix{<:Number}; _mapping_table::MappingTable=mapping_table(T))

    mt = _mapping_table

    f_σv1 = linear_function_for_spatial_operation(C3v_ops[:σv1])
    rep_σv1 = irrep_rep(C3v(), :E)[:σv1] # rep_σv1 is diagonal
    P_1 = find_subspace(T, P_sol, f_σv1; λ= real(rep_σv1[1,1]), _mapping_table=mt)

    f_σv2 = linear_function_for_spatial_operation(C3v_ops[:σv2])
    mat_σv2 = matrix_for_linear_function(T, f_σv2; _mapping_table=mt)
    P_2 = mat_σv2 * P_1 # fix the basis of P_2 to be consistent with P_1

    return P_1, P_2
end
