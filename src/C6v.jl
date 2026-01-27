"""C6v point group (hexagonal lattice with reflections and rotations)."""
struct C6v <: AbstractPointGroup end

const C6v_perm_type = Tuple{Tuple{Int}, NTuple{6, Int}}
const C6v_ops = Dict{Symbol, C6v_perm_type}(
    :Id => ((1, ), (2, 3, 4, 5, 6, 7)),
    :R1 => ((1, ), (3, 4, 5, 6, 7, 2)),
    :R2 => ((1, ), (4, 5, 6, 7, 2, 3)),
    :C2 => ((1, ), (5, 6, 7, 2, 3, 4)),
    :R4 => ((1, ), (6, 7, 2, 3, 4, 5)),
    :R5 => ((1, ), (7, 2, 3, 4, 5, 6)),
    :σd1 => ((1, ), (7, 6, 5, 4, 3, 2)),
    :σd2 => ((1, ), (3, 2, 7, 6, 5, 4)),
    :σd3 => ((1, ), (5, 4, 3, 2, 7, 6)),
    :σv1 => ((1, ), (6, 5, 4, 3, 2, 7)),
    :σv2 => ((1, ), (2, 7, 6, 5, 4, 3)),
    :σv3 => ((1, ), (4, 3, 2, 7, 6, 5)),
)

const C6v_A1_chars = Dict{Symbol, Int}(
    :Id => 1,
    :σd1 => 1, :σd2 => 1, :σd3 => 1,
    :σv1 => 1, :σv2 => 1, :σv3 => 1,
    :R1 => 1, :R2 => 1, :C2 => 1, :R4 => 1, :R5 => 1,
)
const C6v_A2_chars = Dict{Symbol, Int}(
    :Id => 1,
    :σd1 => -1, :σd2 => -1, :σd3 => -1,
    :σv1 => -1, :σv2 => -1, :σv3 => -1,
    :R1 => 1, :R2 => 1, :C2 => 1, :R4 => 1, :R5 => 1,
)
const C6v_B1_chars = Dict{Symbol, Int}(
    :Id => 1,
    :σd1 => -1, :σd2 => -1, :σd3 => -1,
    :σv1 => 1, :σv2 => 1, :σv3 => 1,
    :R1 => -1, :R2 => 1, :C2 => -1, :R4 => 1, :R5 => -1,
)
const C6v_B2_chars = Dict{Symbol, Int}(
    :Id => 1,
    :σd1 => 1, :σd2 => 1, :σd3 => 1,
    :σv1 => -1, :σv2 => -1, :σv3 => -1,
    :R1 => -1, :R2 => 1, :C2 => -1, :R4 => 1, :R5 => -1,
)
const C6v_E1_rep = Dict{Symbol, Matrix{ComplexF64}}(
    :Id => ComplexF64[1 0; 0 1],
    :R1 => ComplexF64[1/2 sqrt(3)/2; -sqrt(3)/2 1/2],
    :R2 => ComplexF64[-1/2 sqrt(3)/2; -sqrt(3)/2 -1/2],
    :C2 => ComplexF64[-1 0; 0 -1],
    :R4 => ComplexF64[-1/2 -sqrt(3)/2; sqrt(3)/2 -1/2],
    :R5 => ComplexF64[1/2 -sqrt(3)/2; sqrt(3)/2 1/2],
    :σv1 => ComplexF64[1 0; 0 -1],
    :σv2 => ComplexF64[-1/2 sqrt(3)/2; sqrt(3)/2 1/2],
    :σv3 => ComplexF64[-1/2 -sqrt(3)/2; -sqrt(3)/2 1/2],
    :σd1 => ComplexF64[1/2 sqrt(3)/2; sqrt(3)/2 -1/2],
    :σd2 => ComplexF64[-1 0; 0 1],
    :σd3 => ComplexF64[1/2 -sqrt(3)/2; -sqrt(3)/2 -1/2],
)
const C6v_E2_rep = Dict{Symbol, Matrix{ComplexF64}}(
    :Id => ComplexF64[1 0; 0 1],
    :R1 => ComplexF64[-1/2 sqrt(3)/2; -sqrt(3)/2 -1/2],
    :R2 => ComplexF64[-1/2 -sqrt(3)/2; sqrt(3)/2 -1/2],
    :C2 => ComplexF64[1 0; 0 1],
    :R4 => ComplexF64[-1/2 sqrt(3)/2; -sqrt(3)/2 -1/2],
    :R5 => ComplexF64[-1/2 -sqrt(3)/2; sqrt(3)/2 -1/2],
    :σv1 => ComplexF64[1 0; 0 -1],
    :σv2 => ComplexF64[-1/2 -sqrt(3)/2; -sqrt(3)/2 1/2],
    :σv3 => ComplexF64[-1/2 sqrt(3)/2; sqrt(3)/2 1/2],
    :σd1 => ComplexF64[-1/2 sqrt(3)/2; sqrt(3)/2 1/2],
    :σd2 => ComplexF64[1 0; 0 -1],
    :σd3 => ComplexF64[-1/2 -sqrt(3)/2; -sqrt(3)/2 1/2],
)

group_elements(::C6v) = C6v_ops
irrep_chars(::C6v, ::Val{:A1}) = C6v_A1_chars
irrep_chars(::C6v, ::Val{:A2}) = C6v_A2_chars
irrep_chars(::C6v, ::Val{:B1}) = C6v_B1_chars
irrep_chars(::C6v, ::Val{:B2}) = C6v_B2_chars
irrep_chars(::C6v, ::Val{:E1}) = Dict(name => tr(mat) for (name, mat) in C6v_E1_rep)
irrep_chars(::C6v, ::Val{:E2}) = Dict(name => tr(mat) for (name, mat) in C6v_E2_rep)
irrep_rep(::C6v, ::Val{:A1}) = Dict(name => [χ;;] for (name, χ) in C6v_A1_chars)
irrep_rep(::C6v, ::Val{:A2}) = Dict(name => [χ;;] for (name, χ) in C6v_A2_chars)
irrep_rep(::C6v, ::Val{:B1}) = Dict(name => [χ;;] for (name, χ) in C6v_B1_chars)
irrep_rep(::C6v, ::Val{:B2}) = Dict(name => [χ;;] for (name, χ) in C6v_B2_chars)
irrep_rep(::C6v, ::Val{:E1}) = C6v_E1_rep
irrep_rep(::C6v, ::Val{:E2}) = C6v_E2_rep

"""
    split_multiplets(::C6v, T::AbstractTensorMap, ::Val{:E1}, P_sol::Matrix{<:Number}; _mapping_table=mapping_table(T))
    split_multiplets(::C6v, T::AbstractTensorMap, ::Val{:E2}, P_sol::Matrix{<:Number}; _mapping_table=mapping_table(T))

Split the E1/E2-irrep subspace `P_sol` into two multiplet components.
"""
function split_multiplets(::C6v, T::AbstractTensorMap, rep::Union{Val{:E1}, Val{:E2}}, P_sol::Matrix{<:Number}; _mapping_table::MappingTable=mapping_table(T))
    mt = _mapping_table

    f_σv1 = linear_function_for_spatial_operation(C6v_ops[:σv1])
    rep_σv1 = irrep_rep(C6v(), rep)[:σv1] # rep_σv1 is diagonal, σv1 is hermitian
    P_1 = find_subspace(T, P_sol, f_σv1; λ= real(rep_σv1[1,1]), is_hermitian=true, _mapping_table=mt)
    P_2_rotated = find_subspace(T, P_sol, f_σv1; λ= real(rep_σv1[2,2]), is_hermitian=true, _mapping_table=mt)

    # `P_2_rotated` is not uniquely determined: it is defined only up to a gauge transformation `Q` within the irrep subspace, i.e. `P_2_rotated = P_2 * Q`.
    # We fix this gauge so that, in the chosen basis, the representation matches the E-irrep matrices in `C6v_E1_rep`/`C6v_E2_rep`.
    f_σd1 = linear_function_for_spatial_operation(C6v_ops[:σd1])
    rep_σd1 = irrep_rep(C6v(), rep)[:σd1] # rep_σd1 not diagonal
    mat_σd1 = matrix_for_linear_function(T, f_σd1; _mapping_table=mt)
    Qdag = inv(rep_σd1[2, 1]) * P_2_rotated' * mat_σd1 * P_1
    P_2 = P_2_rotated * Qdag

    return P_1, P_2
end
