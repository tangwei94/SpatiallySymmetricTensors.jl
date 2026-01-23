"""D2 point group (two reflections and a 180 degree rotation)."""
struct D2 <: AbstractPointGroup end

const D2_perm_type = Tuple{Tuple{Int}, NTuple{4, Int}}
const D2_ops = Dict{Symbol, D2_perm_type}(
    :Id => ((1, ), (2, 3, 4, 5)),
    :σv1 => ((1, ), (4, 3, 2, 5)),
    :σv2 => ((1, ), (2, 5, 4, 3)),
    :C2 => ((1, ), (4, 5, 2, 3)),
)
const D2_A_reps = Dict{Symbol, Int}(
    :Id => 1, :σv1 => 1, :σv2 => 1, :C2 => 1,
)
const D2_B1_reps = Dict{Symbol, Int}(
    :Id => 1, :σv1 => -1, :σv2 => -1, :C2 => 1,
)
const D2_B2_reps = Dict{Symbol, Int}(
    :Id => 1, :σv1 => 1, :σv2 => -1, :C2 => -1,
)
const D2_B3_reps = Dict{Symbol, Int}(
    :Id => 1, :σv1 => -1, :σv2 => 1, :C2 => -1,
)

group_elements(::D2) = D2_ops
irrep_chars(::D2, ::Val{:A}) = D2_A_reps
irrep_chars(::D2, ::Val{:B1}) = D2_B1_reps
irrep_chars(::D2, ::Val{:B2}) = D2_B2_reps
irrep_chars(::D2, ::Val{:B3}) = D2_B3_reps
