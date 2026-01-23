"""C4 point group (rotations only)."""
struct C4 <: AbstractPointGroup end

const C4_perm_type = Tuple{Tuple{Int}, NTuple{4, Int}}
const C4_ops = Dict{Symbol, C4_perm_type}(
    :Id => ((1, ), (2, 3, 4, 5)),
    :R1 => ((1, ), (3, 4, 5, 2)),
    :R3 => ((1, ), (5, 2, 3, 4)),
    :C2 => ((1, ), (4, 5, 2, 3)),
)
const C4_A_reps = Dict{Symbol, Int}(
    :Id => 1, :R1 => 1, :R3 => 1, :C2 => 1,
)
const C4_B_reps = Dict{Symbol, Int}(
    :Id => 1, :R1 => -1, :R3 => -1, :C2 => 1,
)

group_elements(::C4) = C4_ops
irrep_chars(::C4, ::Val{:A}) = C4_A_reps
irrep_chars(::C4, ::Val{:B}) = C4_B_reps
