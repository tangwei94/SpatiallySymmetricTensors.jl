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

group_elements(::C3v) = C3v_ops
irrep_chars(::C3v, ::Val{:A1}) = C3v_A1_reps
irrep_chars(::C3v, ::Val{:A2}) = C3v_A2_reps
