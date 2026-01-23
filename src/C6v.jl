"""C6v point group (hexagonal lattice with reflections and rotations)."""
struct C6v <: AbstractPointGroup end

const C6v_perm_type = Tuple{Tuple{Int}, NTuple{6, Int}}
const C6v_ops = Dict{Symbol, C6v_perm_type}(
    :Id => ((1, ), (2, 3, 4, 5, 6, 7)),
    :R1 => ((1, ), (3, 4, 5, 6, 7, 2)),
    :R2 => ((1, ), (7, 2, 3, 4, 5, 6)),
    :σd1 => ((1, ), (7, 6, 5, 4, 3, 2)),
    :σd2 => ((1, ), (3, 2, 7, 6, 5, 4)),
    :σd3 => ((1, ), (5, 4, 3, 2, 7, 6)),
    :σv1 => ((1, ), (6, 5, 4, 3, 2, 7)),
    :σv2 => ((1, ), (2, 7, 6, 5, 4, 3)),
    :σv3 => ((1, ), (4, 3, 2, 7, 6, 5)),
)

const C6v_A1_reps = Dict{Symbol, Int}(
    :Id => 1,
    :σd1 => 1, :σd2 => 1, :σd3 => 1,
    :σv1 => 1, :σv2 => 1, :σv3 => 1,
    :R1 => 1, :R2 => 1,
)
const C6v_A2_reps = Dict{Symbol, Int}(
    :Id => 1,
    :σd1 => -1, :σd2 => -1, :σd3 => -1,
    :σv1 => -1, :σv2 => -1, :σv3 => -1,
    :R1 => 1, :R2 => 1,
)
const C6v_B1_reps = Dict{Symbol, Int}(
    :Id => 1,
    :σd1 => -1, :σd2 => -1, :σd3 => -1,
    :σv1 => 1, :σv2 => 1, :σv3 => 1,
    :R1 => -1, :R2 => -1,
)
const C6v_B2_reps = Dict{Symbol, Int}(
    :Id => 1,
    :σd1 => 1, :σd2 => 1, :σd3 => 1,
    :σv1 => -1, :σv2 => -1, :σv3 => -1,
    :R1 => -1, :R2 => -1,
)
