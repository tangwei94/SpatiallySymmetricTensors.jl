"""C3v point group (triangular lattice with reflections and 3-fold rotations)."""
struct C3v <: AbstractPointGroup end

const C3v_R1 = ((1, ), (3, 4, 2))
const C3v_R2 = ((1, ), (4, 2, 3))
const C3v_σv1 = ((1, ), (3, 2, 4))
const C3v_σv2 = ((1, ), (2, 4, 3))
const C3v_σv3 = ((1, ), (4, 3, 2))

# in the order of σd, σv, R; for C3v use σv for reflections and ignore σd
const C3v_A1_reps = (0, 1, 1)
const C3v_A2_reps = (0, -1, 1)

function get_reps(::C3v, name::Symbol)
    (name == :A1) && return C3v_A1_reps
    (name == :A2) && return C3v_A2_reps
    throw(ArgumentError("unknown representation name $(name) for C3v"))
end
function get_perm(::C3v, name::Symbol)
    (name == :σd) && return [] # empty diagonal reflections
    (name == :σv) && return [C3v_σv1, C3v_σv2, C3v_σv3]
    (name == :R) && return [C3v_R1, C3v_R2]
    throw(ArgumentError("unknown operation name $(name) for C3v"))
end
