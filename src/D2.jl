"""D2 point group (two reflections)."""
struct D2 <: AbstractPointGroup end

# in the order of σd, σv, R; see http://symmetry.jacobs-university.de/cgi-bin/group.cgi?group=302&option=4
# ignore representation of σd for D2
const D2_A_reps = (0, 1, 1) 
const D2_B1_reps = (0, -1, 1)

const D2_σv1 = ((1, ), (4, 3, 2, 5))
const D2_σv2 = ((1, ), (2, 5, 4, 3))

"""
    get_reps(::D2, name::Symbol)

Return the representation data for `name` in the `D2` point group.

Notes:
- Supported names are `:A` and `:B1`.
- The return value is a tuple `(σd, σv, R)` where `σd` and `R` are unused
  (set to 0 by convention) and `σv` are reflection eigenvalues.
"""
function get_reps(::D2, name::Symbol)
    (name == :A) && return D2_A_reps
    (name == :B1) && return D2_B1_reps
    throw(ArgumentError("unknown representation name $(name) for D2"))
end
"""
    get_perm(::D2, name::Symbol)

Return the list of spatial permutations for the operation family `name`.

Notes:
- Supported names are `:σd`, `:σv`, and `:R`.
- `:σd` and `:R` are empty for `D2` by convention.
"""
function get_perm(::D2, name::Symbol)
    (name == :σd) && return [] # empty diagonal reflections
    (name == :σv) && return [D2_σv1, D2_σv2]
    (name == :R) && return [] # empty rotations
    throw(ArgumentError("unknown operation name $(name) for D2"))
end
