"""C4 point group (rotations only)."""
struct C4 <: AbstractPointGroup end

# in the order of σd, σv, R; see http://symmetry.jacobs-university.de/cgi-bin/group.cgi?group=204&option=4
# ignore representationq of σd and σv for C4
const C4_A_reps = (0, 0, 1) 
const C4_B_reps = (0, 0, -1) 

const C4_R1 = ((1, ), (3, 4, 5, 2))
const C4_R2 = ((1, ), (5, 2, 3, 4))

"""
    get_reps(::C4, name::Symbol)

Return the representation data for `name` in the `C4` point group.

Notes:
- Supported names are `:A` and `:B`.
- The return value is a tuple `(σd, σv, R)` where `σd` and `σv` are unused
  (set to 0 by convention) and `R` is the rotation eigenvalue.
"""
function get_reps(::C4, name::Symbol)
    (name == :A) && return C4_A_reps
    (name == :B) && return C4_B_reps
    throw(ArgumentError("unknown representation name $(name) for C4"))
end
"""
    get_perm(::C4, name::Symbol)

Return the list of spatial permutations for the operation family `name`.

Notes:
- Supported names are `:σd`, `:σv`, and `:R`.
- `:σd` and `:σv` are empty for `C4` by convention.
"""
function get_perm(::C4, name::Symbol)
    (name == :σd) && return [] # empty diagonal reflections
    (name == :σv) && return [] # empty horizontal and vertical reflections
    (name == :R) && return [C4_R1, C4_R2]
    throw(ArgumentError("unknown operation name $(name) for C4"))
end
