"""C4v point group (square lattice with reflections and rotations)."""
struct C4v <: AbstractPointGroup end

const C4v_σd1 = ((1, ), (3, 2, 5, 4))
const C4v_σd2 = ((1, ), (5, 4, 3, 2))
const C4v_σv1 = ((1, ), (4, 3, 2, 5))
const C4v_σv2 = ((1, ), (2, 5, 4, 3))
const C4v_R1 = ((1, ), (3, 4, 5, 2))
const C4v_R2 = ((1, ), (5, 2, 3, 4))

# in the order of σd, σv, R; see http://symmetry.jacobs-university.de/cgi-bin/group.cgi?group=404&option=4
const C4v_A1_reps = (1, 1, 1) 
const C4v_A2_reps = (-1, -1, 1) 
const C4v_B1_reps = (-1, 1, 1) 
const C4v_B2_reps = (1, -1, 1) 

# E irrep in the (x, y) basis
const C4v_E_σd1 = [0 1; 1 0]
const C4v_E_σd2 = [0 -1; -1 0]
const C4v_E_σv1 = [1 0; 0 -1]
const C4v_E_σv2 = [-1 0; 0 1]
const C4v_E_R1 = [0 -1; 1 0]
const C4v_E_R2 = [0 1; -1 0]
const C4v_E_reps = ([C4v_E_σd1, C4v_E_σd2], [C4v_E_σv1, C4v_E_σv2], [C4v_E_R1, C4v_E_R2])

"""
    get_reps(::C4v, name::Symbol)

Return the representation data for `name` in the `C4v` point group.

Notes:
- Supported names are `:A1`, `:A2`, `:B1`, `:B2`, and `:E`.
- The return value is a tuple `(σd, σv, R)` where each entry is either a scalar
  eigenvalue or a list of matrices for the 2D irrep `:E`.
"""
function get_reps(::C4v, name::Symbol)
    (name == :A1) && return C4v_A1_reps
    (name == :A2) && return C4v_A2_reps
    (name == :B1) && return C4v_B1_reps
    (name == :B2) && return C4v_B2_reps
    (name == :E) && return C4v_E_reps
    throw(ArgumentError("unknown representation name $(name) for C4v"))
end
"""
    get_perm(::C4v, name::Symbol)

Return the list of spatial permutations for the operation family `name`.

Notes:
- Supported names are `:σd`, `:σv`, and `:R`.
- The permutations follow the leg ordering convention used by `permute`.
"""
function get_perm(::C4v, name::Symbol)
    (name == :σd) && return [C4v_σd1, C4v_σd2] 
    (name == :σv) && return [C4v_σv1, C4v_σv2]
    (name == :R) && return [C4v_R1, C4v_R2]
    throw(ArgumentError("unknown operation name $(name) for C4v"))
end
