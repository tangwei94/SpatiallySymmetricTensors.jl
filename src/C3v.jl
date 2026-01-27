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
