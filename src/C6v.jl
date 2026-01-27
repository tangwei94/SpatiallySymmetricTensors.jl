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

const C6v_A1_reps = Dict{Symbol, Int}(
    :Id => 1,
    :σd1 => 1, :σd2 => 1, :σd3 => 1,
    :σv1 => 1, :σv2 => 1, :σv3 => 1,
    :R1 => 1, :R2 => 1, :C2 => 1, :R4 => 1, :R5 => 1,
)
const C6v_A2_reps = Dict{Symbol, Int}(
    :Id => 1,
    :σd1 => -1, :σd2 => -1, :σd3 => -1,
    :σv1 => -1, :σv2 => -1, :σv3 => -1,
    :R1 => 1, :R2 => 1, :C2 => 1, :R4 => 1, :R5 => 1,
)
const C6v_B1_reps = Dict{Symbol, Int}(
    :Id => 1,
    :σd1 => -1, :σd2 => -1, :σd3 => -1,
    :σv1 => 1, :σv2 => 1, :σv3 => 1,
    :R1 => -1, :R2 => 1, :C2 => -1, :R4 => 1, :R5 => -1,
)
const C6v_B2_reps = Dict{Symbol, Int}(
    :Id => 1,
    :σd1 => 1, :σd2 => 1, :σd3 => 1,
    :σv1 => -1, :σv2 => -1, :σv3 => -1,
    :R1 => -1, :R2 => 1, :C2 => -1, :R4 => 1, :R5 => -1,
)
const C6v_E1_rep = Dict{Symbol, Matrix{ComplexF64}}(
    :Id => ComplexF64[1 0; 0 1],
    :R1 => ComplexF64[1/2 -sqrt(3)/2; sqrt(3)/2 1/2],
    :R2 => ComplexF64[-1/2 -sqrt(3)/2; sqrt(3)/2 -1/2],
    :C2 => ComplexF64[-1 0; 0 -1],
    :R4 => ComplexF64[-1/2 sqrt(3)/2; -sqrt(3)/2 -1/2],
    :R5 => ComplexF64[1/2 sqrt(3)/2; -sqrt(3)/2 1/2],
    :σv1 => ComplexF64[1 0; 0 -1],
    :σv2 => ComplexF64[-1/2 sqrt(3)/2; sqrt(3)/2 1/2],
    :σv3 => ComplexF64[-1/2 -sqrt(3)/2; -sqrt(3)/2 1/2],
    :σd1 => ComplexF64[1/2 sqrt(3)/2; sqrt(3)/2 -1/2],
    :σd2 => ComplexF64[-1 0; 0 1],
    :σd3 => ComplexF64[1/2 -sqrt(3)/2; -sqrt(3)/2 -1/2],
)
const C6v_E2_rep = Dict{Symbol, Matrix{ComplexF64}}(
    :Id => ComplexF64[1 0; 0 1],
    :R1 => ComplexF64[-1/2 -sqrt(3)/2; sqrt(3)/2 -1/2],
    :R2 => ComplexF64[-1/2 sqrt(3)/2; -sqrt(3)/2 -1/2],
    :C2 => ComplexF64[1 0; 0 1],
    :R4 => ComplexF64[-1/2 -sqrt(3)/2; sqrt(3)/2 -1/2],
    :R5 => ComplexF64[-1/2 sqrt(3)/2; -sqrt(3)/2 -1/2],
    :σv1 => ComplexF64[1 0; 0 -1],
    :σv2 => ComplexF64[-1/2 -sqrt(3)/2; -sqrt(3)/2 1/2],
    :σv3 => ComplexF64[-1/2 sqrt(3)/2; sqrt(3)/2 1/2],
    :σd1 => ComplexF64[-1/2 sqrt(3)/2; sqrt(3)/2 1/2],
    :σd2 => ComplexF64[1 0; 0 -1],
    :σd3 => ComplexF64[-1/2 -sqrt(3)/2; -sqrt(3)/2 1/2],
)

group_elements(::C6v) = C6v_ops
irrep_chars(::C6v, ::Val{:A1}) = C6v_A1_reps
irrep_chars(::C6v, ::Val{:A2}) = C6v_A2_reps
irrep_chars(::C6v, ::Val{:B1}) = C6v_B1_reps
irrep_chars(::C6v, ::Val{:B2}) = C6v_B2_reps
irrep_chars(::C6v, ::Val{:E1}) = Dict(name => tr(mat) for (name, mat) in C6v_E1_rep)
irrep_chars(::C6v, ::Val{:E2}) = Dict(name => tr(mat) for (name, mat) in C6v_E2_rep)
irrep_rep(::C6v, ::Val{:A1}) = Dict(name => [χ;;] for (name, χ) in C6v_A1_reps)
irrep_rep(::C6v, ::Val{:A2}) = Dict(name => [χ;;] for (name, χ) in C6v_A2_reps)
irrep_rep(::C6v, ::Val{:B1}) = Dict(name => [χ;;] for (name, χ) in C6v_B1_reps)
irrep_rep(::C6v, ::Val{:B2}) = Dict(name => [χ;;] for (name, χ) in C6v_B2_reps)
irrep_rep(::C6v, ::Val{:E1}) = C6v_E1_rep
irrep_rep(::C6v, ::Val{:E2}) = C6v_E2_rep
