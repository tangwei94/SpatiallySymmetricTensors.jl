"""
Abstract point group type for spatial symmetry operations.
"""
abstract type AbstractPointGroup end

"""
    group_elements(spg::AbstractPointGroup)

Return a dict of group elements (as permutations) keyed by symbols.
"""
function group_elements(::AbstractPointGroup)
    # if the input does not match any of the defined point groups, throw an error
    throw(ArgumentError("group_elements not implemented for this point group"))
end

"""
    irrep_chars(spg::AbstractPointGroup, ::Val{reps_name})

Return a dict of characters for a 1D irrep keyed by element symbols.
"""
function irrep_chars(::AbstractPointGroup, ::Val{reps_name}) where {reps_name}
    # if the input does not match any of the defined point groups, throw an error
    throw(ArgumentError("irrep_chars not implemented for this point group"))
end

irrep_chars(spg::AbstractPointGroup, reps_name::Symbol) = irrep_chars(spg, Val(reps_name))

"""
    projector_function(spg::AbstractPointGroup, reps_name::Symbol)

Return a function that applies the irrep projector: P = d_reps / |G| ∑_{g ∈ G} conj(χ)(g) ρ(g)
"""
function projector_function(spg::AbstractPointGroup, reps_name::Symbol)
    ops = group_elements(spg)
    reps = irrep_chars(spg, reps_name)

    perms_with_char = Tuple{Any, ComplexF64}[]
    for (op_name, perm) in ops
        χ = get(reps, op_name, nothing)
        χ === nothing && throw(ArgumentError("missing rep value for $(op_name) in $(reps_name)"))
        push!(perms_with_char, (perm, ComplexF64(χ)))
    end
    isempty(perms_with_char) && throw(ArgumentError("no permutations found for $(spg)"))

    norm_factor = 1.0 / length(perms_with_char)

    return function (T)
        acc = zero(T)
        for (perm, χ) in perms_with_char
            acc += conj(χ) * permute(T, perm)
        end
        return norm_factor * acc
    end
end



"""
    find_subspace(spg::AbstractPointGroup, T::AbstractTensorMap, reps_name::Symbol; P_filter=nothing, verbose=false, tol=1e-8)

Find the subspace of `T` transforming under the representation `reps_name` of `spg`.
"""
function find_subspace(spg::AbstractPointGroup, T::AbstractTensorMap, reps_name::Symbol; P_filter=nothing, verbose::Bool=false, tol::Real=1e-8)
    reps = get_reps(spg, reps_name)
    mt = mapping_table(T)

    if isnothing(P_filter)
        num_paras = num_free_parameters(T; _mapping_table=mt)
        P_filter = Matrix{ComplexF64}(I, num_paras, num_paras)
    end

    P_sol = P_filter

    if verbose
        println("init, size(P_sol) = ", size(P_sol))
    end
    for permutation in get_perm(spg, :σd)
        f_op = linear_function_for_spatial_operation(permutation)
        P_sol = find_subspace(T, P_sol, f_op; λ=reps[1], is_hermitian=true, tol=tol, _mapping_table=mt)
        if verbose
            println("operation σd, size(P_sol) = ", size(P_sol))
        end
    end
    for permutation in get_perm(spg, :σv)
        f_op = linear_function_for_spatial_operation(permutation)
        P_sol = find_subspace(T, P_sol, f_op; λ=reps[2], is_hermitian=true, tol=tol, _mapping_table=mt)
        if verbose
            println("operation σv, size(P_sol) = ", size(P_sol))
        end
    end
    for permutation in get_perm(spg, :R)
        f_op = linear_function_for_spatial_operation(permutation)
        P_sol = find_subspace(T, P_sol, f_op; λ=reps[3], is_hermitian=false, tol=tol, _mapping_table=mt)
        if verbose
            println("operation R, size(P_sol) = ", size(P_sol))
        end
    end
    return P_sol
end

"""
    find_solution(spg::AbstractPointGroup, T::AbstractTensorMap, reps_name::Symbol; P_filter=nothing)

Return normalized tensor solutions transforming under `reps_name` of `spg`.
"""
function find_solution(spg::AbstractPointGroup, T::AbstractTensorMap, reps_name::Symbol; P_filter=nothing)
    P_sol = find_subspace(spg, T, reps_name; P_filter=P_filter)
    return find_solution(T, P_sol)
end
