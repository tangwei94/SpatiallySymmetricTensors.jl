"""
Abstract point group type for spatial symmetry operations.
"""
abstract type AbstractPointGroup end

"""
    projector_function(spg::AbstractPointGroup, reps_name::Symbol)

Return a function that applies the irrep projector:  P = d_reps / |G| ∑_{g ∈ G} conj(χ)(g) ρ(g)
"""
function projector_function(spg::AbstractPointGroup, reps_name::Symbol)
    reps = get_reps(spg, reps_name)
    op_names = (:σd, :σv, :R)

    perms_with_char = Tuple{Any, ComplexF64}[]
    for (i, name) in enumerate(op_names)
        χ = reps[i]
        χ == 0 && continue
        for perm in get_perm(spg, name)
            push!(perms_with_char, (perm, χ))
        end
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
