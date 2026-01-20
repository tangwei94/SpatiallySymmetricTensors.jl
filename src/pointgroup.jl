"""
Abstract point group type for spatial symmetry operations.
"""
abstract type AbstractPointGroup end

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
    for (idx, permutation) in enumerate(get_perm(spg, :σd))
        f_op = linear_function_for_spatial_operation(permutation)
        rep = reps[1]
        rep_val = rep isa AbstractVector ? rep[idx] : rep
        if rep_val isa AbstractMatrix
            P_sol = find_subspace_matrixrep(T, P_sol, f_op, rep_val; tol=tol, _mapping_table=mt)
        else
            P_sol = find_subspace(T, P_sol, f_op; λ=rep_val, is_hermitian=true, tol=tol, _mapping_table=mt)
        end
        if verbose
            println("operation σd, size(P_sol) = ", size(P_sol))
        end
    end
    for (idx, permutation) in enumerate(get_perm(spg, :σv))
        f_op = linear_function_for_spatial_operation(permutation)
        rep = reps[2]
        rep_val = rep isa AbstractVector ? rep[idx] : rep
        if rep_val isa AbstractMatrix
            P_sol = find_subspace_matrixrep(T, P_sol, f_op, rep_val; tol=tol, _mapping_table=mt)
        else
            P_sol = find_subspace(T, P_sol, f_op; λ=rep_val, is_hermitian=true, tol=tol, _mapping_table=mt)
        end
        if verbose
            println("operation σv, size(P_sol) = ", size(P_sol))
        end
    end
    for (idx, permutation) in enumerate(get_perm(spg, :R))
        f_op = linear_function_for_spatial_operation(permutation)
        rep = reps[3]
        rep_val = rep isa AbstractVector ? rep[idx] : rep
        if rep_val isa AbstractMatrix
            P_sol = find_subspace_matrixrep(T, P_sol, f_op, rep_val; tol=tol, _mapping_table=mt)
        else
            P_sol = find_subspace(T, P_sol, f_op; λ=rep_val, is_hermitian=false, tol=tol, _mapping_table=mt)
        end
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
