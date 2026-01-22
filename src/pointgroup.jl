"""
Abstract point group type for spatial symmetry operations.

Notes:
- Concrete point groups (e.g., `C4v`, `C3v`) should implement `get_reps` and
  `get_perm`.
"""
abstract type AbstractPointGroup end

"""
    find_subspace(spg::AbstractPointGroup, T::AbstractTensorMap, reps_name::Symbol; P_filter=nothing, verbose=false, tol=1e-8)

Find the subspace of `T` transforming under the representation `reps_name` of `spg`.

Notes:
- `reps_name` is interpreted by `get_reps(spg, reps_name)` and must be supported
  by the specific point group.
- `P_filter` is an optional projector in parameter space applied before symmetry
  constraints (e.g., to enforce occupation constraints).
- The method applies constraints sequentially for the operation families
  `:σd`, `:σv`, and `:R`.
"""
function find_subspace(spg::AbstractPointGroup, T::AbstractTensorMap, reps_name::Symbol; P_filter=nothing, verbose::Bool=false, tol::Real=1e-8)
    reps = get_reps(spg, reps_name)
    mt = mapping_table(T)
    num_paras = num_free_parameters(T; _mapping_table=mt)

    if isnothing(P_filter)
        P_filter = Matrix{ComplexF64}(I, num_paras, num_paras)
    end

    P_sol = P_filter
    d_irrep = nothing
    for rep in reps
        if rep isa AbstractVector && !isempty(rep) && rep[1] isa AbstractMatrix
            d_irrep = size(rep[1], 1)
            break
        elseif rep isa AbstractMatrix
            d_irrep = size(rep, 1)
            break
        end
    end
    if d_irrep !== nothing && size(P_sol, 1) == num_paras
        P_sol = vcat(ntuple(_ -> P_sol, d_irrep)...) # for high-dimensional irreps, stack P_sol d_irrep times
    end

    if verbose
        println("init, size(P_sol) = ", size(P_sol))
    end

    for (idx, permutation) in enumerate(get_perm(spg, :σd))
        f_op = linear_function_for_spatial_operation(permutation)
        rep = reps[1]
        rep_val = rep isa AbstractVector ? rep[idx] : rep # for 2d irreps, rep is a vector of matrices
        if rep_val isa AbstractMatrix # higher dimensional irreps
            d_irrep = size(rep_val, 1)
            if size(P_sol, 1) != num_paras * d_irrep
                throw(ArgumentError("P_filter has $(size(P_sol, 1)) rows, expected $(num_paras) or $(num_paras * d_irrep)"))
            end
            P_sol = find_subspace_matrixrep(T, P_sol, f_op, rep_val; tol=tol, _mapping_table=mt)
        else # 1D irreps
            if size(P_sol, 1) != num_paras
                throw(ArgumentError("P_filter has $(size(P_sol, 1)) rows, expected $(num_paras) for 1D irreps"))
            end
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
            d_irrep = size(rep_val, 1)
            if size(P_sol, 1) != num_paras * d_irrep
                throw(ArgumentError("P_filter has $(size(P_sol, 1)) rows, expected $(num_paras) or $(num_paras * d_irrep)"))
            end
            P_sol = find_subspace_matrixrep(T, P_sol, f_op, rep_val; tol=tol, _mapping_table=mt)
        else
            if size(P_sol, 1) != num_paras
                throw(ArgumentError("P_filter has $(size(P_sol, 1)) rows, expected $(num_paras) for 1D irreps"))
            end
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
            d_irrep = size(rep_val, 1)
            if size(P_sol, 1) != num_paras * d_irrep
                throw(ArgumentError("P_filter has $(size(P_sol, 1)) rows, expected $(num_paras) or $(num_paras * d_irrep)"))
            end
            P_sol = find_subspace_matrixrep(T, P_sol, f_op, rep_val; tol=tol, _mapping_table=mt)
        else
            if size(P_sol, 1) != num_paras
                throw(ArgumentError("P_filter has $(size(P_sol, 1)) rows, expected $(num_paras) for 1D irreps"))
            end
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

Notes:
- This is a convenience wrapper: `find_subspace` followed by `find_solution`.
- For 1D irreps, each solution column yields one tensor.
- For higher-dimensional irreps, each solution column yields `d_irrep` tensors
  by splitting the stacked parameter vector into `d_irrep` blocks.
- The returned vector is flat: "column" refers to a basis vector in `P_sol`,
  and "block" refers to one of the `d_irrep` consecutive `num_paras`-length
  slices of that column. Ordering is:
  (column 1, block 1), (column 1, block 2), ..., (column 1, block d_irrep),
  then (column 2, block 1), ..., (column 2, block d_irrep), etc.
- Returned tensors are normalized individually.
"""
function find_solution(spg::AbstractPointGroup, T::AbstractTensorMap, reps_name::Symbol; P_filter=nothing)
    P_sol = find_subspace(spg, T, reps_name; P_filter=P_filter)
    return find_solution(T, P_sol)
end
