"""
    matrix_for_linear_function(T::AbstractTensorMap, f_op::Function; _mapping_table=mapping_table(T))

Return the matrix representation of `f_op` in the free-parameter basis of `T`.

Notes:
- Constructs the matrix by applying `f_op` to basis tensors obtained from unit
  vectors in parameter space.
- The conversion between the parameter space and the tensor space is defined by `_mapping_table`.
- This is potentially expensive (O(n^2) applications of `f_op`), so cache results
  when used repeatedly.
"""
function matrix_for_linear_function(T::AbstractTensorMap, f_op::Function; _mapping_table::MappingTable=mapping_table(T))
    num_paras = num_free_parameters(T; _mapping_table=_mapping_table)
    M = zeros(ComplexF64, num_paras, num_paras)
    for ix in 1:num_paras
        paras = zeros(ComplexF64, num_paras)
        paras[ix] = 1
    
        T1 = set_data_by_vector(T, paras; _mapping_table=_mapping_table)
        X = vec(f_op(T1))
        M[:, ix] = X 
    end
    return M
end

"""
    find_subspace(T::AbstractTensorMap, P_init::Matrix{<:Number}, f_op::Function; λ=1.0, is_hermitian=false, tol=1e-8, _mapping_table=mapping_table(T))

Project an initial subspace `P_init` onto the eigenspace of `f_op` with eigenvalue `λ`.
Returns the basis matrix for the resulting subspace.

Notes:
- `P_init` must have orthonormal columns to within `tol`.
- `is_hermitian=true` uses a Hermitian eigensolver on the projected operator.
- The returned basis is in the original parameter space; it is orthonormal when
  `is_hermitian=false` (via QR), and inherited otherwise.
"""
function find_subspace(T::AbstractTensorMap, P_init::Matrix{<:Number}, f_op::Function; λ::Real=1.0, is_hermitian::Bool=false, tol::Real=1e-8, _mapping_table::MappingTable=mapping_table(T))

    init_subspace_size = size(P_init, 2)
    if init_subspace_size == 0
        println("input subspace is empty")
        return P_init
    end

    if norm(P_init' * P_init - Matrix{eltype(P_init)}(I, init_subspace_size, init_subspace_size)) > tol
        error("input subspace is not orthonormal")
    end

    M_op = matrix_for_linear_function(T, f_op; _mapping_table=_mapping_table)
    M_op_sub = P_init' * M_op * P_init
    if ! is_hermitian
        Λop, Uop = eigen(M_op_sub)
        Pop = Uop[:, norm.(Λop .- λ) .< tol]
        if size(Pop, 2) == 0
            return P_init[:, 1:0]
        end
        Pop = Matrix(qr(Pop).Q) # orthonormalize the basis
    else
        Λop, Uop = eigen(Hermitian(M_op_sub))
        Pop = Uop[:, norm.(Λop .- λ) .< tol]
        if size(Pop, 2) == 0
            return P_init[:, 1:0]
        end
    end

    P_sol = P_init * Pop
    return P_sol
end

"""
    find_subspace_matrixrep(T::AbstractTensorMap, P_init::Matrix{<:Number}, f_op::Function, D::AbstractMatrix;
        tol=1e-8, _mapping_table=mapping_table(T))

FIXME. docstring for this function is outdated. need to be updated!

Project `P_init` onto the subspace transforming as matrix representation `D` under `f_op`.
Returns a basis matrix for the resulting subspace.

Notes:
- `P_init` columns are an initial basis in parameter space.
- Solves `(I ⊗ M_op) vec(X) - (D^T ⊗ I) vec(X) = 0 ` by looking for the null space of the LHS
- The result is a (not necessarily orthonormal) basis in the original parameter
  space.
"""
function find_subspace_matrixrep(T::AbstractTensorMap, P_init::Matrix{<:Number}, f_op::Function, D::AbstractMatrix;
    tol::Real=1e-8, _mapping_table::MappingTable=mapping_table(T))

    init_subspace_size = size(P_init, 2)
    d_irrep = size(D, 1)
    num_sol_init = size(P_init, 2)
    num_paras = num_free_parameters(T; _mapping_table=_mapping_table)

    if init_subspace_size == 0
        println("input subspace is empty")
        return P_init
    end
    if num_paras != size(P_init, 1) * d_irrep
        throw(ArgumentError("P_init should have $(num_paras)*$(d_irrep) rows"))
    end
    if size(D, 1) != size(D, 2)
        throw(ArgumentError("representation matrix must be square"))
    end

    M_op = matrix_for_linear_function(T, f_op; _mapping_table=_mapping_table)

    Tprom = promote_type(eltype(M_op), eltype(D))
    D = Matrix{Tprom}(D)
    I_paras = Matrix{Tprom}(I, num_paras, num_paras)

    K = zeros(Tprom, num_paras*d_irrep, num_sol_init*d_irrep)
    for ix in 1:d_irrep
        indices_L = (ix-1)*num_paras+(1:num_paras)
        indices_R = (ix-1)*num_sol_init+(1:num_sol_init)
        K[indices_L, indices_R] = P_init[indices_L, :]
    end
    K .= - kron(transpose(D), I_paras) * K
    for ix in 1:d_irrep
        indices_L = (ix-1)*num_paras+(1:num_paras)
        indices_R = (ix-1)*num_sol_init+(1:num_sol_init)
        K[indices_L, indices_R] .+= M_op * P_init[indices_L, :]
    end
    
    _, S, V = svd(K)
    idx = findall(s -> s <= tol, S)
    if isempty(idx)
        return P_init[:, 1:0]
    end

    return V[:, idx]
end

"""
    find_solution(T::AbstractTensorMap, P_sol::Matrix{<:Number}; _mapping_table=mapping_table(T))

Convert a solution subspace `P_sol` into normalized symmetric tensor solutions.

Notes:
- Columns of `P_sol` are interpreted as parameter vectors.
- Uses modified Gram-Schmidt to preserve the original block structure ordering.
- Each returned tensor is normalized by its Frobenius norm.
"""
function find_solution(T::AbstractTensorMap, P_sol::Matrix{<:Number}; _mapping_table::MappingTable=mapping_table(T))
    num_solutions = size(P_sol, 2)
    if num_solutions == 0
        return Vector{typeof(T)}(undef, 0)
    end
    num_paras = num_free_parameters(T; _mapping_table=_mapping_table)
    nrows = size(P_sol, 1)
    if nrows == num_paras
        F = qr(P_sol)
        P_ortho = Matrix(F.Q)[:, 1:num_solutions]
        sols = [set_data_by_vector(T, vec(P_ortho[:, ix]); _mapping_table=_mapping_table) for ix in 1:num_solutions]
        return [sol / norm(sol) for sol in sols]
    elseif nrows % num_paras == 0
        d_irrep = nrows ÷ num_paras
        F = qr(P_sol)
        P_ortho = Matrix(F.Q)[:, 1:num_solutions]
        sols = Vector{typeof(T)}(undef, num_solutions * d_irrep)
        idx = 1
        for ix in 1:num_solutions
            v = P_ortho[:, ix]
            for jx in 1:d_irrep
                block = v[(jx - 1) * num_paras + 1: jx * num_paras]
                sol = set_data_by_vector(T, vec(block); _mapping_table=_mapping_table)
                sols[idx] = sol / norm(sol)
                idx += 1
            end
        end
        return sols
    else
        throw(ArgumentError("P_sol has $(nrows) rows, which is not compatible with num_free_parameters(T)=$(num_paras)"))
    end
end
