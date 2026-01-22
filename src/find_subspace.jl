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
    find_subspace_matrixrep(T::AbstractTensorMap, P_init::Matrix{<:Number}, f_op::Function, D::AbstractMatrix;
        tol=1e-8, _mapping_table=mapping_table(T))

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
    if init_subspace_size == 0
        println("input subspace is empty")
        return P_init
    end
    if size(D, 1) != size(D, 2)
        throw(ArgumentError("representation matrix must be square"))
    end

    M_op = matrix_for_linear_function(T, f_op; _mapping_table=_mapping_table)

    k = size(P_init, 2)
    d = size(D, 1)
    Tprom = promote_type(eltype(M_op), eltype(D))
    D = Matrix{Tprom}(D)
    I_d = Matrix{Tprom}(I, d, d)
    MP = Matrix{Tprom}(M_op * P_init)
    Pp = Matrix{Tprom}(P_init)
    K = kron(I_d, MP) - kron(transpose(D), Pp)
    _, S, V = svd(K)
    idx = findall(s -> s <= tol, S)
    if isempty(idx)
        return P_init[:, 1:0]
    end

    cols = Matrix{Tprom}(undef, size(P_init, 1), 0)
    for j in idx
        X = reshape(V[:, j], k, d)
        cols = hcat(cols, P_init * X)
    end
    return cols
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
        #for ix in 1:init_subspace_size
        #    for iy in 1:init_subspace_size
        #        print(" ", (P_init' * P_init)[ix, iy], " ")
        #    end
        #    println(" ")
        #end
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
    # Orthonormalize columns in order (modified Gram-Schmidt) to preserve block structure.
    max_col_norm = maximum(norm(P_sol[:, ix]) for ix in 1:num_solutions)
    tol = max_col_norm * sqrt(eps(real(eltype(P_sol))))
    P_ortho = Matrix{eltype(P_sol)}(undef, size(P_sol, 1), 0)
    for ix in 1:num_solutions
        v = copy(P_sol[:, ix])
        for jx in 1:size(P_ortho, 2)
            v -= (P_ortho[:, jx]' * v) * P_ortho[:, jx]
        end
        nv = norm(v)
        if nv > tol
            P_ortho = hcat(P_ortho, v / nv)
        end
    end
    if size(P_ortho, 2) == 0
        return Vector{typeof(T)}(undef, 0)
    end
    sols = [set_data_by_vector(T, vec(P_ortho[:, ix]); _mapping_table=_mapping_table) for ix in 1:size(P_ortho, 2)]
    return [sol / norm(sol) for sol in sols]
end
