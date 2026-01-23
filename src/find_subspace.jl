"""
    matrix_for_linear_function(T::AbstractTensorMap, f_op::Function; _mapping_table=mapping_table(T))

Return the matrix representation of `f_op` in the free-parameter basis of `T`.
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
    find_subspace_from_projector(T::AbstractTensorMap, f_proj::Function; P_filter=nothing, tol=1e-8, _mapping_table=mapping_table(T))

Construct a symmetry subspace by QR-decomposing the projector matrix.
"""
function find_subspace_from_projector(T::AbstractTensorMap, f_proj::Function; P_filter=nothing, tol::Real=1e-8, _mapping_table::MappingTable=mapping_table(T))
    num_paras = num_free_parameters(T; _mapping_table=_mapping_table)
    if isnothing(P_filter)
        P_filter = Matrix{ComplexF64}(I, num_paras, num_paras)
    end

    init_subspace_size = size(P_filter, 2)
    if init_subspace_size == 0
        @warn "input subspace is empty, returning empty subspace"
        return P_filter
    end

    if norm(P_filter' * P_filter - Matrix{eltype(P_filter)}(I, init_subspace_size, init_subspace_size)) > tol
        throw(ArgumentError("input subspace is not orthonormal"))
    end

    M_op = matrix_for_linear_function(T, f_proj; _mapping_table=_mapping_table)
    M_sub = P_filter' * M_op * P_filter

    F = qr(M_sub)
    R = F.R
    keep = findall(abs.(diag(R)) .> tol)
    if isempty(keep)
        @warn "output subspace is empty, returning empty subspace"
        return P_filter[:, 1:0]
    end

    Q = Matrix(F.Q)
    return P_filter * Q[:, keep]
end

"""
    find_subspace(T::AbstractTensorMap, P_init::Matrix{<:Number}, f_op::Function; λ=1.0, is_hermitian=false, tol=1e-8, _mapping_table=mapping_table(T))

Project an initial subspace `P_init` onto the eigenspace of `f_op` with eigenvalue `λ`.
Returns the basis matrix for the resulting subspace.
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
    find_solution(T::AbstractTensorMap, P_sol::Matrix{<:Number}; _mapping_table=mapping_table(T))

Convert a solution subspace `P_sol` into normalized symmetric tensor solutions.
"""
function find_solution(T::AbstractTensorMap, P_sol::Matrix{<:Number}; _mapping_table::MappingTable=mapping_table(T))
    num_solutions = size(P_sol, 2)
    sols = [set_data_by_vector(T, vec(P_sol[:, ix]); _mapping_table=_mapping_table) for ix in 1:num_solutions]
    return [sol / norm(sol) for sol in sols]
end
