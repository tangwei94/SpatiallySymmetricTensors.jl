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

function find_subspace(T::AbstractTensorMap, P_init::Matrix{<:Number}, f_op::Function; λ::Real=1.0, is_hermitian::Bool=false, _mapping_table::MappingTable=mapping_table(T))

    init_subspace_size = size(P_init, 2)
    if init_subspace_size == 0
        println("input subspace is empty")
        return P_init
    end

    if norm(P_init' * P_init - Matrix{eltype(P_init)}(I, init_subspace_size, init_subspace_size)) > 1e-8
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
        Pop = Uop[:, norm.(Λop .- λ) .< 1e-8]
        Pop = Matrix(qr(Pop).Q) # orthonormalize the basis
    else
        Λop, Uop = eigen(Hermitian(M_op_sub))
        Pop = Uop[:, norm.(Λop .- λ) .< 1e-8]
    end

    P_sol = P_init * Pop
    return P_sol
end

function find_solution(T::AbstractTensorMap, P_sol::Matrix{<:Number}; _mapping_table::MappingTable=mapping_table(T))
    num_solutions = size(P_sol, 2)
    sols = [set_data_by_vector(T, vec(P_sol[:, ix]); _mapping_table=_mapping_table) for ix in 1:num_solutions]
    return [sol / norm(sol) for sol in sols]
end