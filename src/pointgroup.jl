abstract type AbstractPointGroup end

function find_subspace(spg::AbstractPointGroup, T::AbstractTensorMap, reps_name::Symbol; P_filter=nothing)
    reps = get_reps(spg, reps_name)
    mt = mapping_table(T)

    if isnothing(P_filter)
        num_paras = num_free_parameters(T; _mapping_table=mt)
        P_filter = Matrix{ComplexF64}(I, num_paras, num_paras)
    end

    P_sol = P_filter

    println("init, size(P_sol) = ", size(P_sol))
    for permutation in get_perm(spg, :σd)
        f_op = linear_function_for_spatial_operation(permutation)
        P_sol = find_subspace(T, P_sol, f_op; λ = reps[1], is_hermitian=true, _mapping_table=mt)
        println("operation σd, size(P_sol) = ", size(P_sol))
    end
    for permutation in get_perm(spg, :σv)
        f_op = linear_function_for_spatial_operation(permutation)
        P_sol = find_subspace(T, P_sol, f_op; λ = reps[2], is_hermitian=true, _mapping_table=mt)
        println("operation σv, size(P_sol) = ", size(P_sol))
    end
    for permutation in get_perm(spg, :R)
        f_op = linear_function_for_spatial_operation(permutation)
        P_sol = find_subspace(T, P_sol, f_op; λ = reps[3], is_hermitian=false, _mapping_table=mt)
        println("operation R, size(P_sol) = ", size(P_sol))
    end
    return P_sol
end

function find_solution(spg::AbstractPointGroup, T::AbstractTensorMap, reps_name::Symbol; P_filter=nothing)
    P_sol = find_subspace(spg, T, reps_name; P_filter=P_filter)
    return find_solution(T, P_sol)
end
