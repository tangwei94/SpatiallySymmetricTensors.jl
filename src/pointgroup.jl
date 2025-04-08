abstract type AbstractPointGroup end

function find_solution(spg::AbstractPointGroup, T::AbstractTensorMap, reps_name::Symbol; P_filter=nothing)
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

    # solution 
    num_solutions = size(P_sol, 2)
    sols = [set_data_by_vector(T, vec(P_sol[:, ix]); _mapping_table=mt) for ix in 1:num_solutions]

    #for sol in sols, perm in [get_perm(spg, :σd); get_perm(spg, :σv); get_perm(spg, :R)] 
    #    @show norm(sol - permute(sol, perm))
    #end

    return [sol / norm(sol) for sol in sols]
end