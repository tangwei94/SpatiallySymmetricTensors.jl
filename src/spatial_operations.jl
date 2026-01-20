const MappingEntry = Tuple{Union{FusionTree,Nothing}, Union{FusionTree,Nothing}, Int, Int}
const MappingTable = Vector{MappingEntry}

"""
    mapping_table(T::AbstractTensorMap)

Construct a mapping table from free parameters to tensor blocks.
"""
function mapping_table(T::AbstractTensorMap)
    res = MappingEntry[]
    a = 1
    for (f1, f2) in fusiontrees(T)
        n = length(T[f1, f2])
        push!(res, (f1, f2, a, n))
        a += n
    end
    return res
end

"""
    num_free_parameters(T::AbstractTensorMap)

Return the number of free parameters in a symmetric tensor.
"""
function num_free_parameters(T::AbstractTensorMap; _mapping_table::MappingTable=mapping_table(T))
    return sum([n for (_, _, _, n) in _mapping_table])
end

"""
    set_data_by_vector!(T::AbstractTensorMap, _values; _mapping_table=mapping_table(T))

Set free parameters of `T` from `values` in the mapping-table order.
"""
function set_data_by_vector!(T::AbstractTensorMap, values::Vector{<:Number}; _mapping_table::MappingTable=mapping_table(T))
    for ix in eachindex(_mapping_table)
        f1, f2, a, n = _mapping_table[ix]
        reshape(T[f1, f2], n) .= values[a:a+n-1]
    end
    return T
end

"""
    set_data_by_vector(T::AbstractTensorMap, values::Vector{<:Number}; _mapping_table=mapping_table(T))

Return a new tensor with the same symmetry structure as `T`, with free parameters set by `values`.
"""
function set_data_by_vector(T::AbstractTensorMap, values::Vector{<:Number}; _mapping_table::MappingTable=mapping_table(T))
    T1 = similar(T) 
    set_data_by_vector!(T1, values; _mapping_table=_mapping_table)
    return T1
end

"""
    Base.vec(T::AbstractTensorMap; _mapping_table=mapping_table(T))

Return a vector of free parameters in the same order as `set_data_by_vector`.
"""
function Base.vec(T::AbstractTensorMap; _mapping_table::MappingTable=mapping_table(T))
    num_paras = num_free_parameters(T; _mapping_table=_mapping_table)
    v = Vector{eltype(T)}(undef, num_paras)
    offset = 1
    for (f1, f2, _, n) in _mapping_table
        v[offset:offset+n-1] = vec(T[f1, f2])
        offset += n
    end
    return v
end

"""
    selector(T::AbstractTensorMap, _condition; _mapping_table=mapping_table(T))

Construct a projector onto parameters that satisfy `condition(f1, f2)`.

`condition` takes fusion trees `(f1, f2)` and returns a Bool.
"""
function selector(T::AbstractTensorMap, condition::Function; _mapping_table::MappingTable=mapping_table(T))
    num_paras = num_free_parameters(T; _mapping_table=_mapping_table)
    indices = Int[]
    for (f1, f2, a, n) in _mapping_table
        if condition(f1, f2)
            append!(indices, a:a+n-1)
        end
    end
    P = zeros(ComplexF64, num_paras, length(indices))
    for (j, idx) in enumerate(indices)
        P[idx, j] = 1
    end
    return P
end

function linear_function_for_spatial_operation(permutation)
    return ((T1) -> permute(T1, permutation))
end
"""
    matrix_for_spatial_operation(T::AbstractTensorMap, permutation; _mapping_table=mapping_table(T))

Return the matrix representation of a spatial permutation acting on free parameters of `T`.
"""
function matrix_for_spatial_operation(T::AbstractTensorMap, permutation; _mapping_table::MappingTable=mapping_table(T))
    return matrix_for_linear_function(T, linear_function_for_spatial_operation(permutation); _mapping_table=_mapping_table)
end
