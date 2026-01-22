const MappingEntry = Tuple{Union{FusionTree,Nothing}, Union{FusionTree,Nothing}, Int, Int}
const MappingTable = Vector{MappingEntry}

"""
    mapping_table(T::AbstractTensorMap)

Construct a mapping table from free parameters to tensor blocks.

Notes:
- Each entry is `(f1, f2, a, n)` where `(f1, f2)` are fusion trees for a block,
  `a` is the starting linear index (1-based), and `n` is the number of entries
  in the block `T[f1, f2]`.
- The order follows `fusiontrees(T)` and defines the parameter ordering used by
  `set_data_by_vector`, `set_data_by_vector!`, and `vec`.
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

Notes:
- Equivalent to the total length of all blocks in `mapping_table(T)`.
- If you already have a mapping table, pass it via `_mapping_table` to avoid
  recomputation.
"""
function num_free_parameters(T::AbstractTensorMap; _mapping_table::MappingTable=mapping_table(T))
    return sum([n for (_, _, _, n) in _mapping_table])
end

"""
    set_data_by_vector!(T::AbstractTensorMap, values::Vector{<:Number}; _mapping_table=mapping_table(T))

Set free parameters of `T` from `values` in the mapping-table order. (parameter space -> tensor space)

Notes:
- `values` must have length equal to `num_free_parameters(T)`.
- The write is block-wise: each tensor block is reshaped to a vector and filled
  from the corresponding slice of `values`.
- Mutates `T` and returns it.
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

Return a new tensor with the same symmetry structure as `T`, with free parameters set by `values`. (tensor space -> parameter space)

Notes:
- Uses `similar(T)` to preserve symmetry structure and block layout.
- `values` must match the length implied by `_mapping_table`.
- For in-place updates, use `set_data_by_vector!`.
"""
function set_data_by_vector(T::AbstractTensorMap, values::Vector{<:Number}; _mapping_table::MappingTable=mapping_table(T))
    T1 = similar(T) 
    set_data_by_vector!(T1, values; _mapping_table=_mapping_table)
    return T1
end

"""
    Base.vec(T::AbstractTensorMap; _mapping_table=mapping_table(T))

Return a vector of free parameters in the same order as `set_data_by_vector`. (tensor space -> parameter space)

Notes:
- The ordering is defined by `_mapping_table` (or `mapping_table(T)`).
- This is the inverse operation to `set_data_by_vector` up to reshaping.
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
    selector(T::AbstractTensorMap, condition::Function; _mapping_table=mapping_table(T))

Construct a projector onto parameters that satisfy `condition(f1, f2)`.

`condition` takes fusion trees `(f1, f2)` and returns a Bool.

Notes:
- Returns a dense matrix `P` such that `P' * v` selects the subset of parameters
  whose blocks satisfy `condition`.
- The columns of `P` are standard basis vectors picking the selected indices.
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

Notes:
- `permutation` should be compatible with `permute(T, permutation)` in TensorKit.
- The resulting matrix acts on the parameter vector returned by `vec(T)`.
- For repeated use, cache `_mapping_table` to avoid recomputation.
"""
function matrix_for_spatial_operation(T::AbstractTensorMap, permutation; _mapping_table::MappingTable=mapping_table(T))
    return matrix_for_linear_function(T, linear_function_for_spatial_operation(permutation); _mapping_table=_mapping_table)
end
