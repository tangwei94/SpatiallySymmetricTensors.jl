const U1GradedSpace = GradedSpace{U1Irrep, TensorKit.SortedVectorDict{U1Irrep, Int64}}

"""
    u1_charge_conjugation(f::FusionTree{U1Irrep})

Return the charge-conjugated fusion tree.

Notes:
- Flips the sign of each U(1) charge in both uncoupled and coupled irreps.
- Preserves dual flags and the internal fusion-tree structure.
"""
function u1_charge_conjugation(f::FusionTree{U1Irrep})
    uncoupled = Tuple(U1Irrep(-x.charge) for x in f.uncoupled)
    coupled = U1Irrep(-f.coupled.charge) 
    return FusionTree(uncoupled, coupled, f.isdual, f.innerlines, f.vertices)
end

"""
    u1_charge_conjugation(T::AbstractTensorMap{N, U1GradedSpace}) where N

Apply charge conjugation to a U(1)-graded tensor map.

Notes:
- Each block `(f1, f2)` is copied into the block indexed by the charge-conjugated
  fusion trees `(f1_mapped, f2_mapped)`.
- Throws an `ArgumentError` if the conjugated block is missing or has a size
  mismatch.
- Returns a new tensor; the input tensor is not mutated.
"""
function u1_charge_conjugation(T::AbstractTensorMap{N, U1GradedSpace}) where N
    T_mapped = deepcopy(T)
    for (f1, f2) in fusiontrees(T)
        f1_mapped = u1_charge_conjugation(f1)
        f2_mapped = u1_charge_conjugation(f2)

        mapped_block = try
            T_mapped[f1_mapped, f2_mapped]
        catch err
            throw(ArgumentError("missing charge-conjugated block for $(f1), $(f2): $(err)"))
        end
        if size(T[f1, f2]) != size(mapped_block)
            throw(ArgumentError("size mismatch for charge-conjugated block $(f1), $(f2)"))
        end
        mapped_block .= T[f1, f2]
    end
    return T_mapped
end
  
"""
    find_subspace_for_u1_charge_conjugation(T::AbstractTensorMap{N, U1GradedSpace}, P_init::Matrix{<:Number}; λ=1.0, _mapping_table=mapping_table(T)) where N

Project `P_init` onto the charge-conjugation eigenspace with eigenvalue `λ`.

Notes:
- Equivalent to `find_subspace(T, P_init, u1_charge_conjugation; λ=λ, ...)`.
- Use `λ=1` for charge-conjugation invariance and `λ=-1` for odd parity.
"""
function find_subspace_for_u1_charge_conjugation(T::AbstractTensorMap{N, U1GradedSpace}, P_init::Matrix{<:Number}; λ::Real=1.0, _mapping_table::MappingTable=mapping_table(T)) where N
    return find_subspace(T, P_init, u1_charge_conjugation; λ=λ, _mapping_table=_mapping_table)
end
