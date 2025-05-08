const U1GradedSpace = GradedSpace{U1Irrep, TensorKit.SortedVectorDict{U1Irrep, Int64}}

function u1_charge_conjugation(f::FusionTree{U1Irrep})
    uncoupled = Tuple(U1Irrep(-x.charge) for x in f.uncoupled)
    coupled = U1Irrep(-f.coupled.charge) 
    return FusionTree(uncoupled, coupled, f.isdual, f.innerlines, f.vertices)
end

function u1_charge_conjugation(T::AbstractTensorMap{N, U1GradedSpace}) where N
    T_mapped = deepcopy(T)
    for (f1, f2) in fusiontrees(T)
        f1_mapped = u1_charge_conjugation(f1)
        f2_mapped = u1_charge_conjugation(f2)

        if size(T[f1, f2]) != size(T_mapped[f1_mapped, f2_mapped])
            @error "size mismatch" f1, f2
        end
        T_mapped[f1_mapped, f2_mapped] .= T[f1, f2]
    end
    return T_mapped
end
  
function find_subspace_for_u1_charge_conjugation(T::AbstractTensorMap{N, U1GradedSpace}, P_init::Matrix{<:Number}; λ::Real=1.0, _mapping_table::MappingTable=mapping_table(T)) where N
    return find_subspace(T, P_init, u1_charge_conjugation; λ=λ, _mapping_table=_mapping_table)
end
