const MPOTensor{T,S} = AbstractTensorMap{T,S,2,2} where {T,S}

"""
    mpo_ovlp(A1, A2)

Return the dominant eigenpair of the MPO overlap transfer operator.

Notes:
- Constructs the transfer map using contraction of `A1` with `conj(A2)`.
- Uses `KrylovKit.eigsolve` with a random initial vector; results can depend
  on the seed if the dominant eigenvalue is degenerate.
- Returns `(eigenvalues, eigenvectors)` as in `eigsolve`.
"""
function mpo_ovlp(A1, A2)
    V1 = space(A1, 1)
    V2 = space(A2, 1)

    function mpo_transf(v)
        @tensor Tv[-1; -2] := A1[-1 3; 4 1] * conj(A2[-2 3; 4 2]) * v[1; 2]
        return Tv
    end

    v0 = rand(ComplexF64, V1, V2)
    return eigsolve(mpo_transf, v0, 1, :LM)
end

"""
    mpotensor_dag(T::MPOTensor)

Return the Hermitian conjugate of an MPO tensor.

Notes:
- The data is conjugated and the physical legs are swapped to implement the
  Hermitian adjoint at the tensor level.
- The returned tensor uses the same `space(T)` metadata as the input.
"""
function mpotensor_dag(T::MPOTensor)
    T_data = reshape(T.data, (dims(codomain(T))..., dims(domain(T))...))
    Tdag_data = permutedims(conj.(T_data), (1, 3, 2, 4))
    
    return TensorMap(Tdag_data, space(T))
end

"""
    mpo_hermicity(A)

Return a norm ratio diagnosing Hermiticity of MPO `A`.

Notes:
- Computes dominant eigenvalues of `A * A†` and `A * A` transfer operators.
- The ratio is 1 for exact Hermiticity (up to numerical error).
- Uses a random initial vector for `eigsolve`.
"""
function mpo_hermicity(A)
    v_space = space(A, 1)
    function AA_transf(v)
        @tensor Tv[-1; -2] := A[-1 3; 4 1] * conj(A[-2 3; 4 2]) * v[1; 2]
        return Tv
    end
    function AĀ_transf(v)
        @tensor Tv[-1 -2] := A[-1 3; 4 1] * A[-2 4; 3 2] * v[1 2]
        return Tv
    end

    vaa = rand(ComplexF64, v_space, v_space)
    vaā = Tensor(rand, ComplexF64, v_space*v_space)

    aa, _ = eigsolve(AA_transf, vaa, 1, :LM) 
    aā, _ = eigsolve(AĀ_transf, vaā, 1, :LM) 

    aa = aa[1]
    aā = aā[1]

    return norm(aā / aa) 
end

"""
    mpo_normality(A)

Return a norm ratio diagnosing normality of MPO `A`.

Notes:
- Compares dominant eigenvalues of `A† A A† A` and `A† A A A†` transfer maps.
- The ratio is 1 for exact normality (up to numerical error).
- Uses a random initial vector for `eigsolve`.
"""
function mpo_normality(A)
    v_space = space(A, 1)

    function ĀAĀA_transf(v)
        @tensor Tv[-1 -2; -3 -4] := conj(A[-3 7; 8 6]) * A[-1 7; 5 4] * conj(A[-4 3; 5 2]) * A[-2 3; 8 1] * v[1 4; 2 6]
        return Tv
    end
    function ĀAAĀ_transf(v)
        @tensor Tv[-1 -2; -3 -4] := conj(A[-3 7; 8 6]) * A[-1 7; 5 4] * A[-2 5; 3 2] * conj(A[-4 8; 3 1]) * v[4 2; 6 1]
        return Tv
    end
    
    v0 = rand(ComplexF64, v_space*v_space, v_space*v_space)

    āaāa, _ = eigsolve(ĀAĀA_transf, v0, 1, :LM) 
    āaaā, _ = eigsolve(ĀAAĀ_transf, v0, 1, :LM) 

    āaāa = āaāa[1]
    āaaā = āaaā[1]

    return norm(āaāa / āaaā) 
end
