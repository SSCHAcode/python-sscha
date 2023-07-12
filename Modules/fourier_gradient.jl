using LinearAlgebra
#using LoopVectorization

@doc raw"""
    function get_gradient_fourier!(
        Φ_grad :: Array{Complex{T}, 3},
        Φ_grad_err :: Array{Complex{T}, 3},
        v_tilde :: Array{Complex{T}, 3},
        δf_tilde :: Array{Complex{T}, 3},
        q_opposite_index :: Vector{I})  {where T <: AbstractFloat}

Also a nonbang (!) version of the subroutine exists where the Φ_grad and Φ_grad_err 
are returned instead of beying taken as parameters.

Compute the gradient in fourier space.


$$
    G_\Phi(q) = < Yu(q) δf(-q)>
$$

where


$$
    u_a(q) = \sum_R e^{-i q ⋅ R} u_a(R)
$$

All the average are reported by summing times the weights
without normalization. 
This makes it easier to enable parallelization in the error 
calculation and compute the actual average in a super function. 

## Parameters 

    - Φ_grad : size (nq, 3*nat, 3*nat); the gradient of the free energy 
    - Φ_grad_err : size (nq, 3*nat, 3*nat); the average of squares. Used to evaluate the variance 
    - v_tilde : size(n_random, 3*nat, nq); The result of multiplying Y with the displacements in q space
    - δf_tilde : size (n_random, 3*nat, nq); the forces (Ry / Bohr)
    - R_lat : size (nat_sc, 3) the lattice vectors
    - q : size (nq, 3) the q vectors
    - itau : size (nat_sc); the index in the primitive cell corresponding to the supercell
    - weights : size (n_random); the weights of the configurations
    - q_opposite_index : size(nq) The index of the q = -q + G


$$
    Y_{ab} = \sqrt{m_a m_b} \sum_{\mu} e_{\mu}^a e_{\mu}^b \frac{2\omega_\mu}{\hbar (2 n_\mu + 1)}
$$

where $\omega_\mu$ and $e_\mu$ are the frequencies and eigenmodes.

The method returns the gredient (not symmetrized) in fourier space
"""
function get_gradient_fourier!(Φ_grad :: Array{Complex{T}, 3}, 
        Φ_grad_err :: Array{T, 3},
        v_tilde :: Array{Complex{T}, 3},
        δf_tilde :: Array{Complex{T}, 3},
        weights :: Vector{T},
        minus_q_index :: Vector{N}) where {T <: AbstractFloat, I <: Integer, N <: Integer}
    
    # Get the system size
    nat_sc = size(u_sc, 2) ÷ 3
    n_random = size(u_sc, 1)
    nq = size(q, 2)
    nat = nat_sc ÷ nq

    # Compute the gradient exploiting a new average
    Φ_grad .= 0
    Φ_grad_err .= 0
    # Define helper arrays that avoid memory allocations during the for loop
    
    @time begin 
        for jq in 1:nq
            for k in 1:3*nat
                for j in 1:3*nat
                    for i in 1:n_random
                        tmp = v_tilde[i, j, jq] * conj(δf_tilde[i, k, jq])
                        Φ_grad[j, k, jq] += tmp * weights[i]
                        Φ_grad_err[j, k, jq] += tmp * conj(tmp) * weights[i]

                        # @views mul!(tmp, v_tilde[:, jq, i], δf_tilde[:, jq, i]')  
                        # @. tmp2 = tmp * conj(tmp)

                        # tmp .*= weights[i]
                        # tmp2 .*= weights[i]

                        # Φ_grad[:, :, jq] .+= tmp
                        # Φ_grad_err[:, :, jq] .+= tmp2
                    end
                end
            end
        end
    end

    # Apply the transpose symmetry
    #@time begin 
    #    minus_q_index = zeros(Int, nq)
    #    get_opposit_q!(minus_q_index, q, bg)
    #end

    @time begin
        tmp_grad = zeros(Complex{T}, (3*nat, 3*nat, nq))
        for iq in 1:nq
            @views Φ_grad[:, :, iq] .+= Φ_grad[:, :, iq]'
            @views tmp_grad[:, :, iq] .= Φ_grad[:, :, iq]
        end
        for iq in 1:nq
            @views tmp_grad[:, :, iq] .+= conj.(Φ_grad[:, :, minus_q_index[iq]]')
        end
        Φ_grad .= tmp_grad
        Φ_grad .*= 0.25
    end
end
function get_gradient_fourier(
    v_tilde :: Array{Complex{T}, 3},
    δf_tilde :: Array{Complex{T}, 3},
    weights :: Vector{T},
    minus_q_index :: Vector{N}u_sc :: Matrix{T}
    ) where {T <: AbstractFloat, I <: Integer, N <: Integer}
    
    nq = size(q, 1)
    nat = size(u_sc, 2) ÷ (3*nq) 
    nat_sc = nat * nq

    Φ_grad = zeros(Complex{T}, (nq, 3*nat, 3*nat))
    Φ_grad_err = zeros(T, (nq, 3*nat, 3*nat))
    Φ_grad_jj = zeros(Complex{T}, (3*nat, 3*nat, nq))
    Φ_grad_err_jj = zeros(T, (3*nat, 3*nat, nq))

    get_gradient_fourier!(Φ_grad_jj, Φ_grad_err_jj,
        v_tilde, δf_tilde, weights, minus_q_vector)

    # Reorder the gradient
    for iq in 1:nq
        @views Φ_grad[iq, :, :] = Φ_grad_jj[:, :, iq]
        @views Φ_grad_err[iq, :, :] = Φ_grad_err_jj[:, :, iq]
    end

    return Φ_grad, Φ_grad_err
end


@doc raw"""
    get_opposite_q!(
        opposite_index :: Vector{Int},
        q_list :: Matrix{T},
        bg :: Matrix{T}; far :: Int = 3) where {T <: AbstractFloat}


Identify the opposite q in a list of q points.
The opposite_index array is filled for each q what is its inverse:

```
q_list[i, :] and q_list[opposite_index[i], :]
```

Are the opposite q vectors (q = -q + G)

## Parameters

- opposite_index : The array identifying the opposite
- q_list : (nq, 3) The list of q points
- bg : (3,3) The brilluin zone vector 

bg must be aligned in python way ``bg[i, :]`` is the i-th vector.
"""
function get_opposite_q!(
    opposite_index :: Vector{Int},
    q_list :: Matrix{T},
    bg :: Matrix{T}; far :: Int = 3) where {T <: AbstractFloat}

    opposite_index .= -1
    nq = size(q_list, 1)
    q_vector = zeros(T, 3)
    minus_q = zeros(T, 3)
    G1_tmp = zeros(T, 3)
    G2_tmp = zeros(T, 3)
    G3_tmp = zeros(T, 3)
    δ = zeros(T, 3)
    for i in 1:nq
        @views minus_q .= -q_list[i, :] 
        for j in 1:nq
            skip = false
            for xx in -far : far + 1
                if skip 
                    break
                end
                @views G1_tmp .= bg[1, :]
                G1_tmp .*= xx
                for yy in -far : far + 1
                    if skip 
                        break
                    end
                        @views G2_tmp .= bg[2, :]
                    G2_tmp .*= yy

                    for zz in -far : far + 1
                        @views G3_tmp .= bg[3, :]
                        G3_tmp .*= zz

                        begin 
                            @views q_vector .= q_list[j, :]
                            @views q_vector .+= G1_tmp
                            @views q_vector .+= G2_tmp
                            @views q_vector .+= G3_tmp 
                        end
                        δ .= minus_q
                        δ .-= q_vector
                        δ .= abs.(δ)

                        if max(δ...) < 1e-6
                            opposite_index[i] = j
                            skip = true
                            break
                        end
                    end
                end
            end
        end
    end    

    # Check if all the q have been found
    for i in 1:nq
        if opposite_index[i] == -1
            println("Error on q: $(q_list[i, :])")
            println("All q: $(opposite_index)")

            error("The opposite q for $i has not been found")
        end
    end
end
function get_opposite_q(
    q_list :: Matrix{T},
    bg :: Matrix{T}; far :: Int = 3) :: Vector{Int} where {T <: AbstractFloat}

    opposite_index = zeros(Int, size(q_list, 1))
    get_opposite_q!(opposite_index, q_list, bg; far = far)
    return opposite_index
end


@doc raw"""
    get_upsilon_fourier!(
        Y :: Array{Complex{T}, 3},
        ω :: Matrix{T},
        pols :: Array{Complex{T}, 3},
        masses :: Vector{T},
        temperature :: T;
        w_min_threshold :: T = 1e-6) where {T <: AbstractFloat}



Compute the Υ matrix in the Fourier space.

## Parameters

- Y : The Υ matrix in Fourier space (last index is the q point)
    This is modified by the function
- ω_q : The frequencies in the q point (q is the last index)
- pols_q : The polarization vectors in the q point (q is the last index)
- masses : (length nat) The masses of the atoms (in the primitive cell)
- temperature : The temperature in kelvin
- w_min_threshold : The minimum frequency to consider a mode (in Ry)
"""
function get_upsilon_fourier!(
    Y :: Array{Complex{T}, 3},
    ω_q :: Matrix{T},
    pols_q :: Array{Complex{T}, 3},
    masses :: Vector{T},
    temperature :: T;
    w_min_threshold :: T = 1e-12) where {T <: AbstractFloat}

    K_TO_RY = 6.336857346553283e-06

    nat = size(Y, 1) ÷ 3
    nq = size(Y, 3)

    dyn_mat = zeros(Complex{T}, (3*nat, 3*nat))
    n_ω = zeros(T, 3*nat)
    ω = zeros(T, 3*nat)
    pols = zeros(Complex{T}, (3*nat, 3*nat))
    factor = zeros(T, 3*nat)
    new_pols = zeros(Complex{T}, (3*nat, 3*nat))

    # Convert the temperature in Ry
    temperature = temperature * K_TO_RY
    skip_mode = zeros(Bool, 3*nat)

    for i in 1:nq
        @views ω .= ω_q[:, i]
        @views pols .= pols_q[:, :, i]
        
        # Setup the dynamical matrix
        skip_mode .= false
        for μ in 1:3*nat
            if ω[μ] < w_min_threshold
                ω[μ] = 0
                skip_mode[μ] = true
            end
        end

        # Compute the occupation factor
        if temperature > 0
            @. n_ω = 1.0 / (exp(ω / temperature) - 1.0)
        end

        @. factor = 2 * ω / (1.0 + 2*n_ω)
        
        for μ in 1:3*nat
            if skip_mode[μ]
                factor[μ] = 0
            end
        end
        
        # Get the correct Y matrix in q space
        for j1 in 1:3*nat
            for j2 in 1:3*nat
                new_pols[j1, j2] = pols[j1, j2] * factor[j2] 
            end
        end
        @views mul!(Y[:, :, i], new_pols, pols')

        # Apply the mass rescaling
        for j1 in 1:nat
            for j2 in 1:nat
                @views Y[3*(j1 - 1) + 1 : 3*j1, 3*(j2-1) + 1 : 3*j2, i] .*= √(masses[j1] * masses[j2])
            end
        end
    end
end
function get_upsilon_fourier(
    ω_q :: Matrix{T},
    pols_q :: Array{Complex{T}, 3},
    masses :: Vector{T},
    temperature :: T;
    w_min_threshold :: T = 1e-12) :: Array{Complex{T}, 3} where {T <: AbstractFloat}

    # Create Y and return it
    nq = size(ω_q, 2)
    nat = size(masses, 1)
    Y = zeros(Complex{T}, (3*nat, 3*nat, nq))
    get_upsilon_fourier!(Y, ω_q, pols_q, masses, temperature; w_min_threshold = w_min_threshold)
    return Y
end



@doc raw"""
    vector_r2q!(
        v_q :: Array{Complex{T}, 3},
        v_sc :: Matrix{T},
        q_tot :: Matrix{T})


Fourier transform a vector from real space and q space.

$$
v_k(\vec q) = \frac{1}{\sqrt{N_q}} \sum_{R} e^{-i 2\pi \vec R\cdot \vec q} v_k(\vec R)
$$


## Parameters

- v_q : (n_configs, 3nat, nq) 
    The target vector in Fourier space.
- v_sc : (n_configs, 3*nat_sc)
    The original vector in real space
- q_tot : (3, nq)
    The list of q vectors
- itau : (nat_sc)
    The correspondance for each atom in the supercell with the atom in the primitive cell.
- R_lat : (3, nat_sc)
    The origin coordinates of the supercell in which the atom is
"""
function vector_r2q!(
        v_q :: Array{Complex{T}, 3},
        v_sc :: Matrix{T},
        q :: Matrix{T},
        itau :: Vector{I},
        R_lat :: Matrix{T}
    ) where {T <: AbstractFloat, I <: Integer}

    nq = size(q, 2)
    n_random = size(v_sc, 1)
    nat_sc = size(v_sc, 2) ÷ 3
    nat = size(v_q, 2)

    v_q .= 0

    for jq ∈ 1:nq
        for k ∈ 1:nat_sc
            @views q_dot_R = q[:, jq]' * R_lat[:, k]
            exp_value = exp(- 1im * 2π * q_dot_R)

            for α in 1:3
                index_sc = 3 * (k - 1) + α
                index_uc = 3 * (itau[k] - 1) + α
                for i ∈ 1:n_random
                    v_q[i, index_uc, jq] += exp_value * v_sc[i, index_sc]
                end
            end
        end
    end

    v_q ./= √(nq)
end
function vector_r2q(
        v_sc :: Matrix{T},
        q_jj:: Matrix{T},
        itau :: Vector{I},
        R_jj :: Matrix{T}
    ) :: Array{Complex{T}, 3} where {T <: AbstractFloat, I <: Integer}

    nq = size(q_jj, 1)
    nat_sc = size(R_jj, 1)
    nat = nat_sc ÷ nq
    n_random = size(v_sc, 1)

    # Prepare the input and form the output
    R_lat = zeros(T, (3, nat_sc))
    q = zeros(T, (3, nq))

    R_lat .= R_jj'
    q .= q_jj'

    v_q = zeros(Complex{T}, (n_random, nat*3, nq))

    println("Timing vector_r2q! ")
    vector_r2q!(v_q, v_sc, q, itau, R_lat)
    return v_q
end


@doc raw"""
    vector_q2r!(
        v_sc :: Matrix{T},
        v_q :: Array{Complex{T}, 3},
        q_tot :: Matrix{T},
        itau :: Vector{I},
        R_lat :: Matrix{T}) where {T <: AbstractFloat, I <: Integer}


Fourier transform a vector from q space to real space.

$$
v_k(\vec R) = \frac{1}{\sqrt{N_q}} \sum_{R} e^{+i 2\pi \vec R\cdot \vec q} v_k(\vec q)
$$


## Parameters


- v_sc : (n_configs, 3*nat_sc)
    The target vector in real space
- v_q : (n_configs, nq, 3*nat) 
    The original vector in Fourier space. 
- q_tot : (3, nq)
    The list of q vectors
- itau : (nat_sc)
    The correspondance for each atom in the supercell with the atom in the primitive cell.
- R_lat : (3, nat_sc)
    The origin coordinates of the supercell in which the atom is
"""
function vector_q2r!(
        v_sc :: Matrix{T},
        v_q :: Array{Complex{T}, 3},
        q :: Matrix{T},
        itau :: Vector{I},
        R_lat :: Matrix{T}
    ) where {T <: AbstractFloat, I <: Integer}

    nq = size(q, 2)
    n_random = size(v_sc, 1)
    nat_sc = size(v_sc, 2) ÷ 3
    tmp_vector = zeros(Complex{T}, (n_random, 3*nat_sc))

    v_sc .= 0
    for jq ∈ 1:nq
        for k ∈ 1:nat_sc
            @views q_dot_R = q[:, jq]' * R_lat[:, k]
            exp_value = exp(1im * 2π * q_dot_R)

            for α in 1:3
                index_sc = 3 * (k - 1) + α
                index_uc = 3 * (itau[k] - 1) + α
                for i ∈ 1:n_random
                    tmp_vector[i, index_sc] += exp_value * v_q[i, index_uc, jq]
                end
            end
        end
    end

    v_sc .= real(tmp_vector)
    v_sc ./= √(nq)
end
function vector_q2r(
        v_q :: Array{Complex{T},3},
        q_jj:: Matrix{T},
        itau :: Vector{I},
        R_jj :: Matrix{T}
    ) :: Matrix{T} where {T <: AbstractFloat, I <: Integer}

    nq = size(q_jj, 1)
    nat_sc = size(R_jj, 1)
    nat = nat_sc ÷ nq
    n_random = size(v_q, 1)


    # Prepare the input and form the output
    R_lat = zeros(T, (3, nat_sc))
    q = zeros(T, (3, nq))

    R_lat .= R_jj'
    q .= q_jj'
    v_sc = zeros(T, (n_random, 3*nat_sc))

    println("Timing vector_q2r")
    vector_q2r!(v_sc, v_q, q, itau, R_lat)
    return v_sc
end


@doc raw"""
    multiply_matrix_vector_fourier!(
        v_tilde :: Array{Complex{T}, 3},
        Yq :: Array{Complex{T}, 3},
        u_tilde :: Array{Complex{T}, 3}
    ) where {T <: AbstractFloat}

Multiply the Fourier transform of the Y matrix with the Fourier transform of the u vector.

## Parameters

- v_tilde :: Array{Complex{T}, 3}
    The output vector  
- Yq :: Array{Complex{T}, 3}
    The Fourier transform of the Y matrix
- u_tilde :: Array{Complex{T}, 3}
    The Fourier transform of the u vector

## Returns

- v_tilde :: Array{Complex{T}, 3}
    The output vector
"""
function multiply_matrix_vector_fourier!(
    v_tilde:: Array{Complex{T}, 3},
    Yq :: Array{Complex{T}, 3},
    u_tilde:: Array{Complex{T}, 3}
) where {T <: AbstractFloat}

    nq = size(v_tilde, 3)
    n_random = size(v_tilde, 1)
    nat3 = size(v_tilde, 2)

    v_tilde .= 0

    for jq in 1:nq
        for k in 1:nat3
            for h in 1:nat3
                for i in 1:n_random
                    v_tilde[i, k, jq] += Yq[k, h, jq] *  u_tilde[i, h, jq]
                end
            end
        end
    end
end
function multiply_matrix_vector_fourier(
    Yq :: Array{Complex{T}, 3},
    u_tilde:: Array{Complex{T}, 3}
) :: Array{Complex{T}, 3} where {T <: AbstractFloat}
    n_random = size(u_tilde, 1)

    v_tilde = zeros(Complex{T}, size(u_tilde))
    multiply_matrix_vector_fourier!(v_tilde, Yq, u_tilde)
    return v_tilde
end


# @doc raw"""
#     get_uYu_fourier(
#         u_tilde :: Array{Complex{T}, 3},
#         Yq :: Array{Complex{T}, 3}
#     ) :: Vector{T} where {T <: AbstractFloat}
# 
# Compute the factor 
# 
# $$
# \left< u | Y | u \right>
# $$
# 
# ## Parameters
# 
# - u_tilde :: Array{Complex{T}, 3} (3*nat, nq, n_random)
#     The Fourier transform of the u vector
# - Yq :: Array{Complex{T}, 3} (3*nat, 3*nat, nq)
#     The Fourier transform of the Y matrix
# 
# ## Returns
# 
# - uYu :: Array{T, 1} (n_random)
#     The factor
# """
# function get_uYu_fourier(
#     u_tilde :: Array{Complex{T}, 3},
#     Yq :: Array{Complex{T}, 3}
# ) :: Vector{T} where {T <: AbstractFloat}
#     v_tilde = multiply_matrix_vector_fourier(Yq, u_tilde)
# 
#     n_random = size(u_tilde, 3)
#     nq = size(u_tilde, 2)
#     uYu = zeros(T, n_random)
# 
#     for i in 1:n_random
#         for jq in 1:nq
#             @views uYu[i] += real(u_tilde[:, jq, i]' * v_tilde[:, jq, i])
#         end
#     end
#     return uYu
# end
#     
    
function multiply_vector_vector_fourier!(
        result :: Vector{T},
        vector1 :: Array{Complex{T}, 3},
        vector2 :: Array{Complex{T}, 3}
        ) where {T <: AbstractFloat} 

    nq = size(vector1, 3)
    n_random = size(vector1, 1)
    nat3 = size(vector1, 2)
    
    result .= 0

    for jq in 1:nq
        for k in 1:nat3
            for i in 1:n_random
                result[i] += real(conj(vector1[i, k, jq]) * vector2[i, k, jq])
            end
        end
    end
end
function multiply_vector_vector_fourier(
        vector1 :: Array{Complex{T}, 3},
        vector2 :: Array{Complex{T}, 3}
        ) :: Vector{T} where {T <: AbstractFloat} 

    n_random = size(vector1, 1)
    
    result = zeros(T, n_random)
    multiply_vector_vector_fourier!(result, vector1, vector2)
    return result
end
