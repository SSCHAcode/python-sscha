using LinearAlgebra

@doc raw"""
    function get_gradient_fourier!(
        Φ_grad :: Array{U, 3},
        Φ_grad_err :: Array{U, 3},
        u_sc :: Matrix{T},
        δf_sc :: Matrix{T},
        R_lat :: Matrix{T},
        q :: Matrix{T},
        itau :: Array{Int},
        Y_matrix :: Matrix{T},
        q_opposite_index :: Vector{I})  {where T <: AbstractFloat, U <: AbstractComplex}

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
- u_sc : size (n_random, 3*nat_sc); the atomic displacements (in Bohr)
- δf_sc : size (n_random, 3*nat_sc); the forces (Ry / Bohr)
- R_lat : size (nat_sc, 3) the lattice vectors
- q : size (nq, 3) the q vectors
- itau : size (nat_sc); the index in the primitive cell corresponding to the supercell
- weights : size (n_random); the weights of the configurations
- Y_matrix : (3*nat_sc, 3*nat_sc) The inverse covariance matrix in the supercell.
- bg : (3, 3) The brilluin zone vectors (used to identify the q -> -q + G_)


$$
    Y_{ab} = \sqrt{m_a m_b} \sum_{\mu} e_{\mu}^a e_{\mu}^b \frac{2\omega_\mu}{\hbar (2 n_\mu + 1)}
$$

where $\omega_\mu$ and $e_\mu$ are the frequencies and eigenmodes.

The method returns the gredient (not symmetrized) in fourier space
"""
function get_gradient_fourier!(Φ_grad :: Array{Complex{T}, 3}, 
        Φ_grad_err :: Array{T, 3},
        u_sc :: Matrix{T},
        δf_sc :: Matrix{T},
        R_lat :: Matrix{T},
        q :: Matrix{T},
        itau :: Vector{I},
        weights :: Vector{T},
        Y :: Matrix{T},
        minus_q_index :: Vector{N}) where {T <: AbstractFloat, I <: Integer, N <: Integer}

    @time begin
        # Get the system size
        nat_sc = size(u_sc, 2) ÷ 3
        n_random = size(u_sc, 1)
        nq = size(q, 2)
        nat = nat_sc ÷ nq

        # Evaluate the v vector (product between u and Y matrix).
        v_vectors = zeros(T, (n_random, 3*nat_sc))
        mul!(v_vectors, u_sc, Y)


        # Fourier transform the v and δf vectors
        v_tilde = zeros(Complex{T}, (3*nat, nq, n_random))
        δf_tilde = zeros(Complex{T}, (3*nat, nq, n_random))
        q_dot_R = 0


        for jq ∈ 1:nq
            for k ∈ 1:nat_sc
                @views q_dot_R = q[:, jq]' * R_lat[:, k]
                exp_value = exp(- 1im * 2π * q_dot_R)

                for α in 1:3
                    index_sc = 3 * (k - 1) + α
                    index_uc = 3 * (itau[k] - 1) + α
                    for i ∈ 1:n_random
                        v_tilde[index_uc, jq, i] += exp_value * v_vectors[i, index_sc]
                        δf_tilde[index_uc, jq, i] += exp_value * δf_sc[i, index_sc]
                    end
                end
            end
        end

        v_tilde ./= √(nq)
        δf_tilde ./= √(nq)
    
    end

    # Compute the gradient exploiting a new average
    Φ_grad .= 0
    Φ_grad_err .= 0
    # Define helper arrays that avoid memory allocations during the for loop
    tmp = zeros(Complex{T}, (nat*3, nat*3))
    tmp2 = zeros(T, (nat*3, nat*3))

    @time begin 
        for i in 1:n_random
            for jq in 1:nq
                # The ' also perform the conjugation
                @views mul!(tmp, v_tilde[:, jq, i], δf_tilde[:, jq, i]')  
                @. tmp2 = tmp * conj(tmp)

                tmp .*= weights[i]
                tmp2 .*= weights[i]

                Φ_grad[:, :, jq] .+= tmp
                Φ_grad_err[:, :, jq] .+= tmp2
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
    u_sc :: Matrix{T},
    δf_sc :: Matrix{T},
    R_lat :: Matrix{T},
    q :: Matrix{T},
    itau :: Vector{I},
    weights :: Vector{T},
    Y :: Matrix{T},
    minus_q_vector :: Vector{N}) where {T <: AbstractFloat, I <: Integer, N <: Integer}
    
    nq = size(q, 1)
    nat = size(u_sc, 2) ÷ (3*nq) 
    nat_sc = nat * nq

    Φ_grad = zeros(Complex{T}, (nq, 3*nat, 3*nat))
    Φ_grad_err = zeros(T, (nq, 3*nat, 3*nat))
    Φ_grad_jj = zeros(Complex{T}, (3*nat, 3*nat, nq))
    Φ_grad_err_jj = zeros(T, (3*nat, 3*nat, nq))
    R_lat_jj = zeros(T, (3, nat_sc))
    q_jj = zeros(T, (3, nq))

    R_lat_jj .= R_lat'
    q_jj .= q'

    get_gradient_fourier!(Φ_grad_jj, Φ_grad_err_jj,
        u_sc, δf_sc, R_lat_jj, q_jj, itau, weights, Y, minus_q_vector)

    # Reorder the gradient
    for iq in 1:nq
        Φ_grad[iq, :, :] = Φ_grad_jj[:, :, iq]
        Φ_grad_err[iq, :, :] = Φ_grad_err_jj[:, :, iq]
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