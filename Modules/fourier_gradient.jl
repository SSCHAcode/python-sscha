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
        bg :: Matrix{T})  {where T <: AbstractFloat, U <: AbstractComplex}

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
        bg :: Matrix{T}) where {T <: AbstractFloat, I <: Integer}
    # TODO: Check that U and T are the complex of the other


    # Get the system size
    nat_sc = size(u_sc, 2) ÷ 3
    n_random = size(u_sc, 1)
    nq = size(q, 1)
    nat = nat_sc ÷ nq

    # Evaluate the v vector (product between u and Y matrix).
    v_vectors = zeros(T, (3*nat_sc, n_random))
    v_vectors = Y * u_sc'


    # Fourier transform the v and δf vectors
    v_tilde = zeros(Complex{T}, (n_random, nq, 3*nat))
    δf_tilde = zeros(Complex{T}, (n_random, nq, 3*nat))

    for i ∈  1:n_random
        for jq ∈ 1:nq
            for k ∈ 1:nat_sc
                for α in 1:3
                    index_sc = 3 * (k - 1) + α
                    index_uc = 3 * (itau[k] - 1) + α
                    v_tilde[i, jq, index_uc] += exp(- 1im * 2π * q[jq] * R_lat[k]) * v_vectors[index_sc, i]
                    δf_tilde[i, jq, index_sc] += exp(+ 1im * 2π * q[jq] * R_lat[k]) * δf_sc[i, index_sc]
                end
            end
        end
    end


    # Compute the gradient exploiting a new average
    Φ_grad .= 0
    Φ_grad_err .= 0
    # Define helper arrays that avoid memory allocations during the for loop
    tmp = zeros(Complex{T}, (nat*3, nat*3))
    tmp2 = zeros(Complex{T}, (nat*3, nat*3))



    println("Computing gradient:")
    println("Size comparison:")
    println("Φ : $(size(Φ_grad))")
    println("TMP: $(size(tmp))")
    println()
    for i in 1:n_random
        for jq in 1:nq
            @views tmp .= v_tilde[i, jq, :] * δf_tilde[i, jq, :]'  
            tmp2 .= tmp 
            tmp2 .*= conj(tmp)

            println("JQ: $jq ; tmp modulus = $(sum(abs.(tmp)))")
            @show tmp

            tmp .*= weights[i]
            tmp2 .*= weights[i]

            Φ_grad[jq, :, :] .+= tmp
            Φ_grad_err[jq, :, :] .+= real.(tmp2)
        end
    end

    # Apply the transpose symmetry
    minus_q_index = zeros(Int, nq)
    get_opposit_q!(minus_q_index, q, bg)

    for iq in 1:nq
        @views Φ_grad[iq, :, :] .+= conj(Φ_grad[minus_q_index[iq], :, :])'
    end
    Φ_grad .*= 0.5
end
function get_gradient_fourier(
    u_sc :: Matrix{T},
    δf_sc :: Matrix{T},
    R_lat :: Matrix{T},
    q :: Matrix{T},
    itau :: Vector{I},
    weights :: Vector{T},
    Y :: Matrix{T},
    bg :: Matrix{T}) where {T <: AbstractFloat, I <: Integer}
    
    nq = size(q, 1)
    nat = size(u_sc, 2) ÷ (3*nq) 

    Φ_grad = zeros(Complex{T}, (nq, 3*nat, 3*nat))
    Φ_grad_err = zeros(T, (nq, 3*nat, 3*nat))

    get_gradient_fourier!(Φ_grad, Φ_grad_err,
        u_sc, δf_sc, R_lat, q, itau, weights, Y, bg)

    return Φ_grad, Φ_grad_err
end


@doc raw"""

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
function get_opposit_q!(
    opposite_index :: Vector{Int},
    q_list :: Matrix{T},
    bg :: Matrix{T}; far :: Int = 3) where {T <: AbstractFloat}

    opposite_index .= -1
    nq = size(q_list, 1)
    q_vector = zeros(T, 3)
    minus_q = zeros(T, 3)
    for i in 1:nq
        @views minus_q .= -q_list[i, :] 
        for j in 1:nq
            skip = false
            for xx in -far : far + 1
                if skip 
                    break
                end
                for yy in -far : far + 1
                    if skip 
                        break
                    end
                    for zz in -far : far + 1
                        begin 
                            @views q_vector .= q_list[i, :] + 
                                xx * bg[1, :] +
                                yy * bg[2, :] + 
                                zz * bg[3, :]
                        end
                        if max(abs.(minus_q .- q_vector)...) < 1e-8
                            opposite_index[i] = j
                            skip = true
                            break
                        end
                    end
                end
            end
        end
    end
    println()
    println("q tot:")
    println(q_list)

    println("Opposite q list:")
    println(opposite_index)
end