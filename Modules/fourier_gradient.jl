using LinearAlgebra

@doc raw"""
    function get_gradient_fourier(
        Φ_grad :: Array{U, 3},
        Φ_grad_err :: Array{U, 3},
        u_sc :: Matrix{T},
        δf_sc :: Matrix{T},
        R_lat :: Matrix{T},
        q :: Matrix{T},
        itau :: Array{Int},
        Y_matrix :: Matrix{T})  {where T <: AbstractFloat, U <: AbstractComplex}
 

Compute the gradient in fourier space.

.. math:: 

    G_\Phi(q) = < Yu(q) δf(-q)>


where

.. math::

    u_a(q) = \sum_R e^{-i q ⋅ R} u_a(R)

All the average are reported by summing times the weights
without normalization. 
This makes it easier to enable parallelization in the error 
calculation and compute the actual average in a super function. 

Parameters 
----------

- Φ_grad : size (nq, 3*nat, 3*nat); the gradient of the free energy 
- Φ_grad_err : size (nq, 3*nat, 3*nat); the average of squares. Used to evaluate the variance 
- u_sc : size (n_random, 3*nat_sc); the atomic displacements (in Bohr)
- δf_sc : size (n_random, 3*nat_sc); the forces (Ry / Bohr)
- R_lat : size (nat_sc, 3) the lattice vectors
- q : size (nq, 3) the q vectors
- itau : size (nat_sc); the index in the primitive cell corresponding to the supercell
- weights : size (n_random); the weights of the configurations
- Y_matrix : (3*nat_sc, 3*nat_sc) The inverse covariance matrix in the supercell.


.. math:: 

    Y_{ab} = \sqrt{m_a m_b} \sum_{\mu} e_{\mu}^a e_{\mu}^b \frac{2\omega_\mu}{\hbar (2 n_\mu + 1)}


where :math:`\omega_\mu` and :math:`e_\mu` are the frequencies and eigenmodes.

The method returns the gredient (not symmetrized) in fourier space

Results
-------

"""
function get_gradient_fourier!(Φ_grad :: Array{U, 3}, 
        u_sc :: Matrix{T},
        δf_sc :: Matrix{T},
        R_lat :: Matrix{T},
        q :: Matrix{T},
        itau :: Array{Int},
        weights :: Array{T},
        Y_matrix :: Matrix{T}) {where T <: AbstractFloat, U <: AbstractComplex}
    # TODO: Check that U and T are the complex of the other


    # Get the system size
    nat_sc = size(u_sc), 2) ÷ 3
    n_random = size(u_sc, 1)
    nq = size(q, 1)
    nat = nat_sc ÷ nq

    # Evaluate the v vector (product between u and Y matrix).
    v_vectors = zeros(T, (3*nat_sc, n_random))
    v_vectors = Y * u_sc'


    # Fourier transform the v and δf vectors
    v_tilde = zeros(U, (n_random, nq, 3*nat))
    δf_tilde = zeros(U, (n_random, nq, 3*nat))

    for i ∈  1:n_random
        for jq ∈ 1:nq
            for k ∈ 1:nat_sc
                v_tilde[i, jq, itau[k]] += exp(- 1im * 2π * q[jq] * R_lat[k]) * v_vectors[k, i]
                δf_tilde[i, jq, itau[k]] += exp(+ 1im * 2π * q[jq] * R_lat[k]) * δf_sc[k, i]
            end
        end
    end


    # Compute the gradient exploiting a new average
    Φ_grad .= 0
    Φ_grad_err .= 0
    # Define helper arrays that avoid memory allocations during the for loop
    tmp = zeros(U, (nat*3, nat*3))
    tmp2 = zeros(U, (nat*3, nat*3))

    for i in 1:n_random
        for jq in 1:nq
            @views tmp .= v_tilde[i, jq, :] * δf_tilde[i, jq, :]'  
            @views tmp2 .= tmp .^2

            tmp .*= weights[i]
            tmp2 .*= weights[i]

            Φ_grad[jq, :, :] .+= tmp
            Φ_grad_error[jq, :, :] .+= tmp2
        end
    end
end

