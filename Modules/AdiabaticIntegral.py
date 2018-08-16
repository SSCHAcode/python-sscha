import Ensemble
"""
This file contains the method to get the exact Free energy derivatives.
"""


def GetExactFreeEnergy(ensemble, N_int):
    """
    GET EXACT FREE ENERGY
    =====================

    This function computes the exact free energy of the system performing an adiabatic
    integral over the ensemble:

    .. math::
    
        F = F_0 + \\int_0^1 \\frac{dF}{d\\lambda} d\\lambda

        H_\\lambda = H_0 + (H - H_0)\\lambda

    Where :math:`\\lambda` is an adiabatic parameter, :math:`H` is the real Hamiltonian of the
    system and :math:`H_0` is the hamiltonian used to generate the given ensemble.

    The derivative of the free energy is:

    .. math::
 
        \\frac{dF}{d\\lambda}(\\lambda) = \\left< \\frac{dH_\\lambda}{d\\lambda}\\right>_{\\rho_\\lambda}

    Where :math:`\\rho_\\lambda` is the density matrix of the :math:`H_\\lambda` Hamiltonian.

    This can be trivially obtained as follows


    .. math::

        \\left< \\frac{dH_\\lambda}{d\\lambda}\\right>_{\\rho_\\lambda} = \\sum_{i = 1}^{N_c} \\left[V(r_i) - V_0(r_i)\\right] \\frac{ e^{-\\beta \\left(H_\\lambda(r_i) - H_0\\right)}}{Z_\\lambda}

        Z_\\lambda = \\sum_{i = 1}^{N_c}  e^{-\\beta \\left(H_\\lambda(r_i) - H_0\\right)}

    
    Parameters
    ----------
        ensemble : Ensemble()
            The ensemble on which to perform the calculation
        N_int : int
            The number of step on which to divide the integration.

    Returns
    -------
        free energy : double
        error_on_free_energy : double

    """
    K_to_Ry=6.336857346553283e-06
    
    # Prepare the important sampling
    rho = np.zeros((N_int, ensemble.N))

    Lambda = np.linspace(0, 1, N_int)
    
    on_exponent = np.outer(Lambda, ensemble.energies - ensemble.sscha_energies)

    # TODO: TO BE CONTINUE

    
    
    
