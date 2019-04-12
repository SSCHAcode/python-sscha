from __future__ import print_function

# Common library imports
import numpy as np

# Cell constructor import
import cellconstructor as CC 
import cellconstructor.Phonons 
import cellconstructor.Methods

# Import the sscha modules
import Ensemble

"""
This source contains the subroutines to compute the dynamical spectrum.
They can be used for analysis of dynamical processes.
"""

# TO BE TESTED!!!!
def get_self_consistent_phonons(ensemble, include_v4 = False, compute_lifetimes = False, 
    smearing = 5e-5, n_iterations = 1, n_grid = 10000, w_max_factor = 2): #
    r"""
    GET THE DYNAMICAL PHONONS
    =========================

    This subroutines finds the solution of the equation:

    .. math::

        \omega^2 - D^{(s)} - \Re\Pi(q, \omega) = 0

    Where :math:`\Pi` is the dynamical sscha self-energy, and :math:`D^{(s)}` is the sscha dynamical matrix.
    
    The solution is found by computing the self energy on a coarse grid of w value and it is then interpolated.
    A new value of :math:`Pi` on the solution of the equation is then computed, to get the imaginary part to compute the lifetime.

    Parameters
    ----------
        ensemble : Ensemble()
            The ensemble to use to compute the anharmonic phonons
        include_v4 : bool
            If true the v4 is used in the propagator.
        smearing : float, optional
            The with of the full width half-maximum.
        n_iterations : int, optional
            If specified, the procedure is iterated for the specified number of times.
        n_grid : int, optional
            The number of point to perform the dynmat interpolation in between
        w_max_factor : int, optional
            The maximum searched frequency (max(w) * w_max_factor). By default is 2
    """

    # Get the first attempt with w = 0
    dyn_with_odd = ensemble.get_odd_correction(include_v4, frequency = 0, smearing = smearing)

    super_dyn = ensemble.current_dyn.GenerateSupercellDyn(ensemble.supercell)
    fake_dyn = super_dyn.Copy()
    fake_dyn.dynmats[0] += np.real(dyn_with_odd)
    fake_dyn.Symmetrize()

    freqs, dumb = fake_dyn.DyagDinQ(0)
    tot_freqs = [0]
    tot_dyns = [fake_dyn.dynmats[0]]

    # Delete degenerate mode
    final_freqs = freqs.copy()
    freqs = DeleteReplica(freqs)

    w_max = max(freqs) * w_max_factor
    w_values = np.linspace(0, w_max, n_grid)

    for i in range(n_iterations):
        # Get the fake_dyn for each frequency
        for w in freqs:
            dyn_with_odd = ensemble.get_odd_correction(include_v4, frequency = w, smearing = smearing)
            fake_dyn.dynmats[0] = np.real(dyn_with_odd) + super_dyn.dynmats[0]
            fake_dyn.Symmetrize()
            tot_freqs.append(w)
            tot_dyns.append(np.real(dyn_with_odd))
        
        # Order the frequencies and the matrices
        sort_mask = np.argsort(tot_freqs)
        tot_freqs = [tot_freqs[sort_mask[x]] for x in range(len(tot_freqs))]
        tot_dyns = [tot_dyns[sort_mask[x]] for x in range(len(tot_freqs))]

        # Solve the equation
        def get_ws(w):
            index_0 = np.sum( (np.array(tot_freqs) < w).astype(int)) - 1
            index_1 = index_0 + 1
            if len(tot_freqs) == index_1:
                index_1 = index_0
                index_0 -= 1
            
            # Interpolate the dynamical matrix
            x0 = tot_freqs[index_0]
            x1 = tot_freqs[index_1]
            x = (w - x0) / (x1-x0)
            dynmat = tot_dyns[index_0] + x * (tot_dyns[index_1] - tot_dyns[index_0])
            
            # Compute the frequencies
            fake_dyn.dynmats[0] = dynmat
            #fake_dyn.Symmetrize()
            ws, dumbs = fake_dyn.DyagDinQ(0)
            return ws

        # Find solutions
        new_freqs = []
        for i, w_value in enumerate(w_values):
            f_new = get_ws(w_value) - w_value

            if i == 0:
                f_old = f_new.copy
                w_old = w_value
                continue
            
            # Check if some frequency passed the check
            sol = np.min(f_old*f_new)
            if sol < 0:
                new_freqs.append(.5*w_value + .5*w_old)
                

            f_old = f_new.copy()
            w_old = w_value

        # Delete double frequencies
        freqs = DeleteReplica(new_freqs)
    
    return freqs


def DeleteReplica(array, threshold = 1e-6):
    """
    This simple function deletes the replica from an array if they are closer then
    the given threshold. The final array will be sorted.
    """


    new_array = np.sort(array)
    ret = []
    merge_set = []
    for i in range(1, len(new_array)):
        x_new = new_array[i]
        x_old = new_array[i-1]

        if x_new - x_old < threshold:
            merge_set.append(x_old)
        elif len(merge_set) != 0:
            merge_set.append(x_new)
            ret.append(np.mean(merge_set))
            merge_set = []

    return np.array(ret)


def get_spectral_function(dyn, supercell, self_energy, w_array):
    r"""
    COMPUTE THE SPECTRAL FUNCTION
    =============================

    This method computes the spectral function.

    .. math ::

        A(\omega) = - Tr\left[\Im G\right]

        G(\omega) = \left[\omega^2 - D^{(s)} - \Pi(\omega)\right]^{-1}
    
    where :math:`D^{(s)}` is the dynamical matrix, while :math:`\Pi` is the self-energy.

    Parameters
    ----------
        dyn : Phonons()
            The sscha dynamical matrix
        supercell: list of int
            The supercell of the calculation
        self_energy: list
            The list for each w in w_array of the self-energy matrices
        w_array : ndarray
            The array that contains the frequencies at which the self-energy has been
            computed. [Ry]
    
    Results
    -------
        A : ndarray(size = len(w_array))
            The spectral function.
    """
    # Get the dynamical matrix in the supercell
    superdyn = dyn.GenerateSupercellDyn(supercell)
    nat_sc = superdyn.structure.N_atoms

    # Convert the force constants into a dynamical matrix
    m = superdyn.structure.get_masses_array()
    m = np.tile(m, (3,1)).T.ravel()
    m_mat = np.sqrt(np.einsum("a,b", m, m)) #|m><m|
    dyn_mat = superdyn.dynmats[0] / m_mat
    I = np.eye(3*nat_sc, dtype = np.complex128)

    A = np.zeros(len(w_array))
    for i, w in enumerate(w_array):
        sigma = self_energy[i] / m_mat
        G_inv = w**2*I - dyn_mat - sigma

        # Invert the green function
        G = np.linalg.inv(G_inv)

        # Get the spectral function
        A[i] = np.einsum("aa", -np.imag(G))
    
    return A

def GetRamanResponce(dyn, supercell, self_energies, w_array, in_pol = None, out_pol = None, repeat_times = 10):
    """
    GET RAMAN SIGNAL
    ================

    This subroutines uses the dynamical self-energy to extract the Raman signal.
    We use the raman_tensor of the dyn variable to get the intensity of the scattered light.

    NOTE: The code will rise a ValueException if no raman tensor is found in the dynamical matrix.

    Parameters
    ----------
        dyn : Phonons()
            The sscha dynamical matrix (it must contain the raman tensor).
        supercell : list of 3 ints
            The size of the supercell (must match the q points)
        self_energies: list of (3xnatsc, 3xnatsc) 
            The dynamical w dependent part of the self-energy at many w values.
        w_array : ndarray
            The array of frequencies (in Ry) at which the self energy is computed.
        in_pol : 3d vector or None
            The polarization of the incoming light. If None a random one is chosen
        out_pol : 3d vector or None
            The polarization of the outcoming light. If None a random one is chosen.
        repeat_times : int
            If the polarizations are not provided, the light is supposed unpolarized
            and averaged on many directions. This is the number of averages to be carried
            out on many random polarization vectors

    Results
    -------
        raman_intensity : ndarray(shape = w_array.shape())
            The vector containing the Raman intensity for each frequency in the
            w_array.  
    """

    # Check if the given dynamical matrix has a Raman tensor
    if dyn.raman_tensor is None:
        raise ValueError("Error, the provided dyn must have a defined Raman tensor")

    # Check if the polarization vectors are defined
    repeat = False
    if in_pol is None:
        in_pol = np.random.normal(size=3)
        in_pol /= np.sqrt(in_pol.dot(in_pol))
        repeat = True
    if out_pol is None:
        out_pol = np.random.normal(size=3)
        out_pol /= np.sqrt(out_pol.dot(out_pol))
        repeat = True

    if not repeat:
        repeat_times = 1

    # Get the raman vector(s)
    v_ramans = []
    for i in range(repeat_times):
        v_raman = dyn.GetRamanVector(in_pol, out_pol)
        v_ramans.append(v_raman)

    # Get the dynamical matrix in the supercell
    superdyn = dyn.GenerateSupercellDyn(supercell)
    nat = dyn.structure.N_atoms

    # Convert the force constants into a dynamical matrix
    m = superdyn.structure.get_masses_array()
    m = np.tile(m, (3,1)).T.ravel()
    m_mat = np.sqrt(np.einsum("a,b", m, m)) #|m><m|
    I = np.eye(3*nat, dtype = np.complex128)

    raman_intensity = np.zeros(len(w_array), dtype = np.float64)
    for i, w in enumerate(w_array):

        # Get only the gamma point
        q_selfenergy = CC.Phonons.GetDynQFromFCSupercell(self_energies[i], np.array(dyn.q_tot), 
            dyn.structure, superdyn.structure)
        
        # Add the sscha static diagrams and transform to a dynamical matrix
        gamma_selfenergy = q_selfenergy[0, :, :] + dyn.dynmats[0]
        gamma_selfenergy /= m_mat 

        # Compute the gamma green function
        G_inv = w**2*I - gamma_selfenergy
        G = np.linalg.inv(G_inv)
        G /= m_mat

        # trace the imaginary part of the green function through the raman vector
        # In case of unpolarized light, average on many possible polarizations
        for v_raman in v_ramans:
            raman_intensity[i] += np.einsum("a, ab, b", v_raman, -np.imag(G), v_raman) / repeat_times
    
    return raman_intensity


