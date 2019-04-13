from __future__ import print_function

"""
This module performs the Lanczos algorithm in order to compute the responce function
of a particular perturbation.
"""

import numpy as np

# Import the scipy Lanczos modules
import scipy, scipy.sparse.linalg

import cellconstructor as CC
import cellconstructor.Phonons
import cellconstructor.symmetries



# Define a generic type for the double precision.
TYPE_DP = np.double

try:
    from ase.units import create_units
    units = create_units("2006")#Rydberg, Bohr
    Rydberg = units["Ry"]
    Bohr = units["Bohr"]
    __RyToK__ =  Rydberg / units["kB"]
    
except:
    Rydberg = 13.605698066
    Bohr = 1.889725989
    __RyToK__ = 157887.32400374097

def f_ups(w, T):
    """
    The eigenvalue of the upsilon matrix as a function of the frequency and the
    temperature
    """
    n_w = 0
    if T > 0:
        n_w = 1 / (np.exp(w * __RyToK__ / T) - 1)
    return 2*w / (1 + n_w)




class Lanczos:
    def __init__(self, ensemble):
        """
        INITIALIZE THE LANCZOS
        ======================

        This function extracts the weights, the X and Y arrays for the d3 and d4
        computation as well as the polarization vectors and frequencies.

        Parameters
        ----------
            ensemble : Ensemble.Ensemble()
                The ensemble upon which you want to compute the DynamicalResponce
        """

        self.dyn = ensemble.current_dyn.Copy() 
        superdyn = self.dyn.GenerateSupercellDyn(ensemble.supercell)
        self.uci_structure = ensemble.current_dyn.structure.copy()
        self.super_structure = superdyn.structure

        self.T = ensemble.current_T

        ws, pols = superdyn.DyagDinQ(0)

        self.nat = superdyn.structure.N_atoms

        self.qe_sym = CC.symmetries.QE_Symmetry(self.dyn.structure)
        self.qe_sym.SetupQPoint()

        # Get the masses
        m = superdyn.structure.get_masses_array()
        self.m = np.tile(m, (3,1)).T.ravel()

        # Remove the translations
        trans_mask = CC.Methods.get_translations(pols, m)

        # Get the polarization vectors
        self.w = ws[~trans_mask]
        self.pols = pols[:, ~trans_mask]

        self.n_modes = len(self.w)

        # Ignore v3 or v4. You can set them for testing
        self.ignore_v3 = False
        self.ignore_v4 = False

        # Get the info about the ensemble
        self.N = ensemble.N 
        self.rho = ensemble.rho.copy() 
        self.N_eff = np.sum(self.rho)

        u = ensemble.u_disps
        f = ensemble.forces.reshape(self.N, 3 * self.nat).copy()
        f -= ensemble.sscha_forces.reshape(self.N, 3 * self.nat)

        self.X = np.zeros((self.N, self.n_modes), order = "C", dtype = TYPE_DP)
        self.X[:,:] = np.einsum("a,ia, ab->ib", np.sqrt(self.m), u, self.pols)

        self.Y = np.zeros((self.N, self.n_modes), order = "C", dtype = TYPE_DP)
        self.Y[:,:] = np.einsum("a,ia, ab->ib", 1/np.sqrt(self.m), f, self.pols)

        # Prepare the variable used for the working
        self.psi = np.zeros(self.n_modes + self.n_modes*self.n_modes, dtype = TYPE_DP)

        # Prepare the L as a linear operator
        self.L_linop = scipy.sparse.linalg.LinearOperator(shape = (len(self.psi), len(self.psi)), matvec = self.apply_full_L, dtype = TYPE_DP)

        # Prepare the solution of the Lanczos algorithm
        self.eigvals = None
        self.eigvects = None 

        # Store the basis and the coefficients of the Lanczos procedure
        # In the custom lanczos mode
        self.a_coeffs = [] #Coefficients on the diagonal
        self.b_coeffs = [] # Coefficients close to the diagonal
        self.krilov_basis = [] # The basis of the krilov subspace
        self.arnoldi_matrix = None # If requested, the full quasi-tridiagonal matrix


    def prepare_perturbation(self, vector):
        """
        This function prepares the calculation for the Green function

        <v| G|v>

        Where |v> is the vector passed as input. If you want to compute the
        raman, for istance, it can be the vector of the Raman intensities.

        Parameters
        ----------
            vector: ndarray( size = (3*natoms))
                The vector of the perturbation for the computation of the green function
        """

        self.psi = np.zeros(self.n_modes + self.n_modes*self.n_modes, dtype = TYPE_DP)

        # Convert the vector in the polarization space
        new_v = np.einsum("a, a, ab->b", np.sqrt(self.m), vector, self.pols)
        self.psi[:self.n_modes] = new_v

    def get_vector_dyn_from_psi(self):
        """
        This function returns a standard vector and the dynamical matrix in cartesian coordinates
        
        This can be used to symmetrize the vector.
        """


        vector = self.psi[:self.n_modes]
        dyn = self.psi[self.n_modes:].reshape((self.n_modes, self.n_modes))

        w_a = np.tile(self.w, (self.n_modes, 1))
        w_b = np.tile(self.w, (self.n_modes, 1)).T 

        dyn *= ( 2*(w_a + w_b) * np.sqrt(w_a*w_b*(w_a + w_b)))

        # Get back the real vectors
        real_v = np.einsum("a, b, ab->a", 1/np.sqrt(self.m), vector, self.pols)
        real_dyn = np.einsum("ab, ca, db-> cd", dyn, self.pols, self.pols)
        real_dyn *= np.outer(np.sqrt(self.m), np.sqrt(self.m))

        return real_v, real_dyn 
    
    def set_psi_from_vector_dyn(self, vector, dyn):
        """
        Set the psi vector from a given vector and a force constant matrix.
        Used to reset the psi after the symmetrization.
        """

        new_v = np.einsum("a, a, ab->b",  np.sqrt(self.m), vector, self.pols)
        
        new_dyn = dyn / np.outer(np.sqrt(self.m), np.sqrt(self.m))
        
        w_a = np.tile(self.w, (self.n_modes, 1))
        w_b = np.tile(self.w, (self.n_modes, 1)).T 

        new_dyn /= ( 2*(w_a + w_b) * np.sqrt(w_a*w_b*(w_a + w_b)))

        self.psi[:self.n_modes] = new_v
        self.psi[self.n_modes:] = new_dyn.ravel()

    def symmetrize_psi(self):
        """
        Symmetrize the psi vector.
        """
        
        # First of all, get the vector and the dyn
        vector, dyn = self.get_vector_dyn_from_psi()

        # Symmetrize the vector
        self.qe_sym.SetupQPoint()
        new_v = np.zeros( (self.nat, 3), dtype = np.float64, order = "F")
        new_v[:,:] = vector.reshape((self.nat, 3))
        self.qe_sym.SymmetrizeVector(new_v)
        vector = new_v.ravel()

        # Symmetrize the dynamical matrix
        dyn_q = CC.Phonons.GetDynQFromFCSupercell(dyn, np.array(self.dyn.q_tot), self.uci_structure, self.super_structure)
        self.qe_sym.SymmetrizeFCQ(dyn_q, self.dyn.q_stars, asr = "custom")
        dyn = CC.Phonons.GetSupercellFCFromDyn(dyn_q, np.array(self.dyn.q_tot), self.uci_structure, self.super_structure)

        # Push everithing back into the psi
        self.set_psi_from_vector_dyn(vector, dyn)

    def apply_L1(self):
        """
        APPLY THE L1
        ============

        This is the first part of the application, it involves only harmonic propagation.

        Results
        -------
            out_vect : ndarray(shape(self.psi))
                It returns the application of the harmonic part of the L matrix
        """

        out_vect = np.zeros(np.shape(self.psi), dtype = TYPE_DP)

        # Get the harmonic responce function
        out_vect[:self.n_modes] = self.psi[:self.n_modes] * self.w**2

        # Get the harmonic responce on the propagator
        w_a = np.tile(self.w, (self.n_modes, 1))
        w_b = np.tile(self.w, (self.n_modes, 1)).T 

        new_out = (w_a + w_b)**2
        out_vect[self.n_modes:] = new_out.ravel() * self.psi[self.n_modes:]

        return out_vect

    def apply_L2(self):
        """
        APPLY THE L2
        ============

        L2 is the part of the L operators that mixes the two spaces.
        It involves the phi3 matrix.
        """


        w_a = np.tile(self.w, (self.n_modes, 1)).ravel()
        w_b = np.tile(self.w, (self.n_modes, 1)).T.ravel()

        vector = self.psi[:self.n_modes]
        dyn = self.psi[self.n_modes:]
        new_dyn = -dyn * np.sqrt( (w_a + w_b)/(w_a*w_b)) / 2

        out_v = SlowApplyD3ToDyn(self.X, self.Y, self.rho, self.w, self.T, new_dyn)
        out_d = SlowApplyD3ToVector(self.X, self.Y, self.rho, self.w, self.T, vector)
        out_d *= -np.sqrt( (w_a + w_b)/(w_a*w_b)) / 2

        out_vect = np.zeros(np.shape(self.psi), dtype = TYPE_DP)
        
        out_vect[:self.n_modes] = out_v
        out_vect[self.n_modes:] = out_d
        return out_vect

    def apply_L3(self):
        """
        APPLY THE L3
        ============

        This is the last part of the L matrix, it puts in communication 
        the dyn part of psi with herselfs.
        """

        w_a = np.tile(self.w, (self.n_modes, 1)).ravel()
        w_b = np.tile(self.w, (self.n_modes, 1)).T.ravel()

        dyn = self.psi[self.n_modes:] * np.sqrt((w_a + w_b) / (w_a * w_b)) / 2

        out_dyn = SlowApplyD4ToDyn(self.X, self.Y, self.rho, self.w, self.T, dyn)

        out_dyn *= np.sqrt((w_a + w_b) / (w_a * w_b)) / 2

        output = np.zeros(np.shape(self.psi), dtype = TYPE_DP)
        output[self.n_modes:] = out_dyn

        return output


    def apply_full_L(self, target=None):
        """
        APPLY THE L 
        ===========

        This function applies the L operator to the specified target vector.
        The target vector is first copied into the local psi, and the computed.
        This function will overwrite the current psi with the specified
        target.

        Parameters
        ----------
            target : ndarray ( size = shape(self.psi)), optional
                The garget vector to which you want to apply the
                full L matrix
        """

        # Setup the target vector to the self.psi
        if not target is None:
            self.psi = target 

        # Apply the whole L step by step to self.psi
        output = self.apply_L1()

        if not self.ignore_v3:
            output += self.apply_L2()
        if not self.ignore_v4:
            output += self.apply_L3()

        # Now return the output
        return output

    
    def custom_lanczos(self, iterations, save_dir = ".", correct_arnoldi = True):
        """
        RUN LANCZOS PROCEDURE
        =====================

        This method perform the lanczos procedure to tridiagonalize the L operator.
        It will save the current status after each iteration, to allow the user to simply restart.

        Parameters
        ----------
            iterations : int
                The number of Lanczos iterations to be done. 
                This correspond (in principle) to an exact dyagonalization if the number of iterations
                matches the dimension of the space. In practice the result should converge pretty fast.
            save_dir : string
                Path to the directiory in which to store the status after each iteration.
                This is necessary if you want to restart a calculation or look for convergence after each
                step. Use None if you do not want to use any I/O.
            correct_arnoldi : bool
                If True the Gram-Schmidt procedure will be iterated for all the basis, in this way also
                some element outside the tridiagonal form will be non zero, but this prevents the loss of
                orthogonality, commonly faced in standard Lanczos algorithms.
        """

        # TODO: Implement the full Lanczos procedure
        pass

    def save_status(self, file):
        """
        Save the current data in npz compressed format, in order to reanalyze easily the result (or restart the Lanczos)
        later.

        Parameters
        ----------
            file : string
                Path to where you want to save the data. It must be an npz binary format. The extension
                will be added if it does not match the npz one
        """


        # Add the correct extension
        if not ".npz" in file.lower():
            file += ".npz"
        
        # Save all the data
        np.savez_compressed(file, 
                            T = self.T,
                            nat = self.nat,
                            m = self.m,
                            w = self.w,
                            pols = self.pols,
                            n_modes = self.n_modes,
                            ignore_v3 = self.ignore_v3,
                            ignore_v4 = self.ignore_v4,
                            N = self.N,
                            rho = self.rho
                            X = self.X,
                            Y = self.Y,
                            psi = self.psi,
                            a_coeffs = self.a_coeffs,
                            b_coeffs = self.b_coeffs,
                            krilov_basis = self.krilov_basis,
                            arnoldi_matrix = self.arnoldi_matrix)
            
        


    def run_full_diag(self, number, discard_dyn = True, n_iter = 100):
        r"""
        FULL LANCZOS DIAGONALIZATION
        ============================

        This function runs the standard Lanczos iteration progress.
        It returns the eigenvalues and eigenvectors of the L operator.
        These can be used for computing the spectral function, and the full
        green function as:

        .. math ::

            G_{ab}(\omega) = \sum_{\alpha} \frac{\left<a | \lambda_\alpha\right>\left<\lambda_\alpha|b\right>}{\lambda_\alpha - \omega^2 + i\eta}

        where :math:`\lambda` are eigenvalues and vectors returned by this method, while :math:`\eta` is a
        smearing parameter chosen by the user. 
        Remember the eigenvectors are defined in the polarization basis and they comprend also the dynamical matrix degrees of freedom.
        Since in most application you want to discard the dynamical matrices, you can select discard_din = True.

        The used Lanczos algorithm is the one by ARPACK, as implemented in scipy.sparse module

        Parameters
        ----------
            number = int
                The number of the n highest eigenvalues to be found
            discard_dyn : bool, optional
                If True the dynamical matrix component of the output eigenvectors will be discarded.
            n_iter : int, optional
                The maximum number of Lanczos iterations. Usually must be much higher than the
                number of states you want to describe.
        """

        # Perform the lanczos operation
        eigvals, eigvects = scipy.sparse.linalg.eigsh(self.L_linop, k = number, v0 = self.psi, ncv= n_iter)

        self.eigvals = eigvals
        self.eigvects = eigvects

        # Check if the dynamical part must be discarded
        if discard_dyn:
            eigvects = eigvects[:self.n_modes, :]
    

        return eigvals, eigvects

    def GetSupercellSpectralFunctionFromEig(self, w_array, smearing):
        r"""
        GET SPECTRAL FUNCTION
        =====================

        Get the spectral function from the eigenvalues and eigenvectors.
        The method run_full_diag must already be runned.

        This method returns the spectral function in the supercell.
        The spectral function is computed as:

        .. math ::

            G_{ab}(\omega) = \sum_{\alpha} \frac{\left<a | \lambda_\alpha\right>\left<\lambda_\alpha|b\right>}{\lambda_\alpha - \omega^2 + i\eta}

        where :math:`\lambda` are eigenvalues and vectors returned by this method, while :math:`\eta` is a
        smearing parameter chosen by the user.

        Parameters
        ----------
            w_array : ndarray
                The frequencies to which you want to compute the spectral function.
            smearing : float
                The smearing of the spectral function.

        Returns
        -------
            s(w) : ndarray
                The -ImG(w), the opposite of the imaginary part of the Green function. 
        """

        # Exclude dynamical
        eigvects = self.eigvects[:self.n_modes, :]

        N_w = len(w_array)
        N_alpha = len(self.eigvals)

        # Transform the vectors back in cartesian coordinates
        new_vects = np.einsum("ab, ca, c->cb", eigvects, self.pols, 1 / np.sqrt(self.m))

        spectral_weight = np.einsum("ab, ab -> b", new_vects, np.conj(new_vects))
        spectral_function = np.zeros(N_w, dtype = np.complex128)

        l_alpha = np.tile(self.eigvals, (N_w, 1))
        p_w = np.tile(spectral_weight, (N_w, 1))
        _w_ = np.tile(w_array, (N_alpha, 1)).T 

        big_mat = p_w / (l_alpha - _w_**2 + 1j*smearing)
        spectral_function[:] = np.sum(big_mat, axis = 1)

        return - np.imag(spectral_function)


    def GetFullSelfEnergy(self):
        r"""
        GET SELF ENERGY 
        ===============

        Get the self-energy matrix from the eigenvalues and eigenvectors.
        The method run_full_diag must already be runned.

        This method returns the self energy in the supercell.
        It is computed as

        .. math ::

            \Pi_{ab} = \sum_{\alpha} \lambda_\alpha\left<a | \lambda_\alpha\right>\left<\lambda_\alpha|b\right>

        where :math:`\lambda` are eigenvalues and vectors returned by this method.
        The matrix is in real (cartesian) space.

        Returns
        -------
            s(w) : ndarray
                The -ImG(w), the opposite of the imaginary part of the Green function. 
        """

        # Exclude dynamical
        eigvects = self.eigvects[:self.n_modes, :]


        # Transform the vectors back in cartesian coordinates
        new_vects = np.einsum("ab, ca, c->cb", eigvects, self.pols, 1 / np.sqrt(self.m))

        self_energy = np.einsum("ab, cb, b", new_vects, np.conj(new_vects), self.eigvals)

        return self_energy




def SlowApplyD3ToDyn(X, Y, rho, w, T, input_dyn):
    """
    Apply the D3 vector.

    This is a testing function. It is slow, as it is a pure python implementation.
    """

    new_X = np.einsum("ab,b->ab", X, f_ups(w, T))

    
    n_rand, n_modes = np.shape(X)
    N_eff = np.sum(rho)

    v_out = np.zeros(n_modes, dtype = TYPE_DP)
    for a in range(n_modes):
        for b in range(n_modes):
            for c in range(n_modes):
                # Prepare the D3 calculation
                in_av = new_X[:, a] * new_X[:, b] * Y[:, c]
                in_av +=  new_X[:, a] * new_X[:, c] * Y[:, b]
                in_av +=  new_X[:, c] * new_X[:, b] * Y[:, a]
                in_av *= rho

                # Apply D3
                v_out[a] += - np.sum(in_av) * input_dyn[n_modes*b + c] / (3*N_eff)
    
    return v_out


def SlowApplyD3ToVector(X, Y, rho, w, T, input_vector):
    """
    Apply the D3 vector.

    This is a testing function. It is slow, as it is a pure python implementation.
    """

    new_X = np.einsum("ab,b->ab", X, f_ups(w, T))
    
    n_rand, n_modes = np.shape(X)
    N_eff = np.sum(rho)

    v_out = np.zeros(n_modes*n_modes, dtype = TYPE_DP)
    for a in range(n_modes):
        for b in range(n_modes):
            for c in range(n_modes):
                # Prepare the D3 calculation
                in_av = new_X[:, a] * new_X[:, b] * Y[:, c]
                in_av +=  new_X[:, a] * new_X[:, c] * Y[:, b]
                in_av +=  new_X[:, c] * new_X[:, b] * Y[:, a]
                in_av *= rho

                # Apply D3
                v_out[a*n_modes + b] += - np.sum(in_av) * input_vector[c] / (3*N_eff)
    
    return v_out




def SlowApplyD4ToDyn(X, Y, rho, w, T, input_dyn):
    """
    Apply the D4 matrix.

    This is a testing function. It is slow, as it is a pure python implementation.
    """


    new_X = np.einsum("ab,b->ab", X, f_ups(w, T))

    
    n_rand, n_modes = np.shape(X)
    N_eff = np.sum(rho)

    v_out = np.zeros(n_modes*n_modes, dtype = TYPE_DP)
    for a in range(n_modes):
        for b in range(n_modes):
            for c in range(n_modes):
                for d in range(n_modes):
                    # Prepare the D3 calculation
                    in_av =  new_X[:, a] * new_X[:, b] * new_X[:, c] * Y[:, d]
                    in_av += new_X[:, a] * new_X[:, b] * Y[:, c] * new_X[:, d]
                    in_av += new_X[:, a] * Y[:, b] * new_X[:, c] * new_X[:, d]
                    in_av += Y[:, a] * new_X[:, b] * new_X[:, c] * new_X[:, d]

                    in_av *= rho

                    # Apply D3
                    v_out[a*n_modes + b] += - np.sum(in_av) * input_dyn[n_modes*c + d] / (4*N_eff)
    
    return v_out
