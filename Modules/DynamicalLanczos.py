from __future__ import print_function
from __future__ import division

"""
This module performs the Lanczos algorithm in order to compute the responce function
of a particular perturbation.
"""

import sys, os
import time
import numpy as np

from timeit import default_timer as timer

# Import the scipy Lanczos modules
import scipy, scipy.sparse.linalg

import cellconstructor as CC
import cellconstructor.Phonons
import cellconstructor.symmetries

import Ensemble
import sscha_HP_odd

# Override the print function to print in parallel only from the master
import Parallel
from Parallel import pprint as print

# Define a generic type for the double precision.
TYPE_DP = np.double
__EPSILON__ = 1e-8
N_REP_ORTH = 1

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


__SPGLIB__ = True
try:
    import spglib
except:
    __SPGLIB__ = False
    

def f_ups(w, T):
    """
    The eigenvalue of the upsilon matrix as a function of the frequency and the
    temperature
    """
    n_w = 0
    if T > 0:
        n_w = 1 / (np.exp(w * __RyToK__ / T) - 1)
    return 2*w / (1 + n_w)

# Modes for the calculation
MODE_FAST_MPI = 2
MODE_FAST_SERIAL = 1
MODE_SLOW_SERIAL = 0


class Lanczos:
    def __init__(self, ensemble = None, mode = 1):
        """
        INITIALIZE THE LANCZOS
        ======================

        This function extracts the weights, the X and Y arrays for the d3 and d4
        computation as well as the polarization vectors and frequencies.

        Parameters
        ----------
            ensemble : Ensemble.Ensemble()
                The ensemble upon which you want to compute the DynamicalResponce
            mode : int
                The mode of the speedup.
                   0) Slow python implementation 
                      Use this just for testing
                   1) Fast C serial code
                   2) Fast C parallel (MPI)

        """

        self.mode = mode

        # Define the order
        order = "C"
        if self.mode >= 1:
            order = "F"
    
        self.T = 0
        self.nat = 0
        self.m = []
        self.w = []
        self.pols = []
        self.n_modes = 0
        self.ignore_harmonic = False 
        self.ignore_v3 = False
        self.ignore_v4 = False
        self.N = 0
        self.rho = []
        self.N_eff = 0
        self.X = []
        self.Y = []
        self.psi = []
        self.eigvals = None
        self.eigvects = None
        # In the custom lanczos mode
        self.a_coeffs = [] #Coefficients on the diagonal
        self.b_coeffs = [] # Coefficients close to the diagonal
        self.c_coeffs = [] # Coefficients in the case of the biconjugate Lanczos
        self.krilov_basis = [] # The basis of the krilov subspace
        self.basis_P = [] # The basis of the P vectors for the biconjugate Lanczos
        self.basis_Q = [] # The basis of the Q vectors for the biconjugate Lanczos
        self.arnoldi_matrix = [] # If requested, the upper triangular arnoldi matrix
        self.reverse_L = False
        self.shift_value = 0
        self.symmetrize = False
        self.symmetries = None
        self.degenerate_space = None
        self.N_degeneracy = None
        self.initialized = False
        self.perturbation_modulus = 1
        self.q_vectors = None # The q vectors of each mode

        # Perform a bare initialization if the ensemble is not provided
        if ensemble is None:
            return


        self.dyn = ensemble.current_dyn.Copy() 
        superdyn = self.dyn.GenerateSupercellDyn(ensemble.supercell)
        self.uci_structure = ensemble.current_dyn.structure.copy()
        self.super_structure = superdyn.structure

        self.T = ensemble.current_T

        ws, pols = self.dyn.DiagonalizeSupercell()


        self.nat = superdyn.structure.N_atoms
        n_cell = np.prod(self.dyn.GetSupercell())

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


        # Prepare the list of q point starting from the polarization vectors
        #q_list = CC.symmetries.GetQForEachMode(self.pols, self.uci_structure, self.super_structure, self.dyn.GetSupercell())
        # Store the q vectors in crystal space
        bg = self.uci_structure.get_reciprocal_vectors() / 2* np.pi
        self.q_vectors = np.zeros((self.n_modes, 3), dtype = np.double, order = "C")
        #for iq, q in enumerate(q_list):
        #    self.q_vectors[iq, :] = CC.Methods.covariant_coordinate(bg, q)
        


        # Ignore v3 or v4. You can set them for testing
        self.ignore_v3 = False
        self.ignore_v4 = False

        # Get the info about the ensemble
        self.N = ensemble.N 
        self.rho = ensemble.rho.copy() 
        self.N_eff = np.sum(self.rho)

        u = ensemble.u_disps / Ensemble.Bohr
        f = ensemble.forces.reshape(self.N, 3 * self.nat).copy()
        f -= ensemble.sscha_forces.reshape(self.N, 3 * self.nat)
        # Subtract also the average force to clean more the stochastic noise
        av_force = ensemble.get_average_forces(get_error = False).ravel()
        new_av_force = np.tile(av_force, (n_cell, 1)).ravel()
        f -= np.tile(new_av_force, (self.N, 1)) 
        f *= Ensemble.Bohr

        self.X = np.zeros((self.N, self.n_modes), order = order, dtype = TYPE_DP)
        #self.X[:,:] = np.einsum("a,ia, ab->ib", np.sqrt(self.m), u, self.pols) / Ensemble.Bohr

        self.Y = np.zeros((self.N, self.n_modes), order = order, dtype = TYPE_DP)
        #self.Y[:,:] = np.einsum("a,ia, ab->ib", 1/np.sqrt(self.m), f, self.pols) * Ensemble.Bohr

        ups_mat = self.dyn.GetUpsilonMatrix(self.T)
        v_vec = u.dot(ups_mat)

        # Convert in the polarization basis
        pol_mat = np.einsum("ab, a->ab", self.pols, 1 / np.sqrt(self.m))
        self.X[:,:] = v_vec.dot(pol_mat)
        self.Y[:,:] = f.dot(pol_mat)

        # Prepare the variable used for the working
        len_psi = self.n_modes
        if self.T < __EPSILON__:
            len_psi += self.n_modes**2
        else:
            len_psi += self.n_modes * (self.n_modes + 1)
            print("N MODES:", self.n_modes)
            print("LEN PSI:", len_psi)

        self.psi = np.zeros(len_psi, dtype = TYPE_DP)

        # Prepare the L as a linear operator (Prepare the possibility to transpose the matrix)
        def L_transp(psi):
            return self.apply_full_L(psi, transpose= True)
        self.L_linop = scipy.sparse.linalg.LinearOperator(shape = (len(self.psi), len(self.psi)), matvec = self.apply_full_L, rmatvec = L_transp, dtype = TYPE_DP)

        # Prepare the solution of the Lanczos algorithm
        self.eigvals = None
        self.eigvects = None 

        # Store the basis and the coefficients of the Lanczos procedure
        # In the custom lanczos mode
        self.a_coeffs = [] #Coefficients on the diagonal
        self.b_coeffs = [] # Coefficients close to the diagonal
        self.g_coeffs = [] # Coefficients for the non symmetric Lanczos (finite temperature)
        self.krilov_basis = [] # The basis of the krilov subspace
        self.arnoldi_matrix = [] # If requested, the upper triangular arnoldi matrix

        # These are some options that can be used to properly reverse and shift the L operator to
        # fasten the convergence of low energy modes
        self.reverse_L = False
        self.shift_value = 0
        self.symmetrize = False


    def reset(self):
        """
        RESET THE LANCZOS
        =================

        This function reset the Lanczos algorithm, allowing for a new responce function calculation
        with the same ensemble and the same settings.
        """


        # Prepare the solution of the Lanczos algorithm
        self.eigvals = None
        self.eigvects = None 

        # Store the basis and the coefficients of the Lanczos procedure
        # In the custom lanczos mode
        self.a_coeffs = [] #Coefficients on the diagonal
        self.b_coeffs = [] # Coefficients close to the diagonal
        self.c_coeffs = []

        # The krilov basis for the symmetric and unsymmetric Lanczos
        self.basis_P = []
        self.basis_Q = []
        self.krilov_basis = [] # The basis of the krilov subspace
        self.arnoldi_matrix = [] # If requested, the upper triangular arnoldi matrix


    def init(self):
        """
        INITIALIZE THE CALCULATION
        ==========================

        Perform everithing needed to initialize the calculation
        """
        self.prepare_symmetrization()
        self.initialized = True


    def prepare_symmetrization(self, no_sym = False):
        """
        PREPARE THE SYMMETRIZATION
        ==========================

        This function analyzes the character of the symmetry operations for each polarization vectors.
        This will allow the method do know how many modes are degenerate.

        Parameters
        ----------
            no_sym : bool
                If True, the symmetries are neglected.
        """

        self.initialized = True

        # Generate the dynamical matrix in the supercell
        super_structure = self.dyn.structure.generate_supercell(self.dyn.GetSupercell())
        w, pols = self.dyn.DiagonalizeSupercell()

        # Get the symmetries of the super structure
        if not __SPGLIB__ and not no_sym:
            raise ImportError("Error, spglib module required to perform symmetrization in a supercell. Otherwise, use no_sym")
        
        # Neglect the symmetries
        if no_sym:
            self.symmetries = np.zeros( (1, self.n_modes, self.n_modes), dtype = TYPE_DP)
            self.symmetries[0, :, :] = np.eye(self.n_modes)
            self.N_degeneracy = np.ones(self.n_modes, dtype = np.intc)
            self.degenerate_space = [[i] for i in range(self.n_modes)]
            return

        super_symmetries = CC.symmetries.GetSymmetriesFromSPGLIB(spglib.get_symmetry(super_structure.get_ase_atoms()), False)

        # Get the symmetry matrix in the polarization space
        # Translations are needed, as this method needs a complete basis.
        pol_symmetries = CC.symmetries.GetSymmetriesOnModes(super_symmetries, super_structure, pols)

        Ns, dumb, dump = np.shape(pol_symmetries)
        
        # Now we can pull out the translations
        trans_mask = CC.Methods.get_translations(pols, super_structure.get_masses_array())
        self.symmetries = np.zeros( (Ns, self.n_modes, self.n_modes), dtype = TYPE_DP)
        ptmp = pol_symmetries[:, :,  ~trans_mask] 
        self.symmetries[:,:,:] = ptmp[:, ~trans_mask, :]

        # Get the degeneracy
        w = w[~trans_mask]
        N_deg = np.ones(len(w), dtype = np.intc)
        start_deg = -1
        deg_space = [ [x] for x in range(self.n_modes)]
        for i in range(1, len(w)):
            if np.abs(w[i-1] - w[i]) < __EPSILON__ :
                N_deg[i] = N_deg[i-1] + 1

                if start_deg == -1:
                    start_deg = i - 1

                for j in range(start_deg, i):
                    N_deg[j] = N_deg[i]
                    deg_space[j].append(i)
                    deg_space[i].append(j)
            else:
                start_deg = -1


        self.N_degeneracy = N_deg
        self.degenerate_space = deg_space

    def prepare_raman(self, pol_vec_in= np.array([1,0,0]), pol_vec_out = np.array([1,0,0])):
        """
        PREPARE LANCZOS FOR RAMAN SPECTRUM
        ==================================

        This subroutines prepare the perturbation for the Raman signal.

        The raman tensor is read from the dynamical matrix provided by the original ensemble.

        Parameters
        ----------
            pol_vec_in : ndarray (size =3)
                The polarization vector of the incoming light
            pol_vec_out : ndarray (size = 3)
                The polarization vector for the outcoming light
        """

        # Check if the raman tensor is present
        assert not self.dyn.raman_tensor is None, "Error, no Raman tensor found. Cannot initialize the Raman responce"

        # Get the raman vector
        raman_v = self.dyn.GetRamanVector(pol_vec_in, pol_vec_out)

        # Get the raman vector in the supercell
        n_supercell = np.prod(self.dyn.GetSupercell())
        new_raman_v = np.tile(raman_v.ravel(), n_supercell)

        # Convert in the polarization basis and store the intensity
        self.prepare_perturbation(new_raman_v, masses_exp=-1)


    def prepare_ir(self, effective_charges = None, pol_vec = np.array([1,0,0])):
        """
        PREPARE LANCZOS FOR INFRARED SPECTRUM COMPUTATION
        =================================================

        In this subroutine we prepare the lanczos algorithm for the computation of the
        infrared spectrum signal.

        Parameters
        ----------
            effective_charges : ndarray(size = (n_atoms, 3, 3), dtype = np.double)
                The effective charges. Indices are: Number of atoms in the unit cell,
                electric field component, atomic coordinate. If None, the effective charges
                contained in the dynamical matrix will be considered.
            pol_vec : ndarray(size = 3)
                The polarization vector of the light.
        """

        ec = self.dyn.effective_charges
        if not effective_charges is None:
            ec = effective_charges
        

        n_supercell = np.prod(self.dyn.GetSupercell())

        # Check the effective charges
        assert not ec is None, "Error, no effective charge found. Cannot initialize IR responce"

        ec_size = np.shape(ec)
        MSG = """
        Error, effective charges of the wrong shape: {}
        """.format(ec_size)
        assert len(ec_size) == 3, MSG
        assert ec_size[0] * ec_size[2] * n_supercell == self.n_modes + 3
        assert ec_size[1] == ec_size[2] == 3

        z_eff = np.einsum("abc, b", ec, pol_vec)

        # Get the gamma effective charge
        new_zeff = np.tile(z_eff.ravel(), n_supercell)
        self.prepare_perturbation(new_zeff, masses_exp = -1)


    def prepare_perturbation(self, vector, masses_exp = 1):
        r"""
        This function prepares the calculation for the Green function

        <v| G|v>

        Where |v> is the vector passed as input. If you want to compute the
        raman, for istance, it can be the vector of the Raman intensities.

        The vector can be obtained contracting it with the polarization vectors.
        The contraction can be on the numerator or on the denumerator, depending on the
        observable.

        .. math ::

            v_\mu = \sum_a v_a e_\mu^a \cdot \sqrt{m_a}

            v_\mu = \sum_a v_a \frac{e_\mu^a}{  \sqrt{m_a}}

        Parameters
        ----------
            vector: ndarray( size = (3*natoms))
                The vector of the perturbation for the computation of the green function
            masses_exp : float
                The vector is multiplied by the square root of the masses raised to masses_exp.
                If you want to multiply each component by the square root of the masses use 1,
                if you want to divide by the quare root use -1, use 0 if you do not want to use the
                masses.
        """

        self.psi = np.zeros(self.psi.shape, dtype = TYPE_DP)

        # Convert the vector in the polarization space
        m_on = np.sqrt(self.m) ** masses_exp
        new_v = np.einsum("a, a, ab->b", m_on, vector, self.pols)
        self.psi[:self.n_modes] = new_v

        self.perturbation_modulus = new_v.dot(new_v)

        if self.symmetrize:
            self.symmetrize_psi()

    def get_vector_dyn_from_psi(self):
        """
        This function returns a standard vector and the dynamical matrix in cartesian coordinates
        
        This can be used to symmetrize the vector.

        The vector is returned in [Bohr] and the force constant matrix is returned in [Ry/bohr^2]
        """


        vector = self.psi[:self.n_modes]
        dyn = self.psi[self.n_modes:].reshape((self.n_modes, self.n_modes))

        w_a = np.tile(self.w, (self.n_modes, 1))
        w_b = np.tile(self.w, (self.n_modes, 1)).T 

        dyn *= ( 2*(w_a + w_b) * np.sqrt(w_a*w_b*(w_a + w_b)))

        # Get back the real vectors
        real_v = np.einsum("b, ab->a",  vector, self.pols)  /np.sqrt(self.m)
        real_dyn = np.einsum("ab, ca, db-> cd", dyn, self.pols, self.pols)
        real_dyn *= np.outer(np.sqrt(self.m), np.sqrt(self.m))

        return real_v, real_dyn 
    
    def set_psi_from_vector_dyn(self, vector, dyn):
        """
        Set the psi vector from a given vector of positions [bohr] and a force constant matrix [Ry/bohr^2].
        Used to reset the psi after the symmetrization.
       """

        new_v = np.einsum("a, ab->b",  np.sqrt(self.m) * vector, self.pols)
        
        new_dyn = dyn / np.outer(np.sqrt(self.m), np.sqrt(self.m))
        new_dyn = np.einsum("ab, ai, bj-> ij", new_dyn, self.pols, self.pols) 
        
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

        print ("Vector before symmetries:")
        print (vector)

        # Symmetrize the vector
        self.qe_sym.SetupQPoint()
        new_v = np.zeros( (self.nat, 3), dtype = np.float64, order = "F")
        new_v[:,:] = vector.reshape((self.nat, 3))
        self.qe_sym.SymmetrizeVector(new_v)
        vector = new_v.ravel()

        print ("Vector after symmetries:")
        print (vector)
               
        # Symmetrize the dynamical matrix
        dyn_q = CC.Phonons.GetDynQFromFCSupercell(dyn, np.array(self.dyn.q_tot), self.uci_structure, self.super_structure)
        self.qe_sym.SymmetrizeFCQ(dyn_q, self.dyn.q_stars, asr = "custom")
        dyn = CC.Phonons.GetSupercellFCFromDyn(dyn_q, np.array(self.dyn.q_tot), self.uci_structure, self.super_structure)

        # Push everithing back into the psi
        self.set_psi_from_vector_dyn(vector, dyn)

    def set_max_frequency(self, freq):
        """
        SETUP THE REVERSE LANCZOS
        =========================

        This function prepares the Lanczos algorithm in order to find the lowest eigenvalues
        You should provide the maximum frequencies of the standard spectrum.
        Then the Lanczos is initialized in order to solve the problem

        (Ia - L) x = b

        where a is sqrt(2*freq) so that should match the maximum energy. 
        Since Lanczos is very good in converging the biggest (in magnitude) eigenvectors, this
        procedure should accelerate the convergence of the low energy spectrum.

        NOTE: This method should be executed BEFORE the Lanczos run.

        Parameters
        ----------
            freq : float
               The frequencies (in Ry) of the maximum eigenvalue.
        """

        self.shift_value = 4*freq*freq
        self.reverse_L = True

        print("Shift value:", self.shift_value)

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
        out_vect[:self.n_modes] = (self.psi[:self.n_modes] * self.w) * self.w

        #print("freqsL1: ", self.w)
        #print("out 0:", out_vect[0])
        # Get the harmonic responce on the propagator
        w_a = np.tile(self.w, (self.n_modes, 1))
        w_b = np.tile(self.w, (self.n_modes, 1)).T 

        new_out = (w_a + w_b)**2
        out_vect[self.n_modes:] = new_out.ravel() * self.psi[self.n_modes:]

        #print("out 0 (just end):", out_vect[0])
        return out_vect

    def apply_L1_FT(self, transpose = False):
        """
        APPLY THE L1 AT FINITE TEMPERATURE
        ==================================

        This is the first part of the application, it involves only harmonic propagation.

        Results
        -------
            out_vect : ndarray(shape(self.psi))
                It returns the application of the harmonic part of the L matrix
        """


        # Prepare the free propagator on the positions
        out_vect = np.zeros(np.shape(self.psi), dtype = TYPE_DP)

        if self.ignore_harmonic:
            return out_vect 

        # The elements where w_a and w_b are exchanged are dependent
        # So we must avoid including them
        i_a = np.tile(np.arange(self.n_modes), (self.n_modes,1)).ravel()
        i_b = np.tile(np.arange(self.n_modes), (self.n_modes,1)).T.ravel()

        new_i_a = np.array([i_a[i] for i in range(len(i_a)) if i_a[i] >= i_b[i]])
        new_i_b = np.array([i_b[i] for i in range(len(i_a)) if i_a[i] >= i_b[i]])
        
        w_a = self.w[new_i_a]
        w_b = self.w[new_i_b]

        N_w2 = len(w_a)

        # Get the harmonic responce function
        out_vect[:self.n_modes] = (self.psi[:self.n_modes] * self.w) * self.w


        n_a = np.zeros(np.shape(w_a), dtype = TYPE_DP)
        n_b = np.zeros(np.shape(w_a), dtype = TYPE_DP)
        if self.T > 0:
            n_a = 1 / (np.exp( w_a / np.double(CC.Units.K_B * self.T)) - 1)
            n_b = 1 / (np.exp( w_b / np.double(CC.Units.K_B * self.T)) - 1)


        # Apply the non interacting X operator
        start_Y = self.n_modes
        start_A = self.n_modes + N_w2

        print("start_Y: {} | start_A: {} | end_A: {} | len_psi: {}".format(start_Y, start_A, start_A + N_w2, len(self.psi)))

        ERR_MSG ="""
ERROR,
The initial vector for the Lanczos algorithm has a wrong dimension. 
This may be caused by the Lanczos initialized at the wrong temperature.
"""
        assert len(self.psi) == start_A + N_w2, ERR_MSG

        # Apply the free propagation
        X_ab_NI = -w_a**2 - w_b**2 - (2*w_a *w_b) /( (2*n_a + 1) * (2*n_b + 1))
        out_vect[start_Y: start_A] = - X_ab_NI * self.psi[start_Y: start_A]

        # Perform the same on the A side
        Y_ab_NI = - (8 * w_a * w_b) / ( (2*n_a + 1) * (2*n_b + 1))
        if not transpose:
            out_vect[start_Y : start_A] += - Y_ab_NI * self.psi[start_A: ]
        else:
            out_vect[start_A:] += - Y_ab_NI * self.psi[start_Y : start_A]

        #L_operator[start_Y : start_A, start_A:] = - np.diag(Y_ab_NI) * extra_count
        #L_operator[start_Y + np.arange(self.n_modes**2), start_A + exchange_frequencies] -=  Y_ab_NI / 2

        X1_ab_NI = - (2*n_a*n_b + n_a + n_b) * (2*n_a*n_b + n_a + n_b + 1)*(2 * w_a * w_b) / ( (2*n_a + 1) * (2*n_b + 1))

        if not transpose:
            out_vect[start_A:] += - X1_ab_NI * self.psi[start_Y: start_A]
        else:
            out_vect[start_Y: start_A] += - X1_ab_NI * self.psi[start_A:]
        #L_operator[start_A:, start_Y : start_A] = - np.diag(X1_ab_NI) / 1 * extra_count
        #L_operator[start_A + np.arange(self.n_modes**2), start_Y + exchange_frequencies] -= X1_ab_NI / 2

        Y1_ab_NI = - w_a**2 - w_b**2 + (2*w_a *w_b) /( (2*n_a + 1) * (2*n_b + 1))
        out_vect[start_A:] += - Y1_ab_NI * self.psi[start_A:]
        #L_operator[start_A:, start_A:] = -np.diag(Y1_ab_NI) / 1 * extra_count
        #L_operator[start_A + np.arange(self.n_modes**2),  start_A + exchange_frequencies] -= Y1_ab_NI / 2


        return out_vect

    def apply_L2(self):
        """
        APPLY THE L2
        ============

        L2 is the part of the L operators that mixes the two spaces.
        It involves the phi3 matrix.
        """

        if self.ignore_v3:
            return np.zeros(np.shape(self.psi), dtype = TYPE_DP)



        w_a = np.tile(self.w, (self.n_modes, 1)).ravel()
        w_b = np.tile(self.w, (self.n_modes, 1)).T.ravel()

        vector = self.psi[:self.n_modes]
        dyn = self.psi[self.n_modes:]
        new_dyn = -dyn * np.sqrt( (w_a + w_b)/(w_a*w_b)) / 2

        # Here the time consuming part
        if self.mode == 0:
            # DEBUGGING PYTHON VERSION (SLOW)
            out_v = SlowApplyD3ToDyn(self.X, self.Y, self.rho, self.w, self.T, new_dyn)
            out_d = SlowApplyD3ToVector(self.X, self.Y, self.rho, self.w, self.T, vector)
        elif self.mode >= 1:
            out_v = FastApplyD3ToDyn(self.X, self.Y, self.rho, self.w, self.T, new_dyn, self.symmetries,
                                     self.N_degeneracy, self.degenerate_space, self.mode)
            out_d = FastApplyD3ToVector(self.X, self.Y, self.rho, self.w, self.T, vector, self.symmetries,
                                        self.N_degeneracy, self.degenerate_space, self.mode)
        else:
            print("Error, mode %d not recognized." % self.mode)
            raise ValueError("Mode not recognized %d" % self.mode)
            
        out_d *= -np.sqrt( (w_a + w_b)/(w_a*w_b)) / 2

        out_vect = np.zeros(np.shape(self.psi), dtype = TYPE_DP)
        
        out_vect[:self.n_modes] = out_v
        out_vect[self.n_modes:] = out_d
        return out_vect

    def apply_L2_FT(self, transpose = False):
        """
        Apply the full matrix at finite temperature.
        """ 

        if self.ignore_v3:
            return np.zeros(self.psi.shape, dtype = TYPE_DP)

        return FastD3_FT(self.X, self.Y, self.rho, self.w, self.T, self.psi, self.symmetries, self.N_degeneracy, self.degenerate_space, self.mode, transpose)

    def apply_L3(self):
        """
        APPLY THE L3
        ============

        This is the last part of the L matrix, it puts in communication 
        the dyn part of psi with herselfs.
        """

        w_a = np.tile(self.w, (self.n_modes, 1)).ravel()
        w_b = np.tile(self.w, (self.n_modes, 1)).T.ravel()

        simple_output = np.zeros(np.shape(self.psi), dtype = TYPE_DP)

        #simple_output[self.n_modes:] = self.psi[self.n_modes:] * (w_a + w_b)**2

        if self.ignore_v4:
            return simple_output

        # Apply the D4
        
        dyn = self.psi[self.n_modes:] * np.sqrt((w_a + w_b) / (w_a * w_b)) / 2
        


        # Here the time consuming part [The most of all]!!!
        if self.mode == 0:
            # A very slow implementation
            # Use it just for debugging
            out_dyn = SlowApplyD4ToDyn(self.X, self.Y, self.rho, self.w, self.T, dyn)
        elif self.mode >= 1:
            # The fast C implementation
            #print ("Inside v4 MPI, this will take a while")
            out_dyn = FastApplyD4ToDyn(self.X, self.Y, self.rho, self.w, self.T, dyn,
                                       self.symmetries, self.N_degeneracy, self.degenerate_space, self.mode)       
            
        out_dyn *= np.sqrt((w_a + w_b) / (w_a * w_b)) / 2

        output = np.zeros(np.shape(self.psi), dtype = TYPE_DP)
        output[self.n_modes:] = out_dyn

        output += simple_output

        return output


    def apply_L3_FT(self, transpose = False):
        """
        APPLY THE L3
        ============

        This is the last part of the L matrix, it puts in communication 
        the dyn part of psi with herselfs.
        """

        simple_output = np.zeros(np.shape(self.psi), dtype = TYPE_DP)

        #simple_output[self.n_modes:] = self.psi[self.n_modes:] * (w_a + w_b)**2

        if not self.ignore_v4:
            simple_output[:] = FastD4_FT(self.X, self.Y, self.rho, self.w, self.T, self.psi, self.symmetries, self.N_degeneracy, self.degenerate_space, self.mode, transpose)


        return simple_output

    def apply_full_L(self, target=None, force_t_0 = False, force_FT = False, transpose = False):
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
            force_t_0 : bool
                If False (default) the temperature is looked to chose if use the T = 0 or the finite temperature.
                If True it is forced the T=0 method (This will lead to wrong results at finite temperature).
            force_FT : bool
                If True the finite temperature method is forced even if T = 0.
                The results should be good, use it for testing.
                NOTE: only one between force_t_0 and force_FT should be true

        """

        if force_t_0 and force_FT:
            raise ValueError("Error, only one between force_t_0 and force_FT can be True")

        # Setup the target vector to the self.psi
        if not target is None:
            self.psi = target 

        # Check the initialization
        if not self.initialized:
            raise ValueError("Error, this class must be initialized before lunching the computation:\n Use the .init()")

        #if self.symmetrize:
        #    self.symmetrize_psi()

        #print("Psi before:")
        #print(self.psi[:self.n_modes])
            
            
        # Apply the whole L step by step to self.psi
        t1 = timer()
        if (force_t_0 or self.T < __EPSILON__) and not force_FT:
            output = self.apply_L1()
        else:
            output = self.apply_L1_FT(transpose)
        t2 = timer()
        if (force_t_0 or self.T < __EPSILON__) and not force_FT:
            output += self.apply_L2()
        else:
            output += self.apply_L2_FT(transpose)
        t3 = timer()
        if (force_t_0 or self.T < __EPSILON__) and not force_FT:
            output += self.apply_L3()
        else:
            output += self.apply_L3_FT(transpose)
        t4 = timer()

        print("Time to apply L1: {}".format(t2 - t1))
        print("Time to apply L2: {}".format(t3-t2))
        print("Time to apply L3: {}".format(t4-t3))

        # Apply the shift reverse
        #print ("Output before:")
        #print (output[:self.n_modes])
        if self.reverse_L:
            output *= -1
        output += self.shift_value * self.psi
        #print ("Output after:")
        #print (output[:self.n_modes])

        # Now return the output
        #print ("out just before return:", output[0])
        self.psi = output
        #if self.symmetrize:
        #    self.symmetrize_psi()
        
        return self.psi


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
        if Parallel.am_i_the_master():
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
                                rho = self.rho,
                                X = self.X,
                                Y = self.Y,
                                psi = self.psi,
                                a_coeffs = self.a_coeffs,
                                b_coeffs = self.b_coeffs,
                                krilov_basis = self.krilov_basis,
                                arnoldi_matrix = self.arnoldi_matrix,
                                reverse = self.reverse_L,
                                shift = self.shift_value,
                                symmetries = self.symmetries,
                                N_degeneracy = self.N_degeneracy,
                                initialized = self.initialized,
                                degenerate_space = self.degenerate_space,
                                perturbation_modulus = self.perturbation_modulus,
                                q_vectors = self.q_vectors)
            
    def load_status(self, file):
        """
        Load a previously saved status from the speficied npz file.
        The file must be saved with save_status.
        """

        # Check if the provided file exists
        if not os.path.exists(file):
            print ("Error while loading %s file.\n" % file)
            raise IOError("Error while loading %s" % file)

        
        # Fix the allow pickle error in numpy >= 1.14.4
        try:
            data = np.load(file, allow_pickle = True)
        except:
            data = np.load(file) 

        self.T = data["T"]
        self.nat = data["nat"]
        self.m = data["m"]
        self.w = data["w"]
        self.pols = data["pols"]
        self.n_modes = data["n_modes"]
        self.ignore_v3 = data["ignore_v3"]
        self.ignore_v4 = data["ignore_v4"]
        self.N = data["N"]
        self.rho = data["rho"]
        self.X = data["X"]
        self.Y = data["Y"]
        self.psi = data["psi"]
        self.a_coeffs = data["a_coeffs"]
        self.b_coeffs = data["b_coeffs"]
        self.krilov_basis = data["krilov_basis"]
        self.arnoldi_matrix = data["arnoldi_matrix"]

        if "reverse" in data.keys():
            self.reverse_L = data["reverse"]
            self.shift_value = data["shift"]

        if "symmetries" in data.keys():
            self.symmetries = data["symmetries"]
            self.N_degeneracy = data["N_degeneracy"]
            self.initialized = data["initialized"]
            self.degenerate_space = data["degenerate_space"]
        
        if "perturbation_modulus" in data.keys():
            self.perturbation_modulus = data["perturbation_modulus"]
            self.q_vectors = data["q_vectors"]

        # Rebuild the Linear operator
        self.L_linop = scipy.sparse.linalg.LinearOperator(shape = (len(self.psi), len(self.psi)), matvec = self.apply_full_L, dtype = TYPE_DP)


    def run_conjugate_gradient(self, eigval = 0, n_iters = 100, thr = 1e-5, verbose = True, guess_x = None):
        r"""
        RUN THE CONJUGATE GRADIENT
        ==========================

        The conjugate gradient is a very fast algorithm 
        that allows to compute the (in principle) exact green function. 

        Given the initial vector in the psi direction, it will edit it to obtain:
        
        .. math ::

            \left|x\right> = \left(\mathcal{L} - I\lambda\right)^{-1} \left| {\psi} \right> 

        At the end of the algorithm, the self.psi variable will contain the :math:`\left|x\right>`
        vector.

        Parameters
        ----------
            eigval : float
                The value of the :math:`\lambda` for the inversion problem.
            n_iters : int
                The number of iterations
            thr : float
                The threshold between two iteration after which the algorithm is considered 
                to be converged.
            verbose : bool
                If true print the status of the iteration.
                
        """

        if guess_x is None:
            x = self.psi.copy()
        else:
            x = guess_x.copy()
        A_x = self.apply_full_L()

        r = x - A_x
        p = r.copy()

        mod_r = np.sqrt(r.dot(r))

        if verbose:
            print("   CG step %d : residual = %.4e | threshold = %.4e" % (0,mod_r, thr))

        if mod_r < thr:
            return x

        for i in range(n_iters):
            self.psi = p
            A_p = self.apply_full_L()
            alpha = r.dot(r) / p.dot(A_p)
            x += alpha * p 

            # Update
            r_new = r - alpha * A_p 
            beta = r_new.dot(r_new) / r.dot(r)

            r = r_new
            p = r + beta * p 


            # Check the new iteration
            mod_r = np.sqrt(r.dot(r))
            if verbose:
                print("   CG step %d : residual = %.4e | threshold = %.4e" % (i+1,mod_r, thr))

            if mod_r < thr:
                return x

        print("WARNING: CG ended before the convergence was achieved.") 
        return x

    def get_statical_responce_from_scratch(self, n_iters = 100, thr = 1e-5, verbose = True, sub_block = None, sub_space = None):
        """
        GET STATIC RESPONCE
        ===================

        This algorithm performs the CG minimization to obtain the static self-energy.

        Parameters
        ----------
            n_iters : int
                The number of maximum iteration for a single CG step.
            thr : float
                The threshold for the convergence of the CG algorithm.
            verbose : bool
                If true (default) prints the info during the minimization
            sub_block : list
                A list of indices that identifies the modes id that that you want to select.
                The algorithm is performed only in the reduced space of (N_modes x N_modes)
                given by the length of this list. 
                In this way you will neglect mode interaction, but you can save a lot of time.
                Leave as None if you want the whole space.
            sub_space : ndarray(size = (N_dim, 3*n_at))
                Compute the self-energy only in the subspace given. Leave it as none
                if you do not want to use this option

        Results
        -------
            fc_matrix : ndarray (size=(3*nat, 3*nat))
                The static self-energy.
        """

        n_dim_space = self.n_modes


        PLinvP = np.zeros((n_dim_space, n_dim_space), dtype = TYPE_DP, order = "C")

        if (not sub_block is None) and (not sub_space is None):
            raise ValueError("Error, you cannot specify both sub_block and sub_space.")

        if not sub_block is None:
            n_dim_space = len(sub_block)
        elif not sub_space is None:
            n_dim_space = len(sub_space)
            
        # initialize the algorithm
        for i in range(n_dim_space):
            # Setup the vector
            self.psi = np.zeros((self.n_modes + 1)*self.n_modes, dtype = TYPE_DP)
            guess = self.psi.copy()

            # Create the subbasis
            if not sub_block is None:
                self.psi[sub_block[i]] = 1
                guess[sub_block[i]] = self.w[sub_block[i]] ** 2
            elif not sub_space is None:
                self.psi[:self.n_modes] = sub_block[i].dot(self.pols)
                guess = None
            else:
                self.psi[i] = 1
                guess[i] = self.w[i] ** 2

            if verbose:
                print("")
                print("==== NEW STATIC COMPUTATION ====")
                print("Iteration: %d out of %d" % (i+1, n_dim_space))


            new_v = self.run_conjugate_gradient(n_iters = n_iters, thr = thr, verbose = verbose, guess_x = guess)
            
            if not sub_space is None:
                v_out = self.pols.dot(new_v[:self.n_modes])
                PLinvP[i, :] = np.array(sub_space).dot(v_out)
            else:
                PLinvP[i, :] = new_v[:self.n_modes]

        

        # Invert the P L^-1 P 
        D = np.linalg.inv(PLinvP)
        # Transform to a force constant matrix in cartesian coordinates

        if not sub_space is None:
            fc_matrix = np.einsum("ab, ai, bj->ij", D, np.array(sub_space), np.array(sub_space))
        else:
            fc_matrix = np.einsum("ab, ia, jb->ij", D, self.pols, self.pols)
        fc_matrix *= np.sqrt(np.outer(self.m, self.m))

        return fc_matrix




    def run(self, n_iter, save_dir = ".", verbose = True):
        """
        RUN LANCZOS ITERATIONS
        ======================

        This method performs the Lanczos algorithm to find
        the sequence of a and b coefficients that are the tridiagonal representation 
        of the L matrix to be inverted.

        Parameters
        ----------
            n_iter : int
                The number of iterations to be performed in the Lanczos algorithm.
            save_dir : string
                The directory in which you want to store the results step by step,
                in order to do a preliminar analysis or restart the calculation later.
            verbose : bool
                If true all the info during the minimization will be printed on output.
        """

        # Check if the symmetries has been initialized
        if not self.initialized:
            self.prepare_symmetrization()

        # Get the current step
        i_step = len(self.a_coeffs)

        if verbose:
            header = """
<=====================================>
|                                     |
|          LANCZOS ALGORITHM          |
|                                     |
<=====================================>

Starting the algorithm. It may take a while.
Starting from step %d
""" % i_step
            print(header)

            OPTIONS = """
Should I ignore the third order effect? {}
Should I ignore the fourth order effect? {}
Max number of iterations: {}
""".format(self.ignore_v3, self.ignore_v4, n_iter)
            print(OPTIONS)


        # If this is the current step initialize the algorithm
        if i_step == 0:
            self.krilov_basis = []
            first_vector = self.psi / np.sqrt(self.psi.dot(self.psi))
            self.krilov_basis.append(first_vector)
        else:
            # Convert everything in a list
            self.krilov_basis = list(self.krilov_basis)
            self.a_coeffs = list(self.a_coeffs)
            self.b_coeffs = list(self.b_coeffs)
            self.arnoldi_matrix = list(self.arnoldi_matrix)

            if len(self.krilov_basis) != i_step + 1:
                print("Krilov dim: %d, number of steps perfomed: %d" % (len(self.krilov_basis), i_step))
                print("Error, the krilov basis dimension should be 1 more than the number of steps")
                raise ValueError("Error the starting krilov basis does not matches the matrix, Look stdout.")

        self.psi = self.krilov_basis[-1]

        for i in range(i_step, i_step+n_iter):
            if verbose:
                step_txt = """
 ===== NEW STEP %d =====

 """ % i
                print(step_txt)
                print("Length of the coefficiets: a = {}, b = {}".format(len(self.a_coeffs), len(self.b_coeffs)))
                print()

            # Apply the matrix L
            t1 = time.time()
            #self.psi = self.apply_full_L()
            self.psi = self.L_linop.dot(self.psi)
            t2 = time.time()

            if verbose:
                print("Time to apply the full L: %d s" % (t2 -t1))

            # Get the coefficients for the Lanczos/Arnoldi matrix
            t1 = time.time()
            arnoldi_row = []
            new_vect = self.psi.copy()


            # Lets repeat twice the orthogonalization
            converged = False
            for k_orth in range(N_REP_ORTH):
                for j in range(len(self.krilov_basis)):
                    coeff = new_vect.dot(self.krilov_basis[j])
                    if k_orth == 0:
                        arnoldi_row.append(self.psi.dot(self.krilov_basis[j]))

                    # Gram Schmidt
                    new_vect -= coeff * self.krilov_basis[j]
            
                # Add the new vector to the Krilov Basis
                norm = np.sqrt(new_vect.dot(new_vect))

                if verbose:
                    print("Vector norm after GS number {}: {:16.8e}".format(k_orth, norm))

                # Check the normalization (If zero the algorithm converged)
                if norm < __EPSILON__:
                    converged = True
                    if verbose:
                        print("Obtained a linear dependent vector.")
                        print("The algorithm converged.")
                    break
                
                new_vect /= norm 

            if not converged:
                self.krilov_basis.append(new_vect)
                self.psi = new_vect
            t2 = time.time()

            # Add the coefficients to the variables
            self.a_coeffs.append(arnoldi_row[-1])
            if len(arnoldi_row) > 1:
                self.b_coeffs.append(arnoldi_row[-2])
            self.arnoldi_matrix.append(arnoldi_row)

            if verbose:
                print("Time to perform the Gram-Schmidt and retrive the coefficients: %d s" % (t2-t1))
                print()
                print("a_%d = %.8e" % (i, self.a_coeffs[-1]))
                if i > 0:
                    print("b_%d = %.8e" % (i, self.b_coeffs[-1]))
                print()
            
            # Save the step
            if not save_dir is None:
                self.save_status("%s/LANCZOS_STEP%d" % (save_dir, i))
        
                if verbose:
                    print("Status saved into '%s/LANCZOS_STEP%d'" % (save_dir, i))
            
            if verbose:
                print("Lanczos step %d ultimated." % i)
            

            if converged:
                return


    def build_lanczos_matrix_from_coeffs(self, use_arnoldi=False):
        """
        BUILD THE LANCZOS MATRIX
        ========================

        This method builds the Lanczos matrix from the coefficients. 
        To execute this method correctly you must have already completed the Lanczos algorithm (method run)

        Parameters
        ----------
            use_arnoldi: bool
                If true the full matrix is computed, using all the coefficients from the
                Arnoldi iteration.
        """

        N_size = len(self.a_coeffs)
        matrix = np.zeros((N_size, N_size), dtype = TYPE_DP)
        if not use_arnoldi:
            for i in range(N_size):
                matrix[i,i] = self.a_coeffs[i]
                if i>= 1:
                    # Use the non-symmetric Lanczos if also c_coeffs are present
                    c_coeff = self.b_coeffs[i-1]
                    if len(self.c_coeffs) == len(self.b_coeffs):
                        c_coeff = self.c_coeffs[i-1]
                    matrix[i-1,i] = c_coeff
                    matrix[i,i-1] = self.b_coeffs[i-1]
        else:
            # Check if there are c_coeffs, in this way arnoldi matrix is not computed
            assert len(self.b_coeffs) > len(self.c_coeffs), "Error, cannot Arnoldi with non-symmetric Lanczos not implemented"
            for i in range(N_size):
                matrix[:i+1, i] = self.arnoldi_matrix[i]
                matrix[i, :i+1] = self.arnoldi_matrix[i]


        sign = 1
        if self.reverse_L:
            sign = -1

        matrix =  sign*matrix - sign* np.eye(N_size) * self.shift_value
                    
        return matrix


    def get_green_function_Lenmann(self, w_array, smearing, v_a, v_b, use_arnoldi = False):
        """
        GET GREEN FUNCTION
        ==================

        Compute the green function using the Lemman representation.

        Parameters
        ----------
            w_array : ndarray
                The list of frequencies for which you want to compute the
                dynamical green function.
            smearing : float
                The smearing to take a non zero imaginary part.
            v_a : ndarray(size = 3*self.nat)
                The perturbation operator (on atomic positions)
            v_b : ndarray(size = 3*self.nat)
                The probed responce operator (on atomic positions)
            use_arnoldi: bool
                If true the full arnoldi matrix is used to extract eigenvalues and 
                eigenvectors. Otherwise the tridiagonal Lanczos matrix is used.
                The first one prevents the loss of orthogonality problem.
        """

        # Get the Lanczos matrix
        matrix = self.build_lanczos_matrix_from_coeffs(use_arnoldi)

        assert len(self.c_coeffs) < len(self.b_coeffs), "Lenmann cannot be used with non-symmetric Lanczos"

        # Convert the vectors in the polarization basis
        new_va = np.einsum("a, a, ab->b", 1/np.sqrt(self.m), v_a, self.pols)
        new_vb = np.einsum("a, a, ab->b", 1/np.sqrt(self.m), v_b, self.pols)

        # Dyagonalize the Lanczos matrix
        eigvals, eigvects = np.linalg.eigh(matrix)

        kb = np.array(self.krilov_basis)
        kb = kb[:-1,:]
        #print (np.shape(eigvects), np.shape(kb))
        # Convert in krilov space
        new_eigv = np.einsum("ab, ac->cb", eigvects, kb)


        Na, Nb = np.shape(matrix)
        if Na != Nb:
            raise ValueError("Error, the Lanczos matrix must be square, dim (%d,%d)" % (Na, Nb))
        
        gf = np.zeros(len(w_array), dtype = np.complex128)

        for j in range(Na):
            eig_v = new_eigv[:self.n_modes, j]
            matrix_element = eig_v.dot(new_va) * new_vb.dot(eig_v)
            gf[:] += matrix_element / (eigvals[j]  - w_array**2 + 2j*w_array*smearing)

        return gf

    def get_static_odd_fc(self, use_arnoldi = False):
        """
        GET STATIC FORCE CONSTANT
        =========================

        Get the static force constant matrix

        Parameters
        ----------
            use_arnoldi: bool
                If true the full arnoldi matrix is used, otherwise the Lanczos tridiagonal
                matrix is used.
        """

        # Get the Lanczos matrix
        matrix = self.build_lanczos_matrix_from_coeffs(use_arnoldi)

        # Dyagonalize the Lanczos matrix
        eigvals, eigvects = np.linalg.eigh(matrix)

        Nk = len(self.krilov_basis)

        kb = np.array(self.krilov_basis)
        
        # Lanczos did not converged, discard the last vector
        if Nk > len(eigvals):
            kb = kb[:-1,:]

        #print (np.shape(eigvects), np.shape(kb))
        new_eigv = np.einsum("ab, ac->cb", eigvects, kb)

        Na, Nb = np.shape(matrix)
        if Na != Nb:
            raise ValueError("Error, the Lanczos matrix must be square, dim (%d,%d)" % (Na, Nb))
        

        fc_matrix = np.zeros( (3*self.nat, 3*self.nat), dtype = TYPE_DP)

        # Get the dynamical matrix in the polarization basis
        D = np.einsum("ai, bi, i->ab", new_eigv[:self.n_modes,:], new_eigv[:self.n_modes, :], eigvals)

        # Convert it in the standard basis
        fc_matrix = np.einsum("ab, ia, jb->ij", D, self.pols, self.pols)

        # for i in range(3*self.nat):
        #     # Define the vector
        #     v = np.zeros(3*self.nat, dtype = TYPE_DP)
        #     v[i] = 1

        #     # Convert the vectors in the polarization basis
        #     new_v = np.einsum("a, a, ab->b", np.sqrt(self.m), v, self.pols)
        #     # Convert in the krilov space 
        #     mat_coeff = np.einsum("a, ab", new_v, new_eigv[:self.n_modes, :])
        #     new_w = np.einsum("a, ba, a", mat_coeff, new_eigv[:self.n_modes,:], eigvals)

        #     #v_kb = np.einsum("ab, b", kb[:, :self.n_modes], new_v)
        #     # Apply the L matrix
        #     #w_kb = matrix.dot(v_kb)
        #     # Convert back in the polarization space
        #     #new_w = np.einsum("ab, a", kb[:, :self.n_modes], w_kb)
        #     # Convert back in real spaceDoes anyone know if there is a windows binary or a source code to run QE with GPU enhancement on windows. 
        #     w = np.einsum("a, b, ab ->a", 1/np.sqrt(self.m), new_w, self.pols)

        #     fc_matrix[i, :] = w
            

        # This is the dynamical matrix now we can multiply by the masses
        fc_matrix *= np.sqrt(np.outer(self.m, self.m))

        return fc_matrix


    def get_all_green_functions(self, N_steps = 100, mode_mixing = True, save_step_dir = None, verbose = True):
        """
        GET ALL THE GREEN FUNCTIONS
        ===========================

        This will compute a set of lanczos coefficients for each element of the odd matrix.
        a_n and b_n.
        We will run lanczos for all the elements and all the crosses.
        In this way we have the whole evolution with frequency of the matrix.

        NOTE: This can be a very intensive computation.

        Parameters
        ----------
            N_steps : int
                The number of Lanczos iteration for each green function
            mode_mixing : bool
                If True also non diagonal elements are computed, otherwise the 
                SSCHA eigenvector are supposed to be conserved, and only diagonal
                green functions are considered.
                If False the computation is much less expensive (a factor nat_sc),
                but it is approximated.
            save_step_dir : string
                If not None, the path to the directory in which you want to save 
                each step. So even if stopped the calculation can restart.
            verbose : bool
                If true print all the progress to standard output

        Results
        -------
            a_ns : ndarray( (n_modes, n_modes, N_steps))
                The a coefficients for each element in the mode x mode space
            b_ns : ndarray( (n_modes, n_modes, N_steps-1))
                The b_n coefficients for each mode in the space.
        """

        # Time the function
        t_start = time.time()

        # Check if the save directory exists
        # Otherwise we create it
        if not save_step_dir is None:
            if not os.path.exists(save_step_dir):
                os.makedirs(save_step_dir)

        # Load all the data
        a_ns = np.zeros( (self.n_modes, self.n_modes, N_steps), dtype = np.double)
        b_ns = np.zeros( (self.n_modes, self.n_modes, N_steps-1), dtype = np.double)

        # Incompatible with shift for now
        self.shift_value = 0

        # Compute the diagonal parts
        for i in range(self.n_modes):
            if verbose:
                print("\n")
                print("  ==========================  ")
                print("  |                        |  ")
                print("  |   DIAGONAL ELEMENTS    |  ")
                print("  |       STEP {:5d}       |  ".format(i))
                print("  |                        |  ")
                print("  ==========================  ")
                print()
            
            # Setup the Lanczos
            self.reset()

            # Prepare the perturbation
            self.psi[:] = 0
            self.psi[i] = 1

            # Run the Lanczos perturbation
            self.run(N_steps, save_dir = save_step_dir, verbose = verbose)

            if verbose:
                print()
                print("   ---- > LANCZOS RUN COMPLEATED < ----   ")
                print()

            # Save the status
            if save_step_dir:
                self.save_status("full_lanczos_diagonal_{}".format(i))
            
            # Fill the a_n and b_n
            a_tmp = np.zeros(N_steps, dtype = np.double)
            a_tmp[:len(self.a_coeffs)] = self.a_coeffs
            b_tmp = np.zeros(N_steps-1, dtype = np.double)
            b_tmp[:len(self.b_coeffs)] = self.b_coeffs
            a_ns[i, i, :] = a_tmp
            b_ns[i, i, :] = b_tmp
    
        # If we must compute the mode mixing
        if mode_mixing:
            for i in range(self.n_modes):
                for j in range(i+1, self.n_modes):
                    # TODO: Neglect (i,j) forbidden by symmetries

                    if verbose:
                        print("\n")
                        print("  ============================  ")
                        print("  |                          |  ")
                        print("  |   NON DIAGONAL ELEMENT   |  ")
                        print("  |    STEP ({:5d},{:5d})    |  ".format(i, j))
                        print("  |                          |  ")
                        print("  ============================  ")
                        print()
                    
                    # Setup the Lanczos
                    self.reset()

                    # Prepare the perturbation
                    self.psi[:] = 0
                    self.psi[i] = 1
                    self.psi[j] = 1

                    # Run the Lanczos perturbation
                    self.run(N_steps, save_dir = save_step_dir, verbose = verbose)

                    if verbose:
                        print()
                        print("   ---- > LANCZOS RUN COMPLEATED < ----   ")
                        print()

                    # Save the status
                    if save_step_dir:
                        self.save_status("full_lanczos_off_diagonal_{}_{}".format(i, j))
                    
                    # Fill the a_n and b_n
                    a_tmp = np.zeros(N_steps, dtype = np.double)
                    a_tmp[:len(self.a_coeffs)] = self.a_coeffs
                    b_tmp = np.zeros(N_steps-1, dtype = np.double)
                    b_tmp[:len(self.b_coeffs)] = self.b_coeffs
                    a_ns[i, j, :] = a_tmp
                    b_ns[i, j, :] = b_tmp
                    a_ns[j, i, :] = a_tmp
                    b_ns[j, i, :] = b_tmp

        t_end = time.time()

        total_time = t_end - t_start
        minutes = int(total_time / 60)
        hours = int(minutes / 60)
        minutes -= hours * 60
        seconds = int(total_time - hours*3600 - minutes * 60)

        if verbose:
            print()
            print()
            print("     ======================     ")
            print("     |                    |     ")
            print("     |        DONE        |      ")
            print("     |   In {:3d}:{:02d}:{:02d}s  |     ".format(hours, minutes, seconds))
            print("     ======================     ")
            print()
            print()
            
        return a_ns, b_ns



    def get_spectral_function_from_Lenmann(self, w_array, smearing, use_arnoldi=True):
        """
        GET SPECTRAL FUNCTION
        =====================

        This method computes the spectral function in the supercell
        using the Lenmann representation.

        Parameters
        ----------
            w_array : ndarray
                The list of frequencies for which you want to compute the
                dynamical green function.
            smearing : float
                The smearing to take a non zero imaginary part.
            use_arnoldi: bool
                If true the full arnoldi matrix is used to extract eigenvalues and 
                eigenvectors. Otherwise the tridiagonal Lanczos matrix is used.
                The first one prevents the loss of orthogonality problem.
        """
        # Get the Lanczos matrix
        matrix = self.build_lanczos_matrix_from_coeffs(use_arnoldi)

        # Dyagonalize the Lanczos matrix
        eigvals, eigvects = np.linalg.eigh(matrix)

        Na, Nb = np.shape(matrix)
        if Na != Nb:
            raise ValueError("Error, the Lanczos matrix must be square, dim (%d,%d)" % (Na, Nb))
        
        spectral = np.zeros(len(w_array), dtype = np.complex128)

        kb = np.array(self.krilov_basis)
        if np.shape(kb)[0] > Na:
            kb = kb[:-1,:]
        print ("Shape check: eigvects = {}, kb = {}".format( np.shape(eigvects), np.shape(kb)))
        new_eigv = np.einsum("ab, ac->cb", eigvects, kb)

        for j in range(Na):
            eig_v = new_eigv[:self.n_modes, j]
            matrix_element = eig_v.dot(eig_v)
            spectral[:] += matrix_element / (eigvals[j]  - w_array**2 +2j*w_array*smearing)

        return -np.imag(spectral)


    def get_green_function_continued_fraction(self, w_array, use_terminator = True, last_average = 1, smearing = 0):
        """
        CONTINUED FRACTION GREEN FUNCTION
        =================================

        In this way the continued fraction for the green function is used.
        This should converge faster than the Lenmann representation, and
        has the advantage of adding the possibility to add a terminator.
        This avoids to define a smearing.

        Parameters
        ----------
            w_array : ndarray
                The list of frequencies in which you want to compute the green function
            use_terminator : bool
                If true (default) a standard terminator is used.
            last_average : int
                How many a and be coefficients are averaged to evaluate the terminator?
            smearing : float
                The smearing parameter. If none
        """

        n_iters = len(self.a_coeffs)

        gf = np.zeros(np.shape(w_array), dtype = np.complex128)

        sign =1
        if self.reverse_L:
            sign = -1

        # Get the terminator
        if use_terminator:
            a_av = np.mean(self.a_coeffs[-last_average:])
            b_av = np.mean(self.b_coeffs[-last_average:])
            c_av = b_av
            if len(self.c_coeffs) == len(self.b_coeffs): # Non-symmetric Lanczos
                c_av = np.mean(self.c_coeffs[-last_average:])

            a = a_av * sign - sign* self.shift_value
            b = b_av * sign
            c = c_av * sign

            gf[:] = (a - w_array**2 - np.sqrt( (a - w_array**2)**2 - 4*b*c + 0j))/(2*b*c)
        else:
            a = self.a_coeffs[-1] * sign - sign* self.shift_value
            gf[:] = 1/ (a - w_array**2 + 2j*w_array*smearing)

        for i in range(n_iters-2, -1, -1):
            a = self.a_coeffs[i] * sign - sign* self.shift_value
            b = self.b_coeffs[i] * sign
            c = b
            if len(self.c_coeffs) == len(self.b_coeffs): # Non-symmetric Lanczos
                c = self.c_coeffs[i] * sign
            gf = 1. / (a - w_array**2  + 2j*w_array*smearing - b*c * gf)

        return gf * self.perturbation_modulus

    
    def get_full_L_operator(self, verbose = False, only_pert=False, debug_d3 = None):
        """
        GET THE FULL L OPERATOR
        =======================
        
        Use this method to test everithing. I returns the full L operator as a matrix.
        It is very memory consuming, but it should be fast and practical for small systems.
        

        Results
        -------
           L_op : ndarray(size = (nmodes * (nmodes + 1)), dtype = TYPE_DP)
              The full L operator.
        """


        L_operator = np.zeros( shape = (self.n_modes + self.n_modes * self.n_modes, self.n_modes + self.n_modes * self.n_modes), dtype = TYPE_DP)

        # Fill the first part with the standard dynamical matrix
        if not only_pert:
            L_operator[:self.n_modes, :self.n_modes] = np.diag(self.w**2)

        
        w_a = np.tile(self.w, (self.n_modes,1)).ravel()
        w_b = np.tile(self.w, (self.n_modes,1)).T.ravel()
        chi_beta = -.5 * np.sqrt(w_a + w_b)/(np.sqrt(w_a)*np.sqrt(w_b))


        B_mat = (w_a + w_b)**2
        if not only_pert:
            L_operator[self.n_modes:, self.n_modes:] = np.diag(B_mat)
        

        # Compute the d3 operator
        #new_X = np.einsum("ia,a->ai", self.X, f_ups(self.w, self.T))
        if debug_d3 is None:
            N_eff = np.sum(self.rho)
            Y_weighted = np.einsum("ia, i->ia", self.Y, self.rho)
            #new_Y = np.einsum("ia,i->ai", self.Y, self.rho)

        if not self.ignore_v3:

            if not debug_d3 is None:
                d3 = debug_d3
            else:
                if verbose:
                    print("Computing d3...")
                d3_noperm = np.einsum("ia,ib,ic->abc", self.X, self.X, Y_weighted)
                d3_noperm /= -N_eff 

                # Apply the permuatations
                d3 = d3_noperm.copy()
                d3 += np.einsum("abc->acb", d3_noperm)
                d3 += np.einsum("abc->bac", d3_noperm)
                d3 += np.einsum("abc->bca", d3_noperm)
                d3 += np.einsum("abc->cab", d3_noperm)
                d3 += np.einsum("abc->cba", d3_noperm)
                d3 /= 6

                if verbose:
                    np.save("d3_modes_nosym.npy", d3)

                # Perform the standard symmetrization
                d3 = symmetrize_d3_muspace(d3, self.symmetries)

                if verbose:
                    np.save("d3_modes_sym.npy", d3)
                    np.save("symmetries_modes.npy", self.symmetries)
            

            # Reshape the d3
            d3_reshaped = d3.reshape((self.n_modes, self.n_modes * self.n_modes))

            new_mat = np.einsum("ab,b->ab", d3_reshaped, chi_beta)

            L_operator[:self.n_modes, self.n_modes:] = new_mat
            L_operator[self.n_modes:, :self.n_modes] = new_mat.T
        if not self.ignore_v4:

            if verbose:
                print("Computing d4...")
            d4 =  np.einsum("ai,bi,ci,di", new_X, new_X, new_X, new_Y)
            d4 += np.einsum("ai,bi,ci,di", new_X, new_X, new_Y, new_X)
            d4 += np.einsum("ai,bi,ci,di", new_X, new_Y, new_X, new_X)
            d4 += np.einsum("ai,bi,ci,di", new_Y, new_X, new_X, new_X)
            d4 /= - 4 * N_eff

            if verbose:
                np.save("d4_modes_nosym.npy", d4)

            # Reshape the d4
            d4_reshaped = d4.reshape((self.n_modes*self.n_modes, self.n_modes * self.n_modes))

            new_mat = np.einsum("ab,a,b->ab", d4_reshaped, chi_beta, chi_beta)

            L_operator[self.n_modes:, self.n_modes:] += new_mat

        if verbose:
            print("L superoperator computed.")
        
        self.L_linop = L_operator
        return L_operator


    def get_full_L_operator_FT(self, verbose = False, debug_d3 = None):
        """
        GET THE FULL L OPERATOR (FINITE TEMPERATURE)
        ============================================

        Get the the full matrix L for the biconjugate Lanczos algorithm.
        Use this method to test everithing. I returns the full L operator as a matrix.
        It is very memory consuming, but it should be fast for small systems.

        Maybe we need to drop the exchange between a,b because they are symmetric by definition.

        Results
        -------
           L_op : ndarray(size = (nmodes * (2*nmodes + 1)), dtype = TYPE_DP)
              The full L operator.
        """

        # The elements where w_a and w_b are exchanged are dependent
        # So we must avoid including them
        i_a = np.tile(np.arange(self.n_modes), (self.n_modes,1)).ravel()
        i_b = np.tile(np.arange(self.n_modes), (self.n_modes,1)).T.ravel()

        new_i_a = np.array([i_a[i] for i in range(len(i_a)) if i_a[i] >= i_b[i]])
        new_i_b = np.array([i_b[i] for i in range(len(i_a)) if i_a[i] >= i_b[i]])
        
        w_a = self.w[new_i_a]
        w_b = self.w[new_i_b]

        N_w2 = len(w_a)

        # Prepare the operator
        L_operator = np.zeros( shape = (self.n_modes + 2*N_w2, self.n_modes + 2*N_w2), dtype = TYPE_DP)

        # Set the Z''
        if not self.ignore_harmonic:
            L_operator[:self.n_modes, :self.n_modes] = np.diag(self.w**2)


        #w_a = np.tile(self.w, (self.n_modes,1)).ravel()
        #w_b = np.tile(self.w, (self.n_modes,1)).T.ravel()

        n_a = np.zeros(np.shape(w_a), dtype = TYPE_DP)
        n_b = np.zeros(np.shape(w_a), dtype = TYPE_DP)
        if self.T > 0:
            n_a = 1 / (np.exp( w_a / np.double(CC.Units.K_B * self.T)) - 1)
            n_b = 1 / (np.exp( w_b / np.double(CC.Units.K_B * self.T)) - 1)

        print("NA", n_a[:self.n_modes])

        # Apply the non interacting X operator
        start_Y = self.n_modes
        start_A = self.n_modes + N_w2

        # Since we excluded the w_b < w_a, when w_a = w_b we have a double count
        extra_count = np.ones(N_w2, dtype = np.intc)
        extra_count[new_i_a == new_i_b] = 1.

        # Get the operator that exchanges the frequencies
        # For each index i (a,b), exchange_frequencies[i] is the index that correspond to (b,a)
        #exchange_frequencies = np.array([ (i // self.n_modes) + self.n_modes * (i % self.n_modes) for i in np.arange(self.n_modes**2)])
        #xx = np.tile(np.arange(self.n_modes), (self.n_modes, 1)).T.ravel()
        #yy = np.tile(np.arange(self.n_modes), (self.n_modes, 1)).ravel()
        #all_modes = np.arange(self.n_modes**2)
        #exchange_frequencies = xx + yy
        if not self.ignore_harmonic:
            # Apply the operator to himself and to the exchange on the frequencies.
            X_ab_NI = -w_a**2 - w_b**2 - (2*w_a *w_b) /( (2*n_a + 1) * (2*n_b + 1))
            L_operator[start_Y: start_A, start_Y:start_A] = - np.diag(X_ab_NI)  * extra_count
            #L_operator[start_Y + np.arange(self.n_modes**2) , start_Y + exchange_frequencies] -= X_ab_NI / 2

            # Perform the same on the A side
            Y_ab_NI = - (8 * w_a * w_b) / ( (2*n_a + 1) * (2*n_b + 1))
            L_operator[start_Y : start_A, start_A:] = - np.diag(Y_ab_NI) * extra_count
            #L_operator[start_Y + np.arange(self.n_modes**2), start_A + exchange_frequencies] -=  Y_ab_NI / 2

            X1_ab_NI = - (2*n_a*n_b + n_a + n_b) * (2*n_a*n_b + n_a + n_b + 1)*(2 * w_a * w_b) / ( (2*n_a + 1) * (2*n_b + 1))
            L_operator[start_A:, start_Y : start_A] = - np.diag(X1_ab_NI) / 1 * extra_count
            #L_operator[start_A + np.arange(self.n_modes**2), start_Y + exchange_frequencies] -= X1_ab_NI / 2

            Y1_ab_NI = - w_a**2 - w_b**2 + (2*w_a *w_b) /( (2*n_a + 1) * (2*n_b + 1))
            L_operator[start_A:, start_A:] = -np.diag(Y1_ab_NI) / 1 * extra_count
            #L_operator[start_A + np.arange(self.n_modes**2),  start_A + exchange_frequencies] -= Y1_ab_NI / 2

        # We added all the non interacting propagators

        # Compute the d3 operator
        #new_X = np.einsum("ia,a->ai", self.X, f_ups(self.w, self.T))
        if debug_d3 is None:
            N_eff = np.sum(self.rho)
            Y_weighted = np.einsum("ia, i->ia", self.Y, self.rho)
        #new_Y = np.einsum("ia,i->ai", self.Y, self.rho)

        if not self.ignore_v3:
            if verbose:
                print("Computing d3...")
            if not debug_d3 is None:
                d3 = debug_d3
            else:
                d3_noperm = np.einsum("ia,ib,ic->abc", self.X, self.X, Y_weighted)
                d3_noperm /= -N_eff 

                # Apply the permuatations
                d3 = d3_noperm.copy()
                d3 += np.einsum("abc->acb", d3_noperm)
                d3 += np.einsum("abc->bac", d3_noperm)
                d3 += np.einsum("abc->bca", d3_noperm)
                d3 += np.einsum("abc->cab", d3_noperm)
                d3 += np.einsum("abc->cba", d3_noperm)
                d3 /= 6

                if verbose:
                    np.save("d3_modes_nosym.npy", d3)

                # Perform the standard symmetrization
                d3 = symmetrize_d3_muspace(d3, self.symmetries)

                if verbose:
                    np.save("d3_modes_sym.npy", d3)
                    np.save("symmetries_modes.npy", self.symmetries)
            

            # Reshape the d3
            d3_small_space = np.zeros((N_w2, self.n_modes), dtype = np.double)
            d3_small_space[:,:] = d3[new_i_a, new_i_b, :]

            print("D3 of the following elements:")
            print(new_i_a)
            print(new_i_b)
            print(d3_small_space)

            #d3_reshaped = d3.reshape((self.n_modes* self.n_modes, self.n_modes))
            #d3_reshaped1 = d3.reshape((self.n_modes, self.n_modes* self.n_modes))
            
            # Get the Z coefficient
            Z_coeff = 2 * ((2*n_a + 1)*w_b + (2*n_b + 1)*w_a) / ((2*n_a + 1) * (2*n_b + 1))
            Z_coeff = np.einsum("ab,a -> ab", d3_small_space, Z_coeff)
            L_operator[start_Y: start_A, :start_Y] = -Z_coeff

            # Get the Z' coefficients
            Z1_coeff = 2 * ( (2*n_a + 1)*w_b*n_b*(n_b + 1) + (2*n_b + 1)*w_a*n_a*(n_a+1)) / ((2*n_a + 1) * (2*n_b + 1))
            Z1_coeff = np.einsum("ab,a -> ab", d3_small_space, Z1_coeff)
            L_operator[start_A:, :start_Y] = - Z1_coeff

            # The other coeff between Y and R
            # X''
            extra_count = np.ones(N_w2, dtype = np.intc)
            extra_count[new_i_a != new_i_b] = 2
            X2_coeff = (2*n_b + 1) * (2*n_a +1) / (8*w_a *w_b)
            X2_coeff = np.einsum("ab,a->ba", d3_small_space, X2_coeff * extra_count)
            L_operator[:start_Y, start_Y: start_A] = -X2_coeff

            # The coeff between A and R is zero.
        if not self.ignore_v4:
            raise NotImplementedError("Still d4 not implemented in this feature.")

        if verbose:
            print("L superoperator computed.")
            np.savez_compressed("L_super.npz", L_operator)
        

        def matvec(x):
            return L_operator.dot(x)
        def rmatvec(x):
            return x.dot(L_operator)

        self.L_linop = scipy.sparse.linalg.LinearOperator(L_operator.shape, matvec = matvec, rmatvec = rmatvec)
        return L_operator



            
    def run_FT(self, n_iter, save_dir = ".", verbose = True):
        """
        RUN LANCZOS ITERATIONS FOR FINITE TEMPERATURE
        =============================================

        This method performs the biconjugate Lanczos algorithm to find
        the sequence of a and b and c coefficients that are the tridiagonal representation 
        of the L matrix to be inverted.

        Parameters
        ----------
            n_iter : int
                The number of iterations to be performed in the Lanczos algorithm.
            save_dir : string
                The directory in which you want to store the results step by step,
                in order to do a preliminar analysis or restart the calculation later.
            verbose : bool
                If true all the info during the minimization will be printed on output.
        """

        # Check if the symmetries has been initialized
        if not self.initialized:
            self.prepare_symmetrization()

        # Get the current step
        i_step = len(self.a_coeffs)

        if verbose:
            header = """
<=====================================>
|                                     |
|          LANCZOS ALGORITHM          |
|                                     |
<=====================================>

Starting the algorithm. It may take a while.
Starting from step %d
""" % i_step
            print(header)

            OPTIONS = """
Should I ignore the third order effect? {}
Should I ignore the fourth order effect? {}
Max number of iterations: {}
""".format(self.ignore_v3, self.ignore_v4, n_iter)
            print(OPTIONS)


        # If this is the current step initialize the algorithm
        if i_step == 0:
            self.basis_Q = []
            self.basis_P = []
            first_vector = self.psi / np.sqrt(self.psi.dot(self.psi))
            self.basis_Q.append(first_vector)
            self.basis_P.append(first_vector)
        else:
            # Convert everything in a list
            self.basis_Q = list(self.basis_Q)
            self.basis_P = list(self.basis_P)
            self.a_coeffs = list(self.a_coeffs)
            self.b_coeffs = list(self.b_coeffs)
            self.c_coeffs = list(self.c_coeffs)
            self.arnoldi_matrix = list(self.arnoldi_matrix)

            if len(self.basis_Q) != i_step + 1:
                print("Krilov dim: %d, number of steps perfomed: %d" % (len(self.krilov_basis), i_step))
                print("Error, the krilov basis dimension should be 1 more than the number of steps")
                raise ValueError("Error the starting krilov basis does not matches the matrix, Look stdout.")

        assert len(self.basis_Q) == len(self.basis_P), "Something wrong when restoring the Lanczos."
        assert len(self.b_coeffs) == len(self.c_coeffs), "Something wrong when restoring the Lanczos. {} {}".format(len(self.b_coeffs), len(self.c_coeffs))


        # Select the two vectors for the biconjugate Lanczos iterations
        psi_q = self.basis_Q[-1]
        psi_p = self.basis_P[-1]

        print("Q basis:", self.basis_Q)
        print("P basis:", self.basis_P)
        print("SHAPE PSI Q, P :", psi_q.shape, psi_p.shape)

        next_converged = False
        for i in range(i_step, i_step+n_iter):
            if verbose:
                step_txt = """
 ===== NEW STEP %d =====

 """ % i
                print(step_txt)
                print("Length of the coefficiets: a = {}, b = {}".format(len(self.a_coeffs), len(self.b_coeffs)))
                print()

            # Apply the matrix L
            t1 = time.time()

            L_q = self.L_linop.matvec(psi_q)
            p_L = self.L_linop.rmatvec(psi_p)
            t2 = time.time()

            if verbose:
                print("Time to apply the full L: %d s" % (t2 -t1))

            # Get the a coefficient
            a_coeff = psi_p.dot(L_q)

            # Get the two residual vectors
            rk = L_q - a_coeff * psi_q 
            if len(self.basis_Q) > 1:
                rk -= self.c_coeffs[-1] * self.basis_Q[-2]

            sk = p_L - a_coeff * psi_p 
            if len(self.basis_P) > 1:
                sk -= self.b_coeffs[-1] * self.basis_P[-2]
            
            b_coeff = np.sqrt( rk.dot(rk) )
            c_coeff = sk.dot(rk) / b_coeff

            # Check the convergence
            self.a_coeffs.append(a_coeff)
            if np.abs(b_coeff) < __EPSILON__ or next_converged:
                if verbose:
                    print("Converged (b coefficient is 0)")
                converged = True
                break 
            if np.abs(c_coeff) < __EPSILON__ or next_converged:
                if verbose:
                    print("Converged (c coefficient is 0)")
                converged = True
                break


            # Get the vectors for the next iteration
            psi_q = rk / b_coeff
            psi_p = sk / c_coeff

            t1 = time.time()


            # Lets repeat twice the orthogonalization
            converged = False
            new_q = psi_q.copy()
            new_p = psi_p.copy()
            for k_orth in range(N_REP_ORTH):
                ortho_q = 0
                ortho_p = 0
                for j in range(len(self.basis_P)):
                    coeff1 = self.basis_P[j].dot(new_q) / np.sqrt(self.basis_P[j].dot(self.basis_P[j]))
                    coeff2 = self.basis_Q[j].dot(new_p)

                    # Gram Schmidt
                    new_q -= coeff1 * self.basis_P[j]
                    new_p -= coeff2 * self.basis_Q[j]

                    print("REP {} COEFF {}: scalar: {}".format(k_orth, j, coeff1))

                    ortho_q += np.abs(coeff1)
                    ortho_p += np.abs(ortho_p)

                # Add the new vector to the Krilov Basis
                normq = np.sqrt(new_q.dot(new_q))
                if verbose:
                    print("Vector norm after GS number {}: {:16.8e}".format(k_orth, normq))

                # Check the normalization (If zero the algorithm converged)
                if normq < __EPSILON__:
                    next_converged = True
                    if verbose:
                        print("Obtained a linear dependent Q vector.")
                        print("The algorithm converged.")
                    
                
                new_q /= normq

                # Normalize the p with respect to the scalar product with q
                normp = new_p.dot(new_q)

                # Check the normalization (If zero the algorithm converged)
                if normp < __EPSILON__:
                    next_converged = True
                    if verbose:
                        print("Obtained a linear dependent P vector.")
                        print("The algorithm converged.")
            
                new_p /= normp

                # We have a correctly satisfied orthogonality condition
                if ortho_p < __EPSILON__ and ortho_q < __EPSILON__:
                    break


            if not converged:
                self.basis_Q.append(new_q)
                self.basis_P.append(new_p)
                psi_q = new_q
                psi_p = new_p

                # Add the new coefficients to the Arnoldi matrix
                self.b_coeffs.append(b_coeff)
                self.c_coeffs.append(c_coeff)

            t2 = time.time()


            if verbose:
                print("Time to perform the Gram-Schmidt and retrive the coefficients: %d s" % (t2-t1))
                print()
                print("a_%d = %.8e" % (i, self.a_coeffs[-1]))
                if i > 0:
                    print("b_%d = %.8e" % (i, self.b_coeffs[-1]))
                    print("c_%d = %.8e" % (i, self.c_coeffs[-1]))
                print()
            
            # Save the step
            if not save_dir is None:
                self.save_status("%s/LANCZOS_STEP%d" % (save_dir, i))
        
                if verbose:
                    print("Status saved into '%s/LANCZOS_STEP%d'" % (save_dir, i))
            
            if verbose:
                print("Lanczos step %d ultimated." % i)
            

            if converged:
                return
            

            

        
    
    

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

    # def GetSupercellSpectralFunctionFromEig(self, w_array, smearing):
    #     r"""
    #     GET SPECTRAL FUNCTION
    #     =====================

    #     Get the spectral function from the eigenvalues and eigenvectors.
    #     The method run_full_diag must already be runned.

    #     This method returns the spectral function in the supercell.
    #     The spectral function is computed as:

    #     .. math ::

    #         G_{ab}(\omega) = \sum_{\alpha} \frac{\left<a | \lambda_\alpha\right>\left<\lambda_\alpha|b\right>}{\lambda_\alpha - \omega^2 + i\eta}

    #     where :math:`\lambda` are eigenvalues and vectors returned by this method, while :math:`\eta` is a
    #     smearing parameter chosen by the user.

    #     Parameters
    #     ----------
    #         w_array : ndarray
    #             The frequencies to which you want to compute the spectral function.
    #         smearing : float
    #             The smearing of the spectral function.

    #     Returns
    #     -------
    #         s(w) : ndarray
    #             The -ImG(w), the opposite of the imaginary part of the Green function. 
    #     """

    #     # Exclude dynamical
    #     eigvects = self.eigvects[:self.n_modes, :]

    #     N_w = len(w_array)
    #     N_alpha = len(self.eigvals)

    #     # Transform the vectors back in cartesian coordinates
    #     new_vects = np.einsum("ab, ca, c->cb", eigvects, self.pols, 1 / np.sqrt(self.m))

    #     spectral_weight = np.einsum("ab, ab -> b", new_vects, np.conj(new_vects))
    #     spectral_function = np.zeros(N_w, dtype = np.complex128)

    #     l_alpha = np.tile(self.eigvals, (N_w, 1))
    #     p_w = np.tile(spectral_weight, (N_w, 1))
    #     _w_ = np.tile(w_array, (N_alpha, 1)).T 

    #     big_mat = p_w / (l_alpha - _w_**2 + 1j*smearing)
    #     spectral_function[:] = np.sum(big_mat, axis = 1)

    #     return - np.imag(spectral_function)


    # def GetFullSelfEnergy(self):
    #     r"""
    #     GET SELF ENERGY 
    #     ===============

    #     Get the self-energy matrix from the eigenvalues and eigenvectors.
    #     The method run_full_diag must already be runned.

    #     This method returns the self energy in the supercell.
    #     It is computed as

    #     .. math ::

    #         \Pi_{ab} = \sum_{\alpha} \lambda_\alpha\left<a | \lambda_\alpha\right>\left<\lambda_\alpha|b\right>

    #     where :math:`\lambda` are eigenvalues and vectors returned by this method.
    #     The matrix is in real (cartesian) space.

    #     Returns
    #     -------
    #         s(w) : ndarray
    #             The -ImG(w), the opposite of the imaginary part of the Green function. 
    #     """

    #     # Exclude dynamical
    #     eigvects = self.eigvects[:self.n_modes, :]


    #     # Transform the vectors back in cartesian coordinates
    #     new_vects = np.einsum("ab, ca, c->cb", eigvects, self.pols, 1 / np.sqrt(self.m))

    #     self_energy = np.einsum("ab, cb, b", new_vects, np.conj(new_vects), self.eigvals)

    #     return self_energy




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

def FastApplyD3ToDyn(X, Y, rho, w, T, input_dyn,  symmetries, n_degeneracies, degenerate_space, mode = 1, transpose = False):
    """
    Apply the D3 to dyn
    ======================

    This is a wrapper to the fast C function.


    For details on the mode, look at the parameters list

    Parameters
    ----------
       X : ndarray(size = (n_modes, n_configs), dtype = np.double / np.float32)
           The X array (displacement in mode basis). Note that the dtype should match the mode
       Y : ndarray(size = (n_modes, n_configs))
           The Y array (forces in mode basis).
       rho : ndarray(size = n_configs)
           The weights of the configurations
       w : ndarray(size = n_modes)
           The list of frequencies
       T : float
           The temperature
       input_dyn : ndarray (size = n_modes*n_modes)
           The vector of the input dynamical matrix
       mode : int
           The mode for the execution:
              1) CPU : OpenMP parallelization
       symmetries : ndarray( size =(n_sym, n_modes, n_modes), dtype = np.double)
           The symmetries in the polarization basis.
       n_degeneracies : ndarray( size = n_modes, dtype = np.intc)
           The number of degenerate eigenvalues for each mode
       degenerate_space : list of lists
           The list of modes in the eigen subspace in which that mode belongs to.


    Results
    -------
       output_vector : ndarray (size = n_modes)
           The result of the calculation
    """

    n_modes = len(w)

    transp = 0
    if transpose:
        transp = 1

    output_vector = np.zeros(n_modes, dtype = TYPE_DP)
    #print( "Apply to dyn, nmodes:", n_modes, "shape:", np.shape(output_vector))
    
    deg_space_new = np.zeros(np.sum(n_degeneracies), dtype = np.intc)
    i = 0
    i_mode = 0
    j_mode = 0
    #print("len1 = ", len(deg_space_new), "len2 = ", n_modes)
    #print("Mapping degeneracies:", np.sum(n_degeneracies))
    while i_mode < n_modes:
        #print("i= ", i_mode, "Ndeg:", n_degeneracies[i_mode], "j = ", j_mode, "len = ", len(degenerate_space[i_mode]))
        #print("new_i = ", i, "tot = ", np.sum(n_degeneracies))
        #print("cross_modes: ({}, {}) | deg_imu = {} | i = {}".format(i_mode, j_mode, n_degeneracies[i_mode], i))

        deg_space_new[i] = degenerate_space[i_mode][j_mode]
        j_mode += 1
        i+=1

        if j_mode == n_degeneracies[i_mode]:
            i_mode += 1
            j_mode = 0
    

    sscha_HP_odd.ApplyV3ToDyn(X, Y, rho, w, T, input_dyn, output_vector, mode, symmetries, n_degeneracies, deg_space_new, transp)
    return output_vector


def FastApplyD3ToVector(X, Y, rho, w, T, input_vector, symmetries, n_degeneracies, degenerate_space, mode = 1):
    """
    Apply the D3 to vector
    ======================

    This is a wrapper to the fast C function.


    For details on the mode, look at the parameters list

    Parameters
    ----------
       X : ndarray(size = (n_modes, n_configs), dtype = np.double / np.float32)
           The X array (displacement in mode basis). Note that the dtype should match the mode
       Y : ndarray(size = (n_modes, n_configs))
           The Y array (forces in mode basis).
       rho : ndarray(size = n_configs)
           The weights of the configurations
       w : ndarray(size = n_modes)
           The list of frequencies
       T : float
           The temperature
       input_vector : ndarray (size = n_modes)
           The input dynamical matrix
       mode : int
           The mode for the execution:
              1) CPU : OpenMP parallelization
       symmetries : ndarray( size =(n_sym, n_modes, n_modes), dtype = np.double)
           The symmetries in the polarization basis.
       n_degeneracies : ndarray( size = n_modes, dtype = np.intc)
           The number of degenerate eigenvalues for each mode
       degenerate_space : list of lists
           The list of modes in the eigen subspace in which that mode belongs to.

    Results
    -------
       output_dyn : ndarray (size = n_modes*n_modes)
           The result of the calculation
    """
    n_modes = len(w)
    output_dyn = np.zeros(n_modes*n_modes, dtype = TYPE_DP)
    #print( "Apply to vector, nmodes:", n_modes, "shape:", np.shape(output_dyn))

    deg_space_new = np.zeros(np.sum(n_degeneracies), dtype = np.intc)
    i = 0
    i_mode = 0
    j_mode = 0
    #print("Mapping degeneracies:", np.sum(n_degeneracies))
    while i_mode < n_modes:
        #print("cross_modes: ({}, {}) | deg_i = {}".format(i_mode, j_mode, n_degeneracies[i_mode]))
        deg_space_new[i] = degenerate_space[i_mode][j_mode]
        j_mode += 1
        i += 1
        if j_mode == n_degeneracies[i_mode]:
            i_mode += 1
            j_mode = 0
    
    sscha_HP_odd.ApplyV3ToVector(X, Y, rho, w, T, input_vector, output_dyn, mode, symmetries, n_degeneracies, deg_space_new)
    return output_dyn


def FastD3_FT(X, Y, rho, w, T, input_psi, symmetries, n_degeneracies, degenerate_space, mode = 1, transpose = False):
    """
    Apply the D3 to vector
    ======================

    This is a wrapper to the fast C function.


    For details on the mode, look at the parameters list

    Parameters
    ----------
       X : ndarray(size = (n_modes, n_configs), dtype = np.double / np.float32)
           The X array (displacement in mode basis). Note that the dtype should match the mode
       Y : ndarray(size = (n_modes, n_configs))
           The Y array (forces in mode basis).
       rho : ndarray(size = n_configs)
           The weights of the configurations
       w : ndarray(size = n_modes)
           The list of frequencies
       T : float
           The temperature
       input_psi : ndarray
           The input density matrix
       mode : int
           The mode for the execution:
              1) CPU : OpenMP parallelization
       symmetries : ndarray( size =(n_sym, n_modes, n_modes), dtype = np.double)
           The symmetries in the polarization basis.
       n_degeneracies : ndarray( size = n_modes, dtype = np.intc)
           The number of degenerate eigenvalues for each mode
       degenerate_space : list of lists
           The list of modes in the eigen subspace in which that mode belongs to.

    Results
    -------
       output_psi : ndarray 
           The output density matrix
    """
    n_modes = len(w)


    transp = 0
    if transpose:
        transp = 1

    total_length = len(input_psi)

    output_psi = np.zeros(total_length, dtype = TYPE_DP)
    #print( "Apply to vector, nmodes:", n_modes, "shape:", np.shape(output_dyn))

    # Get the start and end_A
    start_A = ((n_modes + 1) * n_modes) // 2 + n_modes 
    end_A = n_modes + (n_modes + 1) * n_modes


    deg_space_new = np.zeros(np.sum(n_degeneracies), dtype = np.intc)
    i = 0
    i_mode = 0
    j_mode = 0
    #print("Mapping degeneracies:", np.sum(n_degeneracies))
    while i_mode < n_modes:
        #print("cross_modes: ({}, {}) | deg_i = {}".format(i_mode, j_mode, n_degeneracies[i_mode]))
        deg_space_new[i] = degenerate_space[i_mode][j_mode]
        j_mode += 1
        i += 1
        if j_mode == n_degeneracies[i_mode]:
            i_mode += 1
            j_mode = 0
    
    sscha_HP_odd.ApplyV3_FT(X, Y, rho, w, T, input_psi, output_psi, mode, symmetries, n_degeneracies, deg_space_new, start_A, end_A, transp)
    return output_psi



def FastD4_FT(X, Y, rho, w, T, input_psi, symmetries, n_degeneracies, degenerate_space, mode = 1):
    """
    Apply the D4 to vector
    ======================

    This is a wrapper to the fast C function.


    For details on the mode, look at the parameters list

    Parameters
    ----------
       X : ndarray(size = (n_modes, n_configs), dtype = np.double / np.float32)
           The X array (displacement in mode basis). Note that the dtype should match the mode
       Y : ndarray(size = (n_modes, n_configs))
           The Y array (forces in mode basis).
       rho : ndarray(size = n_configs)
           The weights of the configurations
       w : ndarray(size = n_modes)
           The list of frequencies
       T : float
           The temperature
       input_psi : ndarray
           The input density matrix
       mode : int
           The mode for the execution:
              1) CPU : OpenMP parallelization
       symmetries : ndarray( size =(n_sym, n_modes, n_modes), dtype = np.double)
           The symmetries in the polarization basis.
       n_degeneracies : ndarray( size = n_modes, dtype = np.intc)
           The number of degenerate eigenvalues for each mode
       degenerate_space : list of lists
           The list of modes in the eigen subspace in which that mode belongs to.

    Results
    -------
       output_psi : ndarray 
           The output density matrix
    """
    n_modes = len(w)

    total_length = len(input_psi)

    output_psi = np.zeros(total_length, dtype = TYPE_DP)
    #print( "Apply to vector, nmodes:", n_modes, "shape:", np.shape(output_dyn))

    # Get the start and end_A
    start_A = ((n_modes + 1) * n_modes) // 2 + n_modes 
    end_A = n_modes + (n_modes + 1) * n_modes


    deg_space_new = np.zeros(np.sum(n_degeneracies), dtype = np.intc)
    i = 0
    i_mode = 0
    j_mode = 0
    #print("Mapping degeneracies:", np.sum(n_degeneracies))
    # Preparing the symmetry variables for the fast calculation
    while i_mode < n_modes:
        #print("cross_modes: ({}, {}) | deg_i = {}".format(i_mode, j_mode, n_degeneracies[i_mode]))
        deg_space_new[i] = degenerate_space[i_mode][j_mode]
        j_mode += 1
        i += 1
        if j_mode == n_degeneracies[i_mode]:
            i_mode += 1
            j_mode = 0
    
    sscha_HP_odd.ApplyV4_FT(X, Y, rho, w, T, input_psi, output_psi, mode, symmetries, n_degeneracies, deg_space_new, start_A, end_A)
    return output_psi

    

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


def FastApplyD4ToDyn(X, Y, rho, w, T, input_dyn, symmetries, n_degeneracies, degenerate_space, mode = 1):
    """
    Apply the D3 to vector
    ======================

    This is a wrapper to the fast C function.

    Remember to use the correct dtype value:
    if mode == GPU:
       dtype = np.float32
    if mode == CPU:
       dtype = np.float64

    For details on the mode, look at the parameters list

    Parameters
    ----------
       X : ndarray(size = (n_modes, n_configs), dtype = np.double / np.float32)
           The X array (displacement in mode basis). Note that the dtype should match the mode
       Y : ndarray(size = (n_modes, n_configs))
           The Y array (forces in mode basis).
       rho : ndarray(size = n_configs)
           The weights of the configurations
       w : ndarray(size = n_modes)
           The list of frequencies
       T : float
           The temperature
       input_dyn : ndarray (size = n_modes*n_modes)
           The input dynamical matrix
       symmetries : ndarray( size =(n_sym, n_modes, n_modes), dtype = np.double)
           The symmetries in the polarization basis.
       n_degeneracies : ndarray( size = n_modes, dtype = np.intc)
           The number of degenerate eigenvalues for each mode
       degenerate_space : list of lists
           The list of modes in the eigen subspace in which that mode belongs to.
       mode : int
           The mode for the execution:
              1) CPU : OpenMP parallelization

    Results
    -------
       output_dyn : ndarray (size = n_modes*n_modes)
           The result of the calculation
    """
    n_modes = len(w)
    output_dyn = np.zeros(n_modes*n_modes, dtype = TYPE_DP)


    deg_space_new = np.zeros(np.sum(n_degeneracies), dtype = np.intc)
    i = 0
    i_mode = 0
    j_mode = 0
    #print("Mapping degeneracies:", np.sum(n_degeneracies))
    while i_mode < n_modes:
        #print("cross_modes: ({}, {}) | deg_i = {}".format(i_mode, j_mode, n_degeneracies[i_mode]))
        deg_space_new[i] = degenerate_space[i_mode][j_mode]
        j_mode += 1
        i += 1
        if j_mode == n_degeneracies[i_mode]:
            i_mode += 1
            j_mode = 0


    
    #print( "Apply to vector, nmodes:", n_modes, "shape:", np.shape(output_dyn))
    sscha_HP_odd.ApplyV4ToDyn(X, Y, rho, w, T, input_dyn, output_dyn, mode, symmetries, n_degeneracies, deg_space_new)
    return output_dyn





# Here some functions to analyze the data that comes out by a Lanczos
def GetFreeEnergyCurvatureFromContinuedFraction(a_ns, b_ns, pols_sc, masses, mode_mixing = True,\
    use_terminator = True, last_average = 5, smearing = 0):
    """
    GET THE FREE ENERGY CURVATURE FROM MANY LANCZOS
    ===============================================

    This function computes the free energy curvature from the result
    of a full Lanczos computation between all possible perturbations.

    Parameters
    ----------
        a_ns : ndarray(size = (n_modes, n_modes, N_steps))
            The a_n coefficients for each Lanczos perturbation
        b_ns : ndarray(size = (n_modes, n_modes, N_steps-1))
            The b_n coefficients for each Lanczos perturbation
        pols_sc : ndarray(size = (3*nat_sc, n_modes))
            The polarization vectors in the supercell
        masses : ndarray(size = (3*nat_sc))
            The mass associated to each component of pols_sc
        use_terminator : bool
            If true the infinite volume interpolation is performed trought the
            terminator trick
        last_average : int
            Used in combination with the terminator, average the last 'last_average'
            coefficients and replicate them.
        smearing : float
            The smearing for the green function calculation. 
            Usually not needed for this kind of calculation.

    Results
    -------
        odd_fc : ndarray( (3*nat_sc, 3*nat_sc))
            The free energy curvature in the supercell

    """

    n_modes = np.shape(pols_sc)[1]
    nat_sc = int(np.shape(pols_sc)[0] / 3)
    N_steps = np.shape(a_ns)[2]

    assert N_steps -1 == np.shape(b_ns)[2], "Error, an and bn has an incompatible size:\n a_n = {}, b_n = {}".format(np.shape(a_ns), np.shape(b_ns))
    
    
    mat_pol = np.zeros( (n_modes, n_modes), dtype = np.double)
    for i in range(n_modes):

        # Get the number of steps
        n_steps = np.arange(N_steps-1)[b_ns[i, i, :] == 0]
        if len(n_steps) == 0:
            n_steps = N_steps
        else:
            n_steps = n_steps[0] + 1

        
        # Create the Lanczos class
        lanc = Lanczos(None)
        lanc.a_coeffs = a_ns[i, i, :n_steps]
        lanc.b_coeffs = b_ns[i, i, :n_steps - 1]
        lanc.perturbation_modulus = 1


        print("Computing ({},{}) ... n_steps = {}".format(i, i, n_steps))

        # get the green function from continued fraction
        gf = lanc.get_green_function_continued_fraction(np.array([0]), use_terminator = use_terminator, \
            smearing = smearing, last_average = last_average)[0]
        
        mat_pol[i,i] = np.real(gf)

    # If there is the mode-mixing compute also the off-diagonal terms
    if mode_mixing:
        for i in range(n_modes):
            for j in range(i+1, n_modes):
                # Get the number of steps
                n_steps = np.arange(N_steps-1)[b_ns[i, j, :] == 0]
                if len(n_steps) == 0:
                    n_steps = N_steps
                else:
                    n_steps = n_steps[0] + 1


                # Create the Lanczos class)
                lanc = Lanczos(None)
                lanc.a_coeffs = a_ns[i, j, :n_steps]
                lanc.b_coeffs = b_ns[i, j, :n_steps-1]
                lanc.perturbation_modulus = 2

                print("Computing ({},{}) ..., n_steps = {}".format(i, j, n_steps))

                # get the green function from continued fraction
                gf = lanc.get_green_function_continued_fraction(np.array([0]), use_terminator = use_terminator, \
                    smearing = smearing, last_average = last_average)[0]
                
                # Lanczos can compute only diagonal green functions
                # Therefore we need to trick it to get the off-diagonal elements
                # <1|L|2> = 1/2*( <1+2|L|1+2> - <1|L|1>  - <2|L|2>)
                mat_pol[i,j] = (np.real(gf) - mat_pol[i,i] - mat_pol[j,j]) / 2
                mat_pol[j,i] = (np.real(gf) - mat_pol[i,i] - mat_pol[j,j]) / 2

    # The green function is the inverse of the free energy curvature
    np.savetxt("gf_mat.dat", mat_pol)
    fc_pols = np.linalg.inv(mat_pol)
    np.savetxt("fc_pols.dat", fc_pols)

    # Get back into real space
    epols_m = np.einsum("ab, a->ab", pols_sc, np.sqrt(masses)) 
    fc_odd = np.einsum("ab, ca, da ->cd", fc_pols, epols_m, epols_m)

    return fc_odd


def symmetrize_d3_muspace(d3, symmetries):
    """
    SYMMETRIZE D3 IN MODE SPACE
    ===========================

    This function symmetrizes the d3 in the mu space.
    It is quite fast.

    Parameters
    ----------
        d3 : ndarray(n_modes, n_modes, n_modes)
            The d3 tensor to be symmetrized
        symmetries : ndarray(N_sym, n_modes, n_modes)
            The full symmetry matrix

    Results
    -------
        new_d3 : ndarray(n_modes, n_modes, n_modes)
            The d3 tensor symmetrized
    """

    new_d3 = np.zeros(np.shape(d3), dtype = np.double)

    N_sym, nmode, dumb = np.shape(symmetries)

    for i in range(N_sym):
        symmat = symmetries[i, :, :]

        ap = np.einsum("abc, lc ->abl", d3, symmat)
        ap = np.einsum("abc, lb ->alc", ap, symmat)
        ap = np.einsum("abc, la ->lbc", ap, symmat)
        #ap = np.einsum("abc, aa, bb, cc->abc", d3, symmat, symmat, symmat)

        new_d3 += ap 
    
    new_d3 /= N_sym
    return new_d3