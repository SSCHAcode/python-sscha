# -*- coding: utf-8 -*-
import os
import numpy as np

import cellconstructor as CC
import cellconstructor.Structure
import cellconstructor.Phonons
import cellconstructor.Methods
import cellconstructor.Manipulate

import SCHAModules

"""
This source contains the Ensemble class
It is used to Load and Save info about the ensemble.
"""

# The small value considered zero
__EPSILON__ =  1e-6
__A_TO_BOHR__ = 1.889725989


class Ensemble:
    def __init__(self, dyn0, T0):
        """
        PREPARE THE ENSEMBLE
        ====================
        
        This method initializes and prepares the ensemble.
        
        NOTE: To now this works only in the gamma point (dyn0 must be a 1x1x1 supercell)
        
        Parameters
        ----------
            dyn0 : cellconstructor.Phonons.Phonons()
                This is the dynamical matrix used to generate the ensemble.
            T0 : float
                The temperature used to generate the ensemble.
        """
        
        # N is the number of element in the ensemble
        self.N = 0
        self.structures = []
        self.energies = []
        self.forces = []
        self.stresses = []
        
        self.sscha_energies = []
        self.sscha_forces = []
        
        # The original dynamical matrix used to generate the ensemble
        self.dyn_0 = dyn0.Copy()
        self.T0 = T0
        
        # This is the weight of each configuration in the sampling.
        # It is updated with the update_weigths function
        self.rho = []
        self.current_dyn = dyn0.Copy()
        
        self.current_T = T0
        
        
        # To avoid to recompute each time the same variables store something usefull here
        self.q_start = np.zeros( (self.N, self.dyn_0.structure.N_atoms * 3))
        self.current_q = np.zeros( (self.N, self.dyn_0.structure.N_atoms * 3))

        
    def load(self, data_dir, population, N):
        """
        LOAD THE ENSEMBLE
        =================
        
        This function load the ensemble from a standard calculation.
        
        The files need to be organized as follows
        
        data_dir / scf_populationX_Y.dat
        data_dir / energies_supercell_populationX.dat 
        data_dir / forces_populationX_Y.dat
        data_dir / pressures_populationX_Y.dat
        
        X = population
        Y = the configuration id (starting from 1 to N included, fortran convention)
        
        The files scf_population_X_Y.dat must contain the scf file of the structure.
        It should be in alat units, matching the same alat defined in the starting
        dynamical matrix.
        
        The energies_supercell.dat file must contain the total energy in Ry for
        each configuration.
        
        The forces_populationX_Y contains the 
        
        Parameters
        ----------
            data_dir : str
                The path to the directory containing the ensemble. If you used
                the fortran sscha.x code it should match the data_dir option of the
                input file.
            population : int
                The info to distinguish between several ensembles generated in the
                same data_dir. This also should match the correspective property
                of the fortran sscha.x input file.
            N : int
                The dimension of the ensemble. This should match the n_random
                variable from the fortran sscha.x input file.
        """
        A_TO_BOHR = 1.889725989
        
        # Check if the given data_dir is a real directory
        if not os.path.isdir(data_dir):
            raise IOError("Error, the given data_dir %s is not a valid directory." % data_dir)
        
        # Remove the tailoring slash if any
        if data_dir[-1] == "/":
            data_dir = data_dir[:-1]
        
        # Load the structures
        self.N = N
        
        self.forces = np.zeros( (self.N, self.dyn_0.structure.N_atoms, 3))
        self.stresses = np.zeros( (self.N, 3,3))
        
        self.sscha_energies = np.zeros(self.N)
        self.sscha_forces = np.zeros( (self.N, self.dyn_0.structure.N_atoms, 3))
        
        self.structures = []
        
        for i in range(self.N):
            # Load the structure
            structure = CC.Structure.Structure()
            structure.read_scf("%s/scf_population%d_%d.dat" % (data_dir, population, i+1), alat = self.dyn_0.alat)
            structure.has_unit_cell = True
            structure.unit_cell = self.dyn_0.structure.unit_cell
            self.structures.append(structure)
            
            # Load forces (Forces are in Ry/bohr, convert them in Ry /A)
            self.forces[i,:,:] = np.loadtxt("%s/forces_population%d_%d.dat" % (data_dir, population, i+1)) * A_TO_BOHR
            
            # Load stress
            if os.path.exists("%s/pressures_population%d_%d.dat" % (data_dir, population, i+1)):
                self.stresses[i,:,:] =  np.loadtxt("%s/pressures_population%d_%d.dat" % (data_dir, population, i+1)) 
            
            # Setup the sscha energies and forces
            energy, force = self.dyn_0.get_energy_forces(structure)
#            
#            print "Loading: config %d:" % i
#            for j in range(structure.N_atoms):
#                print "Atom %d" % j
#                print "u_disp = ", structure.get_displacement(self.dyn_0.structure)[j,:]  *A_TO_BOHR
#                print "force = ", self.forces[i, j, :] / A_TO_BOHR
#                print "SCHA force = ", force[j, :] / A_TO_BOHR
#            
#            # Debugging stuff
#            u_disp = structure.get_displacement(self.dyn_0.structure).reshape(3 * structure.N_atoms) * A_TO_BOHR
#            print "TEST DOBULE:"
#            print "NORMAL = ", self.dyn_0.dynmats[0].dot(u_disp)
#            print "INVERSE = ", self.dyn_0.dynmats[0].dot(-u_disp)
                
            
            self.sscha_energies[i] = energy 
            self.sscha_forces[i,:,:] = force
            
        # Load the energy
        self.energies = np.loadtxt("%s/energies_supercell_population%d.dat" % (data_dir, population))
        
        # Setup the initial weight
        self.rho = np.ones(self.N)
        
        # Initialize the q_start
        self.q_start = CC.Manipulate.GetQ_vectors(self.structures, self.dyn_0)
        self.current_q = self.q_start.copy()
        
    def save(self, data_dir, population):
        """
        SAVE THE ENSEMBLE
        =================
        
        
        This function saves the ensemble in a way the original fortran SSCHA code can read it.
        Look at the load function documentation to see clearely how it is saved.
        
        NOTE: This method do not save the dynamical matrix used to generate the ensemble (i.e. self.dyn_0)
        remember to save it separately to really save all the info about the ensemble.
        
        Parameters
        ----------
            data_dir : string
                Path to the directory in which the data will be saved. If it does not exists, it will be created
            population : int
                The id of the population, usefull if you want to save more ensemble in the same data_dir without overwriting
                the data.
        """
        A_TO_BOHR = 1.889725989
        
        # Check if the ensemble has really been initialized
        if self.N == 0:
            raise ValueError("Error, the ensemble seems to be not initialized.")
        
        # Check if the given data_dir is a real directory
        if os.path.exists(data_dir) and not os.path.isdir(data_dir):
            raise IOError("Error, the given data_dir %s is not a valid directory." % data_dir)
        if not os.path.exists(data_dir):
            os.mkdir(data_dir)
        
        # Remove the tailoring slash if any
        if data_dir[-1] == "/":
            data_dir = data_dir[:-1]
            
        
        # Save the energies
        np.savetxt("%s/energies_supercell_population%d.dat" % (data_dir, population), self.energies)
        
        self.dyn_0.save_qe("dyn_start_population%d_" % population)
        self.current_dyn.save_qe("dyn_end_population%d_" % population)

        
        # Check if the stresses must be saved
        save_stress = False
        if self.stresses != []:
            print self.stresses
            save_stress = True
            
        for i in range(self.N):
            # Save the forces
            np.savetxt("%s/forces_population%d_%d.dat" % (data_dir, population, i+1), self.forces[i,:,:] / A_TO_BOHR)
            
            # Save the configurations
            struct = self.structures[i]
            struct.save_scf("%s/scf_population%d_%d.dat" % (data_dir, population, i+1), self.dyn_0.alat, True)
            u_disp = struct.get_displacement(self.dyn_0.structure)
            np.savetxt("%s/u_population%d_%d.dat" % (data_dir, population, i+1), u_disp * A_TO_BOHR)
            
            # Save the stress tensors if any
            if save_stress:
                np.savetxt("%s/pressures_population%d_%d.dat" % (data_dir, population, i+1), self.stresses[i,:,:])
            
        

        
    def generate(self, N, evenodd = True):
        """
        GENERATE THE ENSEMBLE
        =====================
        
        This subroutine generates the ensemble from dyn0 and T0 setted when this
        class is created.
        You still need to generate the forces for the configurations.
        
        Parameters
        ----------
            N : int
                The number of random configurations to be extracted
            evenodd : bool, optional
                If true for each configuration also the opposite is extracted
        """
        
        if evenodd and (N % 2 != 0):
            raise ValueError("Error, evenodd allowed only with an even number of random structures")
            
        self.N = N
        self.structures = []
        if evenodd:
            structs = self.dyn_0.ExtractRandomStructures(N / 2, self.T0)
            for i, s in enumerate(structs):
                self.structures.append(s)
                new_s = s.copy()
                # Get the opposite displacement structure
                new_s.coords = self.dyn_0.structure.coords - new_s.get_displacement(self.dyn_0.structure)
                self.structures.append(new_s)
        else:
            self.structures = self.dyn_0.ExtractRandomStructures(N, self.T0)
        
        # Compute the sscha energy and forces
        self.sscha_energies = np.zeros( ( self.N))
        self.sscha_forces = np.zeros((self.N, self.dyn_0.structure.N_atoms, 3))
        self.energies = np.zeros(self.N)
        self.forces = np.zeros( (self.N, self.dyn_0.structure.N_atoms, 3))
        for i, s in enumerate(self.structures):
            energy, force  = self.dyn_0.get_energy_forces(s)
            
            self.sscha_energies[i] = energy
            self.sscha_forces[i,:,:] = force
        
        self.rho = np.ones(self.N)
        self.current_dyn = self.dyn_0.Copy()
        self.current_T = self.T0
        
        
        # Generate the q_start
        self.q_start = CC.Manipulate.GetQ_vectors(self.structures, self.dyn_0)
        self.current_q = self.q_start.copy()
        
        
    def update_weights(self, new_dynamical_matrix, newT):
        """
        IMPORTANCE SAMPLING
        ===================
        
        
        This function updates the importance sampling for the given dynamical matrix.
        The result is written in the self.rho variable
        
        
        Parameters
        ----------
            new_dynamical_matrix : CC.Phonons.Phonons()
                The new dynamical matrix on which you want to compute the averages.
            new_T : float
                The new temperature.
        """
        
        self.current_T = newT
        
        # Get the frequencies of the original dynamical matrix
        w, pols = self.dyn_0.DyagDinQ(0)
        
        # Exclude translations
        w = w[~CC.Methods.get_translations(pols, self.dyn_0.structure.get_masses_array())]

        # Convert from Ry to Ha and in fortran double precision
        w = np.array(w/2, dtype = np.float64)
        
        # Get the a_0
        old_a = SCHAModules.thermodynamic.w_to_a(w, self.T0)
        
        # Now do the same for the new dynamical matrix
        w, pols = new_dynamical_matrix.DyagDinQ(0)
        w = w[~CC.Methods.get_translations(pols, new_dynamical_matrix.structure.get_masses_array())]
        w = np.array(w/2, dtype = np.float64)
        new_a = SCHAModules.thermodynamic.w_to_a(w, newT)
        
        # Get the new q_vectors for the given matrix
        self.current_q = CC.Manipulate.GetQ_vectors(self.structures, new_dynamical_matrix) 
        
        # Convert the q vectors in the Hartree units
        old_q = self.q_start * np.sqrt(2) * __A_TO_BOHR__
        new_q = self.current_q * np.sqrt(2) * __A_TO_BOHR__
        
        
#        # Call the Fortran module to compute rho
#        print "SHAPES:"
#        print "NEW Q:", np.shape(new_q)
#        print "OLD Q:", np.shape(old_q)
#        print "NEW A:", np.shape(new_a)
#        print "OLD A:", np.shape(old_a)
        
        self.rho = SCHAModules.stochastic.get_gaussian_weight(new_q, old_q, new_a, old_a)
        
        for i in range(self.N):
            #print "Weight %d" % i
            #tmp = new_dynamical_matrix.GetRatioProbability(self.structures[i], newT, self.dyn_0, self.T0)
            #print "FORTRAN :", self.rho[i], "PYTHON:", tmp
            self.sscha_energies[i], self.sscha_forces[i, :,:] = new_dynamical_matrix.get_energy_forces(self.structures[i])
            
            
        self.current_dyn = new_dynamical_matrix
        
        
        
    def get_effective_sample_size(self):
        """
        Get the Kong-Liu effective sample size with the given importance sampling.
        """
        
        #return self.N * np.sum(self.rho) / float(np.sum(self.rho**2)) 
        return np.sum(self.rho) * np.sum(self.rho) /  np.float64(np.sum(self.rho**2))
    
    def get_average_energy(self, subtract_sscha = False, return_error = False):
        """
        GET ENERGY
        ==========
        
        This is the average of the energy
        
        .. math::
            
            \\left< E\\right> = \\frac{1}{N} \\sum_{i = 1}^{N} E_i \\rho_i
            
            
        where :math:`\\rho_i` is the ratio between the probability of extracting the configuration $i$
        with the current dynamical matrix and with the dynamical matrix used to extract the ensemble.
        
        Parameters
        ----------
            subtract_sscha : bool, optional, default False
                If true, the average difference of energy respect to the sscha one is returned. This
                is good, because you can compute analytically the sscha energy and sum it on an infinite
                ensembe. Do in this way to suppress the stochastic noise.
            return_error : bool, optional, default False
                If true also the error is returned as a second value
                
        Examples
        --------
        
        
        Example where ensemble is a correctly initialized self variable
        
        >>> energy = ensemble.get_average_energy()
        
        
        The following example return also the stochastic error
        >>> energy, error_on_energy = ensemble.get_average_energy(return_error = True)
        
        """
        
        value = 0
        value2 = 0
        
        norm = np.sum(self.rho)
        
        if subtract_sscha:
            value = np.sum( self.rho * (self.energies - self.sscha_energies)) / norm
            if return_error:
                value2 = np.sum( self.rho * (self.energies - self.sscha_energies)**2) / norm
        else:
            value = np.sum( self.rho * (self.energies)) / norm
            if return_error:
                value2 = np.sum( self.rho * (self.energies)**2) / norm
        
        # Get the error using the equation sqrt(<x^2> - <x>^2)/sqrt(N-1)
        if return_error:
            error = np.sqrt((value2 - value**2)/ (norm))
            return value, error
        return value
     
    def get_average_forces(self):
        """
        GET FORCES
        ==========
        
        This is the average of the forces that acts on the atoms
        
        .. math::
            
            \\left< \\vec F\\right> = \\frac{1}{N} \\sum_{i = 1}^{N}\\vec F_i \\rho_i
            
            
        where :math:`\\rho_i` is the ratio between the probability of extracting the configuration :math:`i`
        with the current dynamical matrix and with the dynamical matrix used to extract the ensemble.
        """
        
        return np.einsum("i, iab ->ab", self.rho, self.forces - self.sscha_forces) / np.sum(self.rho)
    
    
    
    def get_free_energy(self, return_error = False):
        """
        SSCHA FREE ENERGY
        =================
        
        Obtain the SSCHA free energy for the system.
        This is done by integrating the free energy along the hamiltonians, starting
        from current_dyn to the real system.
        
        The result is in Rydberg
        
        .. math::
            
            \\mathcal F = \\mathcal F_0 + \\int_0^1 \\frac{d\\mathcal F_\\lambda}{d\\lambda} d\\lambda
        
        Where :math:`\\lambda` is the parameter for the adiabatic integration of the hamiltonian.
        
        .. math::
            
            H(\\lambda) = H_0 + (H - H_0) \\lambda
        
        here :math:`H_0` is the sscha harmonic hamiltonian, while :math:`H_1` is the real hamiltonian 
        of the system.
        
        
        Parameters
        ----------
            return_error : bool, optional, default False
                If true also the error is returned as a second value.
        
        Returns
        -------
            float
                The free energy in the current dynamical matrix and at the ensemble temperature
        """
        K_to_Ry=6.336857346553283e-06
        
        T = self.current_T

        
        # Dyagonalize the current dynamical matrix
        nq = len(self.current_dyn.dynmats)
        
        # For each q point
        free_energy = 0
        for iq in range(nq):
            w, pols = self.current_dyn.DyagDinQ(iq)
            
            # Remove translations
            if iq == 0:
                tmask = CC.Methods.get_translations(pols, self.current_dyn.structure.get_masses_array())
                w = w[ ~tmask ]
            
            free_energy += np.sum( w / 2)
            if T > 0:
                beta = 1 / (K_to_Ry * T)
                free_energy += np.sum( 1 / beta * np.log(1 - np.exp(-beta * w)))
        
        # We got the F_0 
        # Now we can compute the free energy difference
        anharmonic_free_energy = 0
        error = 0
        if return_error:
            anharmonic_free_energy, error = self.get_average_energy(subtract_sscha = True, return_error = True)
        else:
            anharmonic_free_energy = self.get_average_energy(subtract_sscha = True, return_error = False)

        #print "Free energy harmonic:", free_energy
        #print "Free energy anharmonic:", anharmonic_free_energy
        free_energy += anharmonic_free_energy
        
        if return_error:
            return free_energy, error
        return free_energy

    
    
    def get_free_energy_gradient_respect_to_dyn(self):
        """
        FREE ENERGY GRADIENT
        ====================
        
        Get the free energy gradient respect to the dynamical matrix.
        
        .. math::
            
            \\nabla_\\Phi \\mathcal F = -\\sum_{a\\mu} \\left<\\gamma_\\mu^a q_\\mu\\right>
            
            \\gamma_\\mu^a = \\frac{e_\\mu^a \\nabla_\\Phi \\ln a_\\mu + \\nabla_\\Phi e_\mu^a}{\\sqrt M_a}(f_a - f^{\\Phi}_a)
            
            q_\\mu = \\sum_b \\sqrt M_b e_\\mu^b (R_b - \\mathcal R_b)
            
            \\nabla_\\Phi \\ln a_\\mu = \\frac{1}{2\\omega_\\mu a_\\mu} \\frac{\\partial a_\\mu}{\\partial\\omega_\\mu} \\frac{e_\\mu^a e_\\mu^b}{\\sqrt {M_aM_b}}
            
            \\nabla_\\Phi e_\mu^c  =\\sum_{\\nu \\neq \\mu} \\frac{e_\\nu^a e_\\mu^b}{\\sqrt {M_aM_b} (\\omega_\\mu^2 - \\omega_\\nu^2)} e_\\nu^c
    
    
        NOTE: it works only at gamma.
    
    
        Return
        ------
            A 3Nx3N matrix. The gradient of the free energy (To be symmetrized)
            
        """
        #K_to_Ry=6.336857346553283e-06
        
        #T = self.current_T
        # TODO: TO BE TESTED
        
        
#        # Get the mass vector
#        _m_ = np.zeros(self.dyn_0.structure.N_atoms * 3)
#        for i in range(self.current_dyn.structure.N_atoms):
#            _m_[ 3*i : 3*i + 3] = self.current_dyn.structure.masses[ self.current_dyn.structure.atoms[i]]
#        
#        _m_sqrtinv = 1 / np.sqrt(_m_)
        
        # Get the frequency and polarization vector of the dynamical matrix
        w, pols = self.current_dyn.DyagDinQ(0)
        
        
        # Discard translations and convert in Ha units
        not_trans = ~CC.Methods.get_translations(pols, self.current_dyn.structure.get_masses_array())
        w = np.array(w[not_trans] / 2, dtype = np.float64)
        pols = np.real(pols[:, not_trans])
        
        #n_modes = len(w)
        
        # Convert the q vector into Ha units
        q_new = np.array(self.current_q, dtype = np.float64) * np.sqrt(2) * __A_TO_BOHR__
        
        # Get the ityp variable 
        #ityp = self.current_dyn.structure.get_atomic_types()
        
        # Get the mass and convert in Ha units
        mass = np.array(self.current_dyn.structure.get_masses_array() * 2,
                        dtype = np.float64)
        
        nat = len(mass)
        
        # Prepare the symmetrization
        qe_sym = CC.symmetries.QE_Symmetry(self.current_dyn.structure)
        qe_sym.SetupQPoint(self.current_dyn.q_tot[0])
        
        
#        # Get the a_mu and its derivatives
#        a_mu = np.zeros(n_modes, dtype = np.float64)
#        da_dw = np.zeros(n_modes, dtype = np.float64)
        
        # Use the fortran subroutines
#        if T == 0:
#            a_mu = 1 / np.sqrt(2* w) 
#            da_dw = -1 /  np.sqrt(8 * w**3)
#        else:            
#            beta = 1 / (K_to_Ry*T)
#            a_mu = 1 / np.sqrt( np.tanh(beta*w / 2) *2* w) 
#            da_dw = - (w*beta + np.sinh(w*beta)) / (2 * np.sqrt(2) * w**2 * (np.cosh(beta*w) - 1) * np.sqrt(np.cosh(beta*w / 2) / (np.sinh(beta*w/2) * w)))
#            
#    

        # Print the sscha forces converted
        print "SCHA forces:"
        for i in range(self.N):
            for j in range(self.current_dyn.structure.N_atoms):
                print "Conf\t%d\tAtom\t%d\t" % (i, j), self.sscha_forces[i, j, :]/ (__A_TO_BOHR__)
                
                
        # Convert the forces in Ha / bohr units and in the same type as fortran
        e_forces = np.array( self.forces - self.sscha_forces, dtype = np.float64, order = "F") / (2 * __A_TO_BOHR__)
        
        # Get df_da
        df_da = SCHAModules.anharmonic.get_df_da_nonav(w, w, self.current_T, pols,
                                                       e_forces,
                                                       q_new, mass, "stat_schappp")
        #print np.shape(e_forces)
        # Now get the rest of the derivative
        
        df_dfc = np.zeros( np.shape(self.current_dyn.dynmats[0]), dtype = np.float64)
        err_df_dfc = np.zeros( np.shape(self.current_dyn.dynmats[0]), dtype = np.float64)
        
        # Just to do something good
        da_dcr_mat = np.zeros( (nat * 3, nat * 3, len(w)), dtype = np.float64)
        
        for x_i in range(self.current_dyn.structure.N_atoms * 3):
            for y_i in range(x_i, self.current_dyn.structure.N_atoms * 3):
                da_dcr, de_dcr = SCHAModules.anharmonic.get_da_dcr_and_de_dcr(w, pols, self.current_T,
                                                                              mass, x_i+1, y_i+1)
                
                print "(%d, %d): DA_DCR = " % (x_i+1, y_i+1), da_dcr
                da_dcr_mat[x_i, y_i, :] = da_dcr
                da_dcr_mat[y_i, x_i, :] = da_dcr

                df_dc, delta_df_dc = SCHAModules.anharmonic.get_df_dcoeff_av_new(df_da, da_dcr, e_forces,
                                                                                 q_new, mass, de_dcr, 
                                                                                 self.rho, 1, "err_yesrho")
                # Fill the matrix
                df_dfc[x_i, y_i] = df_dc
                df_dfc[y_i, x_i] = df_dc
                
                err_df_dfc[x_i, y_i] = delta_df_dc
                err_df_dfc[y_i, x_i] = delta_df_dc
        
        # Get the generator
        ghr = np.zeros( (3*nat, 3*nat), dtype = np.float64, order = "F")
        ghr[0,0] = 1
        # Apply the sum rule
        qe_sym.ImposeSumRule(ghr)
        # Apply symmetries
        qe_sym.SymmetrizeDynQ(ghr, self.current_dyn.q_tot[0])
        ghr /= np.sqrt(np.trace(ghr.dot(ghr)))
        print "Generator:"
        print ghr

        print "dA/dGhr = ", np.einsum("ijk, ij", da_dcr_mat, ghr)        
        
        # Force the symmetrization
        qe_sym.ImposeSumRule(df_dfc)
        qe_sym.ImposeSumRule(err_df_dfc)
        qe_sym.SymmetrizeDynQ(df_dfc, self.current_dyn.q_tot[0])
        qe_sym.SymmetrizeDynQ(err_df_dfc, self.current_dyn.q_tot[0])
        
        

#        # Prepare the w as a matrix
#        _w_ = np.tile(w, (n_modes, 1))
#        # 1 / (w_mu^2 - w_nu^2)
#        one_over_omegamunu = 1 / (_w_**2 - _w_.transpose()**2)
#        #one_over_omegamunu *= 1 - np.eye(n_modes) # Remove the therms for mu equal to nu
#        one_over_omegamunu[ (_w_ - _w_.transpose()) < __EPSILON__] = 0
#        
#        #print "freqs:", w
#        #print "w", _w_
#        #print "one_over_omega:", one_over_omegamunu
#                                        
#        # Get the derivative of the lna_mu respect to the dynamical matrix
#        # Inner product
#        d_lna_d_dyn = np.einsum("i, ai, bi, ci, a, b, c->abic", da_dw/(2 * w * a_mu), pols, pols, pols, _m_sqrtinv, _m_sqrtinv, _m_sqrtinv)
#        
#        # Get the derivative respect to the polarization vector
#        d_pol_d_dyn = np.einsum("ai,bj,ci,ji,a,b,c->abjc", pols, pols, pols, one_over_omegamunu, _m_sqrtinv, _m_sqrtinv, _m_sqrtinv)
#        
#        #print "d_lna:", d_lna_d_dyn
#        #print "d_pol:", d_pol_d_dyn
#        
#        pre_sum = d_lna_d_dyn + d_pol_d_dyn
#        
#        # Get the q vector
#        d_F_d_dyn = np.zeros(np.shape(self.current_dyn.dynmats[0]))
#        for i in range(self.N):
#            # Get the displacements of the structure
#            u_disp = self.structures[i].get_displacement(self.current_dyn.structure).reshape(3 * self.current_dyn.structure.N_atoms)
#            
#            # Get the forces on the configuration
#            delta_f = (self.forces[i,:,:] - self.sscha_forces[i,:,:]).reshape(3 * self.current_dyn.structure.N_atoms)
#            
#            # Get the q vector
#            q = np.einsum("i, ij, i", np.sqrt(_m_), pols, u_disp)
#            
#            # Get gamma matrix
#            gamma = np.einsum("abcd, d", pre_sum, delta_f)
#            
#            #print "%d) delta_f = " % (i+1), delta_f
#            #print "%d) q = " % (i+1), q
#            #print "%d) gamma = " % (i+1), gamma
#            
#            # Contract the gamma matrix and multiply it for the weight
#            partial_gradient = - np.einsum("abc, c", gamma, q)
#            d_F_d_dyn += partial_gradient * self.rho[i]
#            
#            #print "conf %d | weight %.4e | partial gradient:" % (i, self.rho[i]), partial_gradient
#            
#            
#        # Normalization
#        d_F_d_dyn /= np.sum(self.rho)
        #print "Grad:"
        #for i in range(np.shape(d_F_d_dyn)[0]):
        #    print " ".join(["%7.2e" % x for x in list(d_F_d_dyn[i,:])])
        
        #TODO: apply symmetries
            
        return df_dfc, err_df_dfc
