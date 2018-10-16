# -*- coding: utf-8 -*-
import os
import numpy as np
import time

"""
This is part of the program python-sscha
Copyright (C) 2018  Lorenzo Monacelli

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>. 
"""

import cellconstructor as CC
import cellconstructor.Structure
import cellconstructor.Phonons
import cellconstructor.Methods
import cellconstructor.Manipulate

import SCHAModules

# Try to load the parallel library if any
try:
    from mpi4py import MPI
    __MPI__ = True
except:
    __MPI__ = False

try:
    from ase.units import Rydberg, Bohr
except:
    pass

"""
This source contains the Ensemble class
It is used to Load and Save info about the ensemble.
"""

# The small value considered zero
__EPSILON__ =  1e-6
__A_TO_BOHR__ = 1.889725989


class Ensemble:
    def __init__(self, dyn0, T0, supercell = (1,1,1)):
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
            supercell: optional, list of int
                The supercell dimension
        """
        
        # N is the number of element in the ensemble
        self.N = 0
        self.structures = []
        self.energies = []
        self.forces = []
        self.stresses = []
        self.xats = []
        
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
        
        # Supercell size
        self.supercell = np.ones(3, dtype = np.intc)
        self.supercell[:] = supercell
        
        # How many atoms in the supercell
        Nsc = np.prod(self.supercell) * self.dyn_0.structure.N_atoms 
        
        # To avoid to recompute each time the same variables store something usefull here
        self.q_start = np.zeros( (self.N, Nsc * 3))
        self.current_q = np.zeros( (self.N, Nsc * 3))
        
        # Store also the displacements
        self.u_disps = np.zeros( (self.N, Nsc * 3))
        
        # A flag that memorize if the ensemble has also the stresses
        self.has_stress = False

        
    def load(self, data_dir, population, N, verbose = False):
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
            verbose : bool, optional
                If true (default false) prints the real timing of the different part
                during the loading.
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
        
        Nat_sc = np.prod(self.supercell) * self.dyn_0.structure.N_atoms
        
        self.forces = np.zeros( (self.N, Nat_sc, 3), order = "F", dtype = np.float64)
        self.xats = np.zeros( (self.N, Nat_sc, 3), order = "C", dtype = np.float64)

        self.stresses = np.zeros( (self.N, 3,3), order = "F", dtype = np.float64)
        
        self.sscha_energies = np.zeros(self.N, dtype = np.float64)
        self.energies = np.zeros(self.N, dtype = np.float64)
        self.sscha_forces = np.zeros( (self.N, Nat_sc, 3), order = "F", dtype = np.float64)
        
        self.u_disps = np.zeros( (self.N, Nat_sc * 3), order = "F", dtype = np.float64)
        
        # Add a counter to check if all the stress tensors are present
        count_stress = 0 
        
        # Superstructure
        dyn_supercell = self.dyn_0.GenerateSupercellDyn(self.supercell)
        super_structure = dyn_supercell.structure
        super_fc = dyn_supercell.dynmats[0]

        self.structures = []
        
        total_t_for_loading = 0
        total_t_for_sscha_ef = 0
        t_before_for = time.time()
        for i in range(self.N):
            # Load the structure
            structure = CC.Structure.Structure()
            if os.path.exists("%s/scf_population%d_%d.dat" % (data_dir, population, i+1)):
                t1 = time.time()
                structure.read_scf("%s/scf_population%d_%d.dat" % (data_dir, population, i+1), alat = self.dyn_0.alat)
                t2 = time.time()
                total_t_for_loading += t2 - t1
                
                structure.has_unit_cell = True
                structure.unit_cell = super_structure.unit_cell
            else:
                structure = super_structure.copy()
                t1 = time.time()
                disp =np.loadtxt("%s/u_population%d_%d.dat" % (data_dir, population, i+1)) /__A_TO_BOHR__
                t2 = time.time()
                total_t_for_loading += t2 - t1
                
                structure.coords += disp
                
            self.xats[i, :, :] = structure.coords
            self.structures.append(structure)
            
            # Get the displacement [ANGSTROM]
            self.u_disps[i,:] = structure.get_displacement(super_structure).reshape( 3 * Nat_sc)
            
            # Load forces (Forces are in Ry/bohr, convert them in Ry /A)
            t1 = time.time()
            self.forces[i,:,:] = np.loadtxt("%s/forces_population%d_%d.dat" % (data_dir, population, i+1)) * A_TO_BOHR
            
            # Load stress
            if os.path.exists("%s/pressures_population%d_%d.dat" % (data_dir, population, i+1)):
                self.stresses[i,:,:] =  np.loadtxt("%s/pressures_population%d_%d.dat" % (data_dir, population, i+1)) 
                count_stress += 1
            t2 = time.time()
            total_t_for_loading += t2 - t1
            
            # Setup the sscha energies and forces
            t1 = time.time()
            energy, force = self.dyn_0.get_energy_forces(structure, supercell = self.supercell, real_space_fc=super_fc)
            t2 = time.time()
            total_t_for_sscha_ef += t2 - t1
            
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
        
        if verbose:
            print "[LOAD ENSEMBLE]: time elapsed for the cycle over the configurations:", time.time() - t_before_for
        
        t1 = time.time()
        # Load the energy
        total_energies = np.loadtxt("%s/energies_supercell_population%d.dat" % (data_dir, population))
        t2 = time.time()
        total_t_for_sscha_ef += t2 - t1
        self.energies = total_energies[:N]
        
        # Setup the initial weight
        self.rho = np.ones(self.N, dtype = np.float64)
        
        # Initialize the q_start
        
        t1 = time.time()
        self.q_start = CC.Manipulate.GetQ_vectors(self.structures, dyn_supercell)
        t2 = time.time()
        self.current_q = self.q_start.copy()
        
        if verbose:
            print "[LOAD ENSEMBLE]: time elapsed to compute the current q vectors:", t2 - t1
            print "[LOAD ENSEMBLE]: time elapsed while loading with numpy:", total_t_for_loading
            print "[LOAD ENSEMBLE]: time elapsed for computing sscha energy and forces:", total_t_for_sscha_ef
        
        # If all the stress tensors have been found, set the stress flag
        if count_stress == self.N:
            self.has_stress = True
        else:
            self.has_stress = False
            
            # Allow the garbage collector to free the memory
            del self.stresses
        
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
            #print self.stresses
            save_stress = True
            
        super_dyn = self.dyn_0.GenerateSupercellDyn(self.supercell)
            
        for i in range(self.N):
            # Save the forces
            np.savetxt("%s/forces_population%d_%d.dat" % (data_dir, population, i+1), self.forces[i,:,:] / A_TO_BOHR)
            
            # Save the configurations
            struct = self.structures[i]
            struct.save_scf("%s/scf_population%d_%d.dat" % (data_dir, population, i+1), self.dyn_0.alat, True)
            u_disp = struct.get_displacement(super_dyn.structure)
            np.savetxt("%s/u_population%d_%d.dat" % (data_dir, population, i+1), u_disp * A_TO_BOHR)
            
            # Save the stress tensors if any
            if save_stress:
                np.savetxt("%s/pressures_population%d_%d.dat" % (data_dir, population, i+1), self.stresses[i,:,:])
            
        

    def save_bin(self, data_dir, population_id = 1):
        """
        FAST SAVE OF THE ENSEMBLE
        =========================
        
        This function is a fast way of saving the ensemble.
        It is faster and make use of less disk space than the save.
        The drawback is that can only be opened with numpy
        
        Parameters
        ----------
            data_dir : string
                path to the folder in which the ensemble is saved
            population_id : int
                The id of the population. This can be used to save
                several ensembles in the same data_dir
        """
        
        np.save("%s/energies_pop%d.npy" % (data_dir, population_id), self.energies)
        np.save("%s/forces_pop%d.npy" % (data_dir, population_id), self.forces)
        
        # Save the structures
        np.save("%s/xats_pop%d.npy" % (data_dir, population_id), self.xats)
        
        if self.has_stress:
            np.save("%s/stresses_pop%d.npy" % (data_dir, population_id), self.stresses)
        
        self.dyn_0.save_qe("%s/dyn_gen_pop%d_" % (data_dir, population_id))
        
        
        
    def load_bin(self, data_dir, population_id = 1):
        """
        LOAD THE BINARY ENSEMBLE
        ========================
        
        This function loads the ensemble saved with save_bin(...)
        
        Parameters
        ----------
            data_dir : string
                The directory containing the ensemble
            population_id : int
                The esnemble population identifier.
        """
        
        self.energies = np.load("%s/energies_pop%d.npy" % (data_dir, population_id))
        self.forces = np.load("%s/forces_pop%d.npy" % (data_dir, population_id))
        self.xats = np.load("%s/xats_pop%d.npy" % (data_dir, population_id))
        
        # Load the number of configurations
        self.N = len(self.energies)
        
        stress_path = "%s/stresses_pop%d.npy" % (data_dir, population_id)
        if os.path.exists(stress_path):
            self.stresses = np.load(stress_path)
            self.has_stress = True
        else:
            self.has_stress = False
            
        # Load the original dynamical matrix
        self.dyn_0 = CC.Phonons.Phonons("%s/dyn_gen_pop%d_" % (data_dir, population_id))
        self.current_dyn = self.dyn_0.Copy()
        
        dyn_supercell = self.dyn_0.GenerateSupercellDyn(self.supercell)
        super_structure = dyn_supercell.structure
        super_fc = dyn_supercell.dynmats[0]
        Nat_sc = super_structure.N_atoms
        
        self.sscha_energies = np.zeros(self.N, dtype = np.float64)
        self.sscha_forces = np.zeros( (self.N, Nat_sc, 3), order = "F", dtype = np.float64)
        self.u_disps = np.zeros( (self.N, 3 * Nat_sc), order = "F", dtype = np.float64)
        
        # Build the structures
        self.structures = [None] * self.N
        for i in range(self.N):
            self.structures[i] = super_structure.copy()
            self.structures[i].coords = self.xats[i,:,:]
            self.u_disps[i, :] = (self.xats[i, :, :] - super_structure.coords).reshape( 3*Nat_sc )
            
            energy, force = self.dyn_0.get_energy_forces(self.structures[i], supercell = self.supercell, 
                                                         real_space_fc=super_fc)
            
            self.sscha_energies[i] = energy
            self.sscha_forces[i, :, :] = force


        # Setup the initial weights
        self.rho = np.ones(self.N, dtype = np.float64)
        
        # Initialize the q_start
        self.q_start = CC.Manipulate.GetQ_vectors(self.structures, dyn_supercell)
        self.current_q = self.q_start.copy()

        
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
        Nat_sc = np.prod(self.supercell) * self.dyn_0.structure.N_atoms
        self.structures = []
        super_dyn = self.dyn_0.GenerateSupercellDyn(self.supercell)
        if evenodd:
            structs = super_dyn.ExtractRandomStructures(N / 2, self.T0)
            for i, s in enumerate(structs):
                self.structures.append(s)
                new_s = s.copy()
                # Get the opposite displacement structure
                new_s.coords = super_dyn.structure.coords - new_s.get_displacement(super_dyn.structure)
                self.structures.append(new_s)
        else:
            self.structures = super_dyn.ExtractRandomStructures(N, self.T0)
        
        # Compute the sscha energy and forces
        self.sscha_energies = np.zeros( ( self.N), dtype = np.float64)
        self.sscha_forces = np.zeros((self.N, Nat_sc, 3), dtype = np.float64, order = "F")
        self.energies = np.zeros(self.N, dtype = np.float64)
        self.forces = np.zeros( (self.N, Nat_sc, 3), dtype = np.float64, order = "F")
        self.stresses = np.zeros( (self.N, 3, 3), dtype = np.float64, order = "F")
        self.u_disps = np.zeros( (self.N, Nat_sc * 3), dtype = np.float64, order = "F")
        self.xats = np.zeros((self.N, Nat_sc, 3), dtype = np.float64, order = "C")
        for i, s in enumerate(self.structures):
            energy, force  = self.dyn_0.get_energy_forces(s, supercell = self.supercell, 
                                                         real_space_fc=super_dyn.dynmats[0])
            
            self.sscha_energies[i] = energy
            self.sscha_forces[i,:,:] = force
            
            # Get the displacements
            self.u_disps[i, :] = s.get_displacement(super_dyn.structure).reshape((3* Nat_sc))
            self.xats[i, :, :] = s.coords
        
        self.rho = np.ones(self.N, dtype = np.float64)
        self.current_dyn = self.dyn_0.Copy()
        self.current_T = self.T0
        
        
        # Generate the q_start
        self.q_start = CC.Manipulate.GetQ_vectors(self.structures, super_dyn)
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
        
        t1 = time.time()
        # Get the frequencies of the original dynamical matrix
        super_dyn = self.dyn_0.GenerateSupercellDyn(self.supercell)
        w, pols = super_dyn.DyagDinQ(0)
        
        # Exclude translations
        w = w[~CC.Methods.get_translations(pols, super_dyn.structure.get_masses_array())]

        # Convert from Ry to Ha and in fortran double precision
        w = np.array(w/2, dtype = np.float64)
        
        # Get the a_0
        old_a = SCHAModules.thermodynamic.w_to_a(w, self.T0)
        
        # Now do the same for the new dynamical matrix
        new_super_dyn = new_dynamical_matrix.GenerateSupercellDyn(self.supercell)
        w, pols = new_super_dyn.DyagDinQ(0)
        w = w[~CC.Methods.get_translations(pols, new_super_dyn.structure.get_masses_array())]
        w = np.array(w/2, dtype = np.float64)
        new_a = SCHAModules.thermodynamic.w_to_a(w, newT)
        
        super_structure = new_super_dyn.structure
        Nat_sc = super_structure.N_atoms
        
        self.current_q = CC.Manipulate.GetQ_vectors(self.structures, new_super_dyn) 
        
        # Convert the q vectors in the Hartree units
        old_q = self.q_start * np.sqrt(2) * __A_TO_BOHR__
        new_q = self.current_q * np.sqrt(2) * __A_TO_BOHR__
        
        t2 = time.time()
        print "Time elapsed to prepare the rho update:", t2 - t1, "s"
        
        t1 = time.time()
        self.rho = SCHAModules.stochastic.get_gaussian_weight(new_q, old_q, new_a, old_a)
        t2 = time.time()
        
        print "Time elapsed to update the stochastic weights:", t2 - t1, "s"
        #print "RHO:", self.rho
        
        for i in range(self.N):
            self.sscha_energies[i], self.sscha_forces[i, :,:] = new_super_dyn.get_energy_forces(self.structures[i], real_space_fc = new_super_dyn.dynmats[0])
            
            # Get the new displacement
            #self.u_disps[i, :] = self.structures[i].get_displacement(new_super_dyn.structure).reshape(3 * new_super_dyn.structure.N_atoms)
            self.u_disps[i, :] = (self.xats[i, :, :] - super_structure.coords).reshape( 3*Nat_sc )
        t1 = time.time()
        
        print "Time elapsed to update the sscha energies, forces and displacements:", t1 - t2, "s"
        self.current_dyn = new_dynamical_matrix.Copy()
        
        
        
    def get_effective_sample_size(self):
        """
        Get the Kong-Liu effective sample size with the given importance sampling.
        """
        
        sum_rho = np.float64(np.sum(self.rho))
        sum_rho2 = np.float64(np.sum(self.rho**2))
        return sum_rho**2 / sum_rho2
    
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
     
    def get_average_forces(self, get_error, in_unit_cell = True):
        """
        GET FORCES
        ==========
        
        This is the average of the forces that acts on the atoms
        
        .. math::
            
            \\left< \\vec F\\right> = \\frac{1}{N} \\sum_{i = 1}^{N}\\vec F_i \\rho_i
            
            
        where :math:`\\rho_i` is the ratio between the probability of extracting the configuration :math:`i`
        with the current dynamical matrix and with the dynamical matrix used to extract the ensemble.

        Parameters
        ----------
            - get_errors : bool
                If true the error is also returned (as get_free_energy).
            - in_unit_cell : bool, optional
                If True (default True) the mean force is averaged on all the atoms in the supercell,
                then it returns the forces that acts on the unit cell atoms only.
        """
        
        eforces = self.forces - self.sscha_forces
        
        if in_unit_cell and not np.prod(self.supercell) == 1:
            # Refold the forces in the unit cell
            super_structure = self.current_dyn.structure.generate_supercell(self.supercell)
            itau = super_structure.get_itau(self.current_dyn.structure) - 1 # Fort -> Py
            
            nat = self.dyn_0.structure.N_atoms
            new_forces = np.zeros((self.N, nat, 3), dtype  =np.float64, order = "C")
            
            # Project in the unit cell the forces
            for i in range(nat):
                #print "%d) ITAU LIST:" % i, itau == i
                new_forces[:, i, :] = np.sum(eforces[:, itau==i,:], axis = 1) / np.prod(self.supercell)
                #new_forces[:, i, :] = 
            
            eforces = new_forces

        force = np.einsum("i, iab ->ab", self.rho, eforces) / np.sum(self.rho)
        if get_error:
            f2 = np.einsum("i, iab ->ab", self.rho, (eforces)**2) / np.sum(self.rho)
            err = np.sqrt( f2 - force**2 )
            return force, err
        return force
    
    
    
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
                
            if len(w[w < 0]) >= 1:
                raise ValueError("Error, the dynamical matrix has imaginary frequencies")
            
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

    def get_fc_from_self_consistency(self, subtract_sscha = False, return_error = False):
        """
        SELF CONSISTENT SCHA EQUATION
        =============================
        
        This function evaluate the self consistent scha equation. This can be used
        to evaluate the goodness of the minimization procedure, as well as an
        independent minimizer.
        
        
        .. math::
            
            \\Phi_{ab} = \\frac 12 \\sum_c \\Upsilon_{ac} \\left< u_c f_a\\right>_{\\Phi}
            
        The previous equation is true only if the :math:`\\Phi` matrix is the solution
        of the SCHA theory. Here :math:`\vec u` are the displacements of the configurations
        and :math:`f` are the forces of the real system acting on the simulation.
        
        Parameters
        ----------
            subtract_sscha : bool, optional
                This is an optional parameter, if true the forces used to evaluate the 
                new force constant matrix are subtracted by the sscha forces. 
                This means that the result is a gradient of the new matrix with respect 
                to the old one.
            return_error : bool, optional
                If true also the stochastic error is returned.
                
        Results
        -------
            fc : ndarray (3*nat x 3*nat)
                The real space force constant matrix obtained by the
                self-consistent equation.
        """
        
        
        # Get the upsilon matrix
        ups_mat = self.current_dyn.GetUpsilonMatrix(self.current_T)
        
        

        # Get the pseudo-displacements obtained as
        # v = Upsilon * u = u * Upsilon^T  = u * Upsilon (we use the last to exploit fast indexing array)
        #vs = np.einsum("ij,jk", self.u_disps, ups_mat) 
        vs = self.u_disps.dot(ups_mat) # This should be faster if BLAS and MKL libraries are available (it is executed in parallel)
        
        # Build the force vector
        if subtract_sscha:
            f_vector = (self.forces - self.sscha_forces).reshape( (self.N, 3 * self.current_dyn.structure.N_atoms))
        else:
            f_vector = self.forces.reshape( (self.N, 3 * self.current_dyn.structure.N_atoms))
            
        # Average the ensemble
        new_phi = np.einsum("i, ij, ik", self.rho, vs, f_vector) / np.sum(self.rho)
        new_phi = (new_phi + np.transpose(new_phi)) * .5

        # DEBUGGING
        #np.savetxt("uf_py.dat", np.einsum("i, ij, ik", self.rho, self.u_disps, f_vector) / np.sum(self.rho), header=" <UF> matrix created by python")

        
        # Compute the stochastic error
        if (return_error):
            delta_new_phi = np.einsum("i, ij, ik", self.rho, vs**2, f_vector**2) / np.sum(self.rho) 
            delta_new_phi = (delta_new_phi + np.transpose(delta_new_phi)) * .5 
            delta_new_phi -= new_phi**2
            delta_new_phi = np.sqrt(delta_new_phi)
            return new_phi, delta_new_phi
        
        return new_phi
    
    def get_preconditioned_gradient(self, subtract_sscha = True, return_error = False, 
                                    use_ups_supercell = True, preconditioned = 1):
        """
        SELF CONSISTENT SCHA EQUATION
        =============================
        
        This function evaluate the self consistent scha equation. This can be used
        to evaluate the goodness of the minimization procedure, as well as an
        independent minimizer. This is the same as get_fc_from_self_consistency,
        but works also with supercell
        
        
        .. math::
            
            \\Phi_{ab} = \\sum_c \\upsilon_{ac} \\left< u_c f_a\\right>_{\\Phi}
            
        The previous equation is true only if the :math:`\\Phi` matrix is the solution
        of the SCHA theory. Here :math:`\\vec u` are the displacements of the configurations
        and :math:`f` are the forces of the real system acting on the simulation.
        
        NOTE: It does not takes into account for the symmetrization. 
        
        Parameters
        ----------
            subtract_sscha : bool, optional
                This is an optional parameter, if true the forces used to evaluate the 
                new force constant matrix are subtracted by the sscha forces. 
                This means that the result is a gradient of the new matrix with respect 
                to the old one.
            return_error : bool, optional
                If true also the stochastic error is returned.
            use_ups_supercell : bool, optional
                If true the gradient is computed enterely in real space, and then transformed
                with fourier in q space. This is computationally heavier, but can be used
                to test if everything is working correctly. For now this flag 
                is ignored and always True.
            preconitioned : int, optional
                If 1 (default) the gradient is returned multiplied by the preconditioned,
                otherwise it is returned as it should be.
                
        Results
        -------
            fc : ndarray (nq x 3*nat x 3*nat)
                The real space force constant matrix obtained by the
                self-consistent equation.
        """
        
        supercell_dyn = self.current_dyn.GenerateSupercellDyn(self.supercell)
        
        # Dyagonalize
        w, pols = supercell_dyn.DyagDinQ(0)
        trans = CC.Methods.get_translations(pols, supercell_dyn.structure.get_masses_array())
        ityp = supercell_dyn.structure.get_ityp() + 1 # Py to fortran convertion
        mass = np.array(supercell_dyn.structure.masses.values())
        
        log_err = "err_yesrho"
        
        mass *= 2
        w /= 2

        nat = supercell_dyn.structure.N_atoms
        eforces = np.zeros((self.N, nat, 3), dtype = np.float64, order = "F")
        u_disp = np.zeros((self.N, nat, 3), dtype = np.float64, order = "F")
        #print nat
        if subtract_sscha:
            eforces[:,:,:] = self.forces - self.sscha_forces
        else:
            eforces[:,:,:] = self.forces
        for i in range(self.N):
            u_disp[i, :, :] = np.reshape(self.u_disps[i,:], (nat, 3)) 
            

        
        grad, grad_err = SCHAModules.get_gradient_supercell(self.rho, u_disp, eforces, w, pols, trans,
                                                            self.current_T, mass, ityp, log_err, self.N,
                                                            nat, 3*nat, len(mass), preconditioned)

        # Perform the fourier transform
        q_grad = CC.Phonons.GetDynQFromFCSupercell(grad, np.array(self.current_dyn.q_tot),
                                                   self.current_dyn.structure, supercell_dyn.structure)
        q_grad_err = CC.Phonons.GetDynQFromFCSupercell(grad_err, np.array(self.current_dyn.q_tot),
                                                       self.current_dyn.structure, supercell_dyn.structure)
        
        
        
        if return_error:
            return q_grad, q_grad_err
        else:
            return q_grad
        
#        
#        nat = self.current_dyn.structure.N_atoms
#        natsc = np.prod(self.supercell) * nat
#        
#        f_vector = self.forces
#        if subtract_sscha:
#            f_vector -= self.sscha_forces
#        
#        f_vector = f_vector.reshape((self.N, 3 * natsc), order = "F")
#        
#        sum_rho = np.sum(self.rho)
#        
#        # Get the <u F> matrix
#        uf_supercell = np.einsum("i, ij, ik", self.rho, self.u_disps, f_vector) / sum_rho
#        
#        superstructure = self.dyn_0.structure.generate_supercell(self.supercell)
#        
#        # Project the <uF> matrix in q space
#        if not use_ups_supercell:
#            uf_q = CC.Phonons.GetDynQFromFCSupercell(uf_supercell, np.array(self.dyn_0.q_tot), self.dyn_0.structure, superstructure)
#        
#        if return_error:
#            uf_delta = np.einsum("i, ij, ik", self.rho, self.u_disps**2, f_vector**2) / sum_rho
#            uf_delta -= uf_supercell**2
#            if not use_ups_supercell:
#                uf_q_delta = CC.Phonons.GetDynQFromFCSupercell(uf_delta, np.array(self.dyn_0.q_tot), self.dyn_0.structure, superstructure)
#            
#        
#        
#        # For each q point, get the gradient
#        nq = len(self.dyn_0.q_tot)
#        new_phi = np.zeros( (nq, 3 * nat, 3*nat), dtype = np.complex128, order = "C")
#            
#        if return_error:
#            error_phi = np.zeros( (nq, 3 * nat, 3*nat), dtype = np.complex128, order = "C")
#        
#        if use_ups_supercell:
#            # Perform the calculation in the supercell
#            ups_mat = self.current_dyn.GenerateSupercellDyn(self.supercell).GetUpsilonMatrix(self.current_T)
#            new_phi_sc = ups_mat.dot(uf_supercell)
#            
#            # Convert in q space
#            new_phi = CC.Phonons.GetDynQFromFCSupercell(new_phi_sc, np.array(self.dyn_0.q_tot), self.dyn_0.structure, superstructure)
#            
#            if return_error:
#                error_new_phi_sc = ups_mat.dot(uf_delta)
#                error_phi = CC.Phonons.GetDynQFromFCSupercell(error_new_phi_sc, np.array(self.dyn_0.q_tot), self.dyn_0.structure, superstructure)
#        else:
#            # Perform the calculation in the q space
#            for iq in range(nq):
#                ups_mat = self.current_dyn.GetUpsilonMatrix(self.current_T, iq)
#                
#                new_phi[iq, :, :] = ups_mat.dot(uf_q[iq,:,:])
#                if return_error:
#                    error_phi[iq, :, :] = ups_mat.dot(uf_q_delta[iq,:,:])
#        
#        if return_error:
#            error_phi = np.sqrt(error_phi)
#            return new_phi, error_phi
        
 #
    
    def get_covmat_from_ensemble(self):
        """
        GET THE COVARIANCE STOCASTICALLY
        ================================
        
        This method is for testing, allows to use the ensemble to
        evaluate the covariance matrix stochastically. It should be equal
        to the matrix Upsilon^-1 that is obtained with the GetUpsilonMatrix method
        from the Phonons package.
        
        .. math::
            
            \Upsilon^{-1}_{ab} = \\left< u_a u_b\\right>
        
        Results
        -------
            cov_mat : 3nat x 3nat, ndarray
                A numpy matrix of the covariance matrix.
        """
        
        # A C style matrix of double precision real values
        cov_mat = np.einsum("i, ij, ik", self.rho, self.u_disps, self.u_disps) / np.sum(self.rho)
        
        return cov_mat
    
    
    def get_stress_tensor(self, offset_stress = None):
        """
        GET STRESS TENSOR
        =================
        
        The following subroutine computes the anharmonic stress tensor
        calling the fortran code get_stress_tensor
        
        NOTE: unit of measure is Ry/bohr^3 to match the quantum espresso one
        
        Parameters
        ----------
            offset_stress : 3x3 matrix, optional
                An offset stress to be subtracted to the real stress tensor.
                Usefull if you want to compute just the anharmonic contribution.
        
        Results
        -------
            stress_tensor : 3x3 matrix
                The anharmonic stress tensor obtained by averaging both the ab-initio
                stresses and correcting with the sscha non-linearity.
            err_stress : 3x3 matrix
                The matrix of the error on the stress tensor.
        """
        
        if not self.has_stress:
            raise ValueError("Error, the stress tensors are not present in the current ensemble.")
        
        # Volume bohr^3
        volume = np.linalg.det(self.current_dyn.structure.unit_cell) * __A_TO_BOHR__**3
            
        
        # Get frequencies and polarization vectors
        super_dyn = self.current_dyn.GenerateSupercellDyn(self.supercell)
        wr, pols = super_dyn.DyagDinQ(0)
        trans = ~ CC.Methods.get_translations(pols, super_dyn.structure.get_masses_array())
        wr = np.real( wr[trans])
        pols = np.real( pols[:, trans])
        
        nat = super_dyn.structure.N_atoms 
        
        # Get the correctly shaped polarization vectors
        er = np.zeros( (nat, len(wr), 3), dtype = np.float64, order = "F")
        
        for i in range(len(wr)):
            for j in range(nat):
                er[j, i, 0] = pols[3*j, i] 
                er[j, i, 1] = pols[3*j+1, i] 
                er[j, i, 2] = pols[3*j+2, i] 
                
        # Prepare the displacement in fortran order
        u_disps = np.zeros((self.N, nat, 3), dtype = np.float64, order = "F")
        for i in range(self.N):
            u_disps[i,:,:] = np.reshape(self.u_disps[i,:], (nat, 3))
        
        abinit_stress = np.einsum("abc -> cba", self.stresses, order = "F")
        
        stress, err_stress = SCHAModules.get_stress_tensor(volume, self.forces / __A_TO_BOHR__, u_disps * __A_TO_BOHR__, 
                                                           abinit_stress, wr, er, self.current_T, self.rho, "err_yesrho", 
                                                           self.N, nat, len(wr))
        
        
        # Check the offset
        if not offset_stress is None:
            stress -= offset_stress
        
        return stress, err_stress
    
    def get_free_energy_gradient_respect_to_dyn(self):
        """
        FREE ENERGY GRADIENT
        ====================
        
        Get the free energy gradient respect to the dynamical matrix.
        The result is in [Ry/bohr^3] as the dynamical matrix are stored
        in [Ry/bohr^2].
        
        NOTE: Not working
        
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
        
        # Convert from [Ha/bohr] in [Ry/bohr]
        df_dfc *= 2
        err_df_dfc *=  2
        

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



    def get_energy_forces(self, ase_calculator, compute_stress = True):
        """
        GET ENERGY AND FORCES FOR THE CURRENT ENSEMBLE
        ==============================================
        
        This subroutine uses the ase calculator to compute the abinitio energies and forces
        of the self ensemble.
        This subroutine requires to have ASE installed and properly configured to
        interface with your favourite ab-initio software.
        
                
        Parameters
        ----------
            ase_calculator : ase.calculator
                The ASE interface to the calculator to run the calculation.
            compute_stress : bool
                If true, the stress is requested from the ASE calculator. Be shure
                that the calculator you provide supports stress calculation
            
        """
        
        # Setup the calculator for each structure
        if __MPI__:
            comm = MPI.COMM_WORLD
            size = comm.Get_size()
            rank = comm.Get_rank()
            
            # Broad cast to all the structures
            structures = comm.bcast(self.structures, root = 0)            
            nat3 = comm.bcast(self.current_dyn.structure.N_atoms* 3, root = 0)
            N_rand = comm.bcast(self.N, root=0)
            
            # Setup the label of the calculator
            ase_calculator = comm.bcast(ase_calculator, root = 0)
            ase_calculator.set_label("esp_%d" % rank) # Avoid overwriting the same file
            
            compute_stress = comm.bcast(compute_stress, root = 0)
            
            # Check if the parallelization is correct        
            if N_rand % size != 0:
                raise ValueError("Error, for paralelization the ensemble dimension must be a multiple of the processors")
            
        else:
            size = 1
            rank = 0
            structures = self.structures
            nat3 = self.current_dyn.structure.N_atoms* 3
            N_rand = self.N
            
        # Only for the master
        
        # Prepare the energy, forces and stress array
        energies = np.zeros( N_rand / size, dtype = np.float64)
        forces = np.zeros( ( N_rand / size) * nat3 , dtype = np.float64)
        if compute_stress:
            stress = np.zeros( (N_rand / size) * 9, dtype = np.float64)

        if rank == 0:
            total_forces = np.zeros( N_rand * nat3, dtype = np.float64)
            total_stress = np.zeros( N_rand * 9, dtype = np.float64)
        else:
            total_forces = np.empty( N_rand * nat3, dtype = np.float64)
            total_stress = np.empty( N_rand * 9, dtype = np.float64)
            

        # If an MPI istance is running, split the calculation
        for i0 in range(N_rand / size):
            i = i0 + size * rank
            
            
            struct = structures[i]
            atms = struct.get_ase_atoms()
            
            # Setup the ASE calculator
            atms.set_calculator(ase_calculator)
            
            # Avoid for errors
            run = True
            count_fails = 0
            while run:
                try:
                    energy = atms.get_total_energy() / Rydberg # eV => Ry
                    run = False
                except:
                    print "Rerun the job %d" % i
                    count_fails += 1
                    if count_fails >= 5:
                        run = False
                        raise ValueError("Error in the ASE calculator for more than 5 times")
            
            # Get energy, forces (and stress)
            energy = atms.get_total_energy() / Rydberg # eV => Ry
            forces_ = atms.get_forces() / Rydberg # eV / A => Ry / A
            if compute_stress:
                stress[9*i0 : 9*i0 + 9] = atms.get_stress(False).reshape(9) * Bohr**3 / Rydberg  # ev/A^3 => Ry/bohr
            
            # Copy into the ensemble array
            energies[i0] = energy
            forces[nat3*i0 : nat3*i0 + nat3] = forces_.reshape( nat3 )
            
            # Print the status
            if rank == 0:
                print "conf %d / %d" % (i0, N_rand / size)

            
        
        # Collect all togheter
        
        if __MPI__:
            comm.Allgather([energies, MPI.DOUBLE], [self.energies, MPI.DOUBLE])
            comm.Allgather([forces, MPI.DOUBLE], [total_forces, MPI.DOUBLE])
            
            if compute_stress:
                comm.Allgather([stress, MPI.DOUBLE], [total_stress, MPI.DOUBLE])
            
        else:
            self.energies = energies
            total_forces = forces
            if compute_stress:
                total_stress = stress
        
        # Reshape the arrays
        self.forces[:, :, :] = np.reshape(total_forces, (N_rand, self.current_dyn.structure.N_atoms, 3), order = "C")
        
        if compute_stress:
            self.stresses[:,:,:] = np.reshape(total_stress, (N_rand, 3, 3), order = "C")
            self.has_stress = True
        else:
            self.has_stress = False
            
            
        
        