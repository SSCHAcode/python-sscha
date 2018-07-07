# -*- coding: utf-8 -*-
import os
import numpy as np

import cellconstructor as CC
import cellconstructor.Structure
import cellconstructor.Phonons

"""
This source contains the Ensemble class
It is used to Load and Save info about the ensemble.
"""

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
        self.dyn_0 = dyn0
        self.T0 = T0
        
        # This is the weight of each configuration in the sampling.
        # It is updated with the update_weigths function
        self.rho = []
        self.current_dyn = dyn0
        
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
        self.sscha_energies = np.zeros( (self.N, self.dyn_0.structure.N_atoms, 3))
        
        for i in range(self.N):
            # Load the structure
            structure = CC.Structure.Structure()
            structure.read_scf("%s/scf_population%d_%d.dat" % (data_dir, population, i+1), alat = self.dyn_0.alat)
            structure.has_unit_cell = True
            structure.unit_cell = self.dyn_0.structure.unit_cell
            
            # Load forces (Forces are in Ry/bohr, convert them in Ry /A)
            self.forces[i,:,:] = np.loadtxt("%s/forces_population%d_%d.dat" % (data_dir, population, i+1)) * A_TO_BOHR
            
            # Load stress
            self.stresses[i,:,:] =  np.loadtxt("%s/pressures_population%d_%d.dat" % (data_dir, population, i+1)) 
            
            # Setup the sscha energies and forces
            self.sscha_energies[i], self.sscha_forces[i,:,:] = self.dyn_0.get_energy_forces(structure)
            
        # Load the energy
        self.energies = np.loadtxt("%s/energies_supercell_population%d.dat" % (data_dir, population))
        
        # Setup the initial weight
        self.rho = np.ones(self.N)
        
        # Setup the sscha 
        
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
        
        
        for i in range(self.N):
            self.rho[i] = new_dynamical_matrix.GetRatioProbability(self.structures[i], newT, self.T0, self.dyn_0)
            self.sscha_energies[i], self.sscha_forces[i, :,:] = new_dynamical_matrix.get_energy_forces(self.structures[i])
            
        self.current_dyn = new_dynamical_matrix
        
    def get_effective_sample_size(self):
        """
        Get the Kong-Liu effective sample size with the given importance sampling.
        """
        
        return self.N * np.sum(self.rho) / float(np.sum(self.rho**2)) 
    
    def get_average_energy(self, subtract_sscha = False):
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
        """
        
        if subtract_sscha:
            return np.sum( self.rho * (self.energies - self.sscha_energies)) / self.N
        return np.sum( self.rho * (self.energies)) / self.N
     
    def get_average_forces(self):
        """
        GET FORCES
        ==========
        
        This is the average of the forces that acts on the atoms
        
        .. math::
            
            \\left< \\vec F\\right> = \\frac{1}{N} \\sum_{i = 1}^{N}\\vec F_i \\rho_i
            
            
        where :math:`\\rho_i` is the ratio between the probability of extracting the configuration $i$
        with the current dynamical matrix and with the dynamical matrix used to extract the ensemble.
        """
        
        new_rho = np.tile(self.rho, (self.dyn_0.structure.N_atoms, 3, 1))
        return np.sum( new_rho  * (self.forces - self.sscha_forces)) / self.N
    
    