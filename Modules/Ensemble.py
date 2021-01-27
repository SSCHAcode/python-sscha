# -*- coding: utf-8 -*-
from __future__ import print_function
import sys, os
import numpy as np
import time

# Fix the xrange python3 problem
try:
    xrange = xrange
    # We have Python 2
except:
    xrange = range
    # We have Python 3


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
import cellconstructor.Settings

import sscha.Parallel as Parallel


import SCHAModules
import sscha_HP_odd

_SSCHA_ODD_ = False 
#try:
#    import sscha_HP_odd
#    _SSCHA_ODD_ = True
#except:
#    _SSCHA_ODD_ = False 


# Try to load the parallel library if any
try:
    from mpi4py import MPI
    __MPI__ = True
except:
    __MPI__ = False
    
# Check if you can load spglib
try:
    import spglib
    __SPGLIB__ = True
except: 
    __SPGLIB__ = False
    
__ASE__ = True 
try:
    import ase, ase.io
except:
    __ASE__ = False

# The small value considered zero
__EPSILON__ =  1e-6
__A_TO_BOHR__ = 1.889725989

try:
    from ase.units import create_units
    units = create_units("2006")#Rydberg, Bohr
    Rydberg = units["Ry"]
    Bohr = units["Bohr"]
    __RyToK__ =  Rydberg / units["kB"]
except:
    Rydberg = 13.605698066
    Bohr = 1/__A_TO_BOHR__
    __RyToK__ = 157887.32400374097


__DEBUG_RHO__ = False

"""
This source contains the Ensemble class
It is used to Load and Save info about the ensemble.
"""


UNITS_DEFAULT = "default"
UNITS_HARTREE = "hartree"
SUPPORTED_UNITS = [UNITS_DEFAULT, UNITS_HARTREE]


class Ensemble:
    __debug_index__ = 0
    
    def __init__(self, dyn0, T0, supercell = None):
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
                The supercell dimension. If not provided, it will be determined by dyn0
        """
        
        # N is the number of element in the ensemble
        self.N = 0
        self.structures = []
        self.energies = []
        self.forces = []
        self.stresses = []
        self.xats = []
        self.units = UNITS_DEFAULT
        
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

        if supercell is not None:
            self.supercell[:] = supercell
            
            # Check if there are inconsistencies
            for i in range(3):
                if self.supercell[i] != dyn0.GetSupercell()[i]:
                    raise ValueError("""Error, you specified a supercell of {},
    while from the dynamical matrix provided I expect a supercell of {}
""".format(self.supercell, dyn0.GetSupercell()))
        else:
            self.supercell[:] = dyn0.GetSupercell()
        
        # How many atoms in the supercell
        Nsc = np.prod(self.supercell) * self.dyn_0.structure.N_atoms 
        
        # To avoid to recompute each time the same variables store something usefull here
        self.q_start = np.zeros( (self.N, Nsc * 3))
        self.current_q = np.zeros( (self.N, Nsc * 3))
        
        # Store also the displacements
        self.u_disps = np.zeros( (self.N, Nsc * 3))
        
        # A flag that memorize if the ensemble has also the stresses
        self.has_stress = True

        # A flag for each configuration that check if it possess a force and a stress
        self.force_computed = None 
        self.stress_computed = None

    def convert_units(self, new_units):
        """
        CONVERT ALL THE VARIABLE IN A COHERENT UNIT OF MEASUREMENTS
        ===========================================================

        This function is used to jump between several unit of measurement.
        You should always call this function before processing data assuming
        a particular kind of units.

        Supported units are:
           - "default" : 
               This is the default units. Here the forces are Ry/A displacements and structure are in A
               Dynamical matrix is in Ry/bohr^2. Mass is in Ry units
           - "hartree" :
               Here, everything is stored in Ha units.

        Parameters
        ----------
           - new_units : string
               The target units
        """

        # Check if the input is ok
        assert new_units in SUPPORTED_UNITS, "Error, {} unit is unknown. Try one of {}".format(new_units, SUPPORTED_UNITS)

        # If we already are in the correct units, ignore it
        if new_units == self.units:
            return
        
        if new_units == UNITS_HARTREE:
            if self.units == UNITS_DEFAULT:
                # Convert the dynamical matrix
                for iq, q in enumerate(self.dyn_0.q_tot):
                    self.dyn_0.dynmats[iq] /= 2
                    self.current_dyn.dynmats[iq] /= 2
                    self.dyn_0.q_tot[iq] /= __A_TO_BOHR__ 
                    self.current_dyn.q_tot[iq] /= __A_TO_BOHR__

                for k in self.dyn_0.structure.masses.keys():
                    self.dyn_0.structure.masses[k] *= 2
                    self.current_dyn.structure.masses[k] *= 2

                    
                # Convert the cell shape and the coordinates
                self.dyn_0.structure.coords *= __A_TO_BOHR__
                self.dyn_0.structure.unit_cell *= __A_TO_BOHR__
                self.current_dyn.structure.coords *= __A_TO_BOHR__
                self.current_dyn.structure.unit_cell *= __A_TO_BOHR__

                self.forces /= 2 * __A_TO_BOHR__ #Ry/A -> Ha/bohr
                self.sscha_forces /= 2 * __A_TO_BOHR__
                self.xats *= __A_TO_BOHR__
                self.sscha_energies /= 2 # Ry -> Ha
                self.energies /= 2
                self.u_disps *= __A_TO_BOHR__

                if self.has_stress:
                    self.stresses /= 2
            else:
                raise NotImplementedError("Error, I do not know how to convert between {} and {}.".format(self.units, new_units))    
                
        elif new_units == UNITS_DEFAULT:
            if self.units == UNITS_HARTREE:
                # Convert the dynamical matrix
                for iq, q in enumerate(self.dyn_0.q_tot):
                    self.dyn_0.dynmats[iq] *= 2
                    self.current_dyn.dynmats[iq] *= 2
                    self.dyn_0.q_tot[iq] *= __A_TO_BOHR__ 
                    self.current_dyn.q_tot[iq] *= __A_TO_BOHR__

                for k in self.dyn_0.structure.masses.keys():
                    self.dyn_0.structure.masses[k] /= 2
                    self.current_dyn.structure.masses[k] /= 2

                    
                # Convert the cell shape and the coordinates
                self.dyn_0.structure.coords /= __A_TO_BOHR__
                self.dyn_0.structure.unit_cell /= __A_TO_BOHR__
                self.current_dyn.structure.coords /= __A_TO_BOHR__
                self.current_dyn.structure.unit_cell /= __A_TO_BOHR__

                self.forces *= 2 * __A_TO_BOHR__ # Ha/bohr -> Ry/A
                self.sscha_forces *= 2 * __A_TO_BOHR__
                self.xats /= __A_TO_BOHR__
                self.sscha_energies *= 2 # Ha -> Ry
                self.energies *= 2
                self.u_disps /= __A_TO_BOHR__

                if self.has_stress:
                    self.stresses *= 2

            else:
                raise NotImplementedError("Error, I do not know how to convert between {} and {}.".format(self.units, new_units))
        else:
            raise NotImplementedError("Error, I do not know anything about this conversion")
                    

        # Update the units flag
        self.units = new_units

        
    def load(self, data_dir, population, N, verbose = False, load_displacements = True, raise_error_on_not_found = False, load_noncomputed_ensemble = False):
        """
        LOAD THE ENSEMBLE
        =================
        
        This function load the ensemble from a standard calculation.
        
        The files need to be organized as follows
        
        data_dir / scf_populationX_Y.dat
        data_dir / energies_supercell_populationX.dat 
        data_dir / forces_populationX_Y.dat
        data_dir / pressures_populationX_Y.dat
        data_dir / u_populationX_Y.dat
        
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
            load_displacement: bool
                If true the structures are loaded from the u_populationX_Y.dat files,
                otherwise they are loaded from the scf_populationX_Y.dat files.
            raise_error_on_not_found : bool
                If true, raises an error if one force file is missing
            load_noncomputed_ensemble: bool
                If True, it allows for loading an ensemble where some of the configurations forces and stresses are missing.
                Note that it must be compleated before running a SCHA minimization
        """
        A_TO_BOHR = 1.889725989
        
        # Check if the given data_dir is a real directory
        if not os.path.isdir(data_dir):
            raise IOError("Error, the given data_dir '%s' is not a valid directory." % data_dir)
        
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

        # Initialize the computation of energy and forces
        self.force_computed = np.zeros( self.N, dtype = bool)
        self.stress_computed = np.zeros(self.N, dtype = bool)
        
        # Add a counter to check if all the stress tensors are present
        count_stress = 0 
        
        # Superstructure
        dyn_supercell = self.dyn_0.GenerateSupercellDyn(self.supercell)
        super_structure = dyn_supercell.structure
        super_fc = np.real(dyn_supercell.dynmats[0])

        self.structures = []
        
        total_t_for_loading = 0
        total_t_for_sscha_ef = 0
        t_before_for = time.time()
        for i in xrange(self.N):
            # Load the structure
            structure = CC.Structure.Structure()
            if os.path.exists(os.path.join(data_dir, "scf_population%d_%d.dat" % (population, i+1))) and not load_displacements:
                t1 = time.time()
                structure.read_scf(os.path.join(data_dir, "scf_population%d_%d.dat" % (population, i+1)), alat = self.dyn_0.alat)
                t2 = time.time()
                total_t_for_loading += t2 - t1
                
                structure.has_unit_cell = True
                structure.unit_cell = super_structure.unit_cell

                # Get the displacement [ANGSTROM]
                self.u_disps[i,:] = structure.get_displacement(super_structure).reshape( 3 * Nat_sc)

            else:
                structure = super_structure.copy()
                t1 = time.time()
                disp =np.loadtxt(os.path.join(data_dir, "u_population%d_%d.dat" % (population, i+1))) /__A_TO_BOHR__
                t2 = time.time()
                total_t_for_loading += t2 - t1
                
                structure.coords += disp
                self.u_disps[i, :] = disp.ravel()
                
            self.xats[i, :, :] = structure.coords
            self.structures.append(structure)
            
            
            # Load forces (Forces are in Ry/bohr, convert them in Ry /A)
            t1 = time.time()
            force_path = os.path.join(data_dir, "forces_population%d_%d.dat" % (population, i+1))

            if os.path.exists(force_path):
                self.forces[i,:,:] = np.loadtxt(force_path) * A_TO_BOHR
                self.force_computed[i] = True
            else:
                if raise_error_on_not_found:
                    ERROR_MSG = """
Error, the file '{}' is missing from the ensemble
       data_dir = '{}'
       please, check better your data.
""".format(force_path, data_dir)
                    print(ERROR_MSG)
                    raise IOError(ERROR_MSG)
                else:
                    self.force_computed[i] = False
            
            # Load stress
            if os.path.exists(os.path.join(data_dir, "pressures_population%d_%d.dat" % (population, i+1))):
                self.stresses[i,:,:] =  np.loadtxt(os.path.join(data_dir, "pressures_population%d_%d.dat" % (population, i+1)))
                self.stress_computed[i] = True
            else:
                self.stress_computed[i] = False
            t2 = time.time()
            total_t_for_loading += t2 - t1
            
           
            
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
                
            
        
        if verbose:
            print( "[LOAD ENSEMBLE]: time elapsed for the cycle over the configurations:", time.time() - t_before_for)
        
        t1 = time.time()
        # Load the energy
        total_energies = np.loadtxt(os.path.join(data_dir, "energies_supercell_population%d.dat" % (population)))
        t2 = time.time()
        total_t_for_sscha_ef += t2 - t1
        self.energies = total_energies[:N]

        # Compute the SSCHA energies and forces
        self.sscha_energies[:], self.sscha_forces[:,:,:] = self.dyn_0.GenerateSupercellDyn(self.supercell).get_energy_forces(None, displacement = self.u_disps)

        
        # Setup the initial weight
        self.rho = np.ones(self.N, dtype = np.float64)
        
        # Initialize the q_start
        
        t1 = time.time()
        self.q_start = CC.Manipulate.GetQ_vectors(self.structures, dyn_supercell, self.u_disps)
        t2 = time.time()
        self.current_q = self.q_start.copy()
        
        if verbose:
            print( "[LOAD ENSEMBLE]: time elapsed to compute the current q vectors:", t2 - t1)
            print( "[LOAD ENSEMBLE]: time elapsed while loading with numpy:", total_t_for_loading)
            print( "[LOAD ENSEMBLE]: time elapsed for computing sscha energy and forces:", total_t_for_sscha_ef)
        

        p_count = np.sum(self.stress_computed.astype(int))
        if p_count > 0:
            self.has_stress = True
        else:
            self.has_stress = False

        # Check if the forces and stresses are present
        if not load_noncomputed_ensemble:
            if np.sum(self.force_computed.astype(int)) != self.N:
                ERROR_MSG = """
Error, the following force files are missing from the ensemble:
{}""".format(np.arange(self.N)[~self.force_computed])
                print(ERROR_MSG)
                raise IOError(ERROR_MSG)
            
            if p_count > 0 and p_count != self.N:
                ERROR_MSG = """
Error, the following stress files are missing from the ensemble:
{}""".format(np.arange(self.N)[~self.stress_computed])
                print(ERROR_MSG)
                raise IOError(ERROR_MSG)


    def load_from_calculator_output(self, directory, out_ext = ".pwo"):
        """
        LOAD THE ENSEMBLE FROM A CALCULATION
        ====================================

        This subroutine allows to directly load the ensemble from the output files
        of a calculation. This works and has been tested for quantum espresso,
        however in principle any output file from an ase supported format 
        should be readed.

        NOTE: This subroutine requires ASE to be correctly installed.

        Parameters
        ----------
            directory : string
                Path to the directory that contains the output of the calculations
            out_ext : string
                The extension of the files that will be readed.
        """

        assert __ASE__, "ASE library required to load from the calculator output file."

        # Get all the output file
        output_files = ["{}/{}".format(directory, x) for x in os.listdir(directory) if x.endswith(out_ext)]

        self.N = len(output_files)
        nat_sc = np.prod(self.supercell) * self.dyn_0.structure.N_atoms

        self.forces = np.zeros( (self.N, nat_sc, 3), order = "F", dtype = np.float64)
        self.xats = np.zeros( (self.N, nat_sc, 3), order = "C", dtype = np.float64)

        self.stresses = np.zeros( (self.N, 3,3), order = "F", dtype = np.float64)
        
        self.sscha_energies = np.zeros(self.N, dtype = np.float64)
        self.energies = np.zeros(self.N, dtype = np.float64)
        self.sscha_forces = np.zeros( (self.N, nat_sc, 3), order = "F", dtype = np.float64)
        
        self.u_disps = np.zeros( (self.N, nat_sc * 3), order = "F", dtype = np.float64)
        
        # Add a counter to check if all the stress tensors are present
        count_stress = 0 
        
        # Superstructure
        dyn_supercell = self.dyn_0.GenerateSupercellDyn(self.supercell)
        super_structure = dyn_supercell.structure
        super_fc = dyn_supercell.dynmats[0]

        self.structures = []

        for i, outf in enumerate(output_files):
            ase_struct = ase.io.read(outf)

            # Get the structure
            structure = CC.Structure.Structure()
            structure.generate_from_ase_atoms(ase_struct)

            self.xats[i, :, :] = structure.coords
            self.structures.append(structure)

            # Get the displacement [ANGSTROM]
            self.u_disps[i,:] = structure.get_displacement(super_structure).reshape( 3 * nat_sc)

            # Get the energy
            energy = ase_struct.get_potential_energy()
            energy /= Rydberg
            self.energies[i] = energy

            # Get the forces [eV/A -> Ry/A]
            forces = ase_struct.get_forces() / Rydberg 
            self.forces[i, :, :] = forces

            # Get the stress if any
            try:
                stress = ase_struct.get_stress(voigt=False)
                # eV/A^3 -> Ry/bohr^3
                stress /= Rydberg / Bohr**3
                count_stress += 1
                self.stresses[i, :, :] = stress
            except:
                pass

            # Get the sscha energy and forces            
            energy, force = self.dyn_0.get_energy_forces(structure, supercell = self.supercell, real_space_fc=super_fc)

            self.sscha_energies[i] = energy 
            self.sscha_forces[i,:,:] = force

        self.rho = np.ones(self.N, dtype = np.float64)

        t1 = time.time()
        self.q_start = CC.Manipulate.GetQ_vectors(self.structures, dyn_supercell, self.u_disps)
        t2 = time.time()
        self.current_q = self.q_start.copy()

        if count_stress == self.N:
            self.has_stress = True
        else:
            self.has_stress = False
            
            

        
    def save(self, data_dir, population, use_alat = False):
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
            use_alat : bool
                If true the scf_populationX_Y.dat files will be saved in alat units, as specified by the dynamical matrix.
                Also the unit cell will be omitted. This is done to preserve retrocompatibility with ensembles generated by
                older versions of the sscha code
        """
        A_TO_BOHR = 1.889725989

        # Check if the data dir exists
        if not os.path.exists(data_dir):
            os.makedirs(data_dir)
        
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
        np.savetxt(os.path.join(data_dir, "energies_supercell_population%d.dat" % (population)), self.energies)
        
        self.dyn_0.save_qe("dyn_start_population%d_" % population)
        self.current_dyn.save_qe("dyn_end_population%d_" % population)

        
        super_dyn = self.dyn_0.GenerateSupercellDyn(self.supercell)
            
        for i in xrange(self.N):
            # Save the forces
            if self.force_computed[i]:
                np.savetxt("%s/forces_population%d_%d.dat" % (data_dir, population, i+1), self.forces[i,:,:] / A_TO_BOHR)
            
            # Save the configurations
            struct = self.structures[i]
            if use_alat:
                struct.save_scf("%s/scf_population%d_%d.dat" % (data_dir, population, i+1), self.dyn_0.alat, True)
            else:
                struct.save_scf("%s/scf_population%d_%d.dat" % (data_dir, population, i+1))

            u_disp = struct.get_displacement(super_dyn.structure)
            np.savetxt("%s/u_population%d_%d.dat" % (data_dir, population, i+1), u_disp * A_TO_BOHR)
            
            # Save the stress tensors if any
            if self.has_stress and self.stress_computed[i]:
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

        # Check if the data dir exists
        if not os.path.exists(data_dir):
            os.makedirs(data_dir)
        

        if Parallel.am_i_the_master():
            np.save("%s/energies_pop%d.npy" % (data_dir, population_id), self.energies)
            np.save("%s/forces_pop%d.npy" % (data_dir, population_id), self.forces)
            
            # Save the structures
            np.save("%s/xats_pop%d.npy" % (data_dir, population_id), self.xats)
            
            if self.has_stress:
                np.save("%s/stresses_pop%d.npy" % (data_dir, population_id), self.stresses)
            
            self.dyn_0.save_qe("%s/dyn_gen_pop%d_" % (data_dir, population_id))
        
    def load_bin(self, data_dir, population_id = 1, avoid_loading_dyn = False):
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
            avoid_loading_dyn : bool
                If true, the dynamical matrix is not loaded.
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
        if not avoid_loading_dyn:
            self.dyn_0 = CC.Phonons.Phonons("%s/dyn_gen_pop%d_" % (data_dir, population_id), self.dyn_0.nqirr)
            self.current_dyn = self.dyn_0.Copy()
            self.supercell = self.dyn_0.GetSupercell()
            
        super_structure = self.dyn_0.structure.generate_supercell(self.supercell)
        Nat_sc = super_structure.N_atoms
        
        self.sscha_energies = np.zeros(self.N, dtype = np.float64)
        self.sscha_forces = np.zeros( (self.N, Nat_sc, 3), order = "F", dtype = np.float64)
        self.u_disps = np.zeros( (self.N, 3 * Nat_sc), order = "F", dtype = np.float64)
        
        # Build the structures
        self.structures = [None] * self.N
        for i in xrange(self.N):
            self.structures[i] = super_structure.copy()
            self.structures[i].coords = self.xats[i,:,:]
            self.u_disps[i, :] = (self.xats[i, :, :] - super_structure.coords).reshape( 3*Nat_sc )
            
            
        self.sscha_energies[:], self.sscha_forces[:,:,:] = self.dyn_0.get_energy_forces(None, displacement = self.u_disps)

        # Setup the initial weights
        self.rho = np.ones(self.N, dtype = np.float64)
        

    def init_from_structures(self, structures):
        """
        Initialize the ensemble from the given list of structures

        Parameters
        ----------
            structures : list of structures
                The list of structures used to initialize the ensemble
        """

        # Perform the standard initialization

        self.N = len(structures)
        Nat_sc = np.prod(self.supercell) * self.dyn_0.structure.N_atoms

        self.structures = [x for x in structures]

        self.sscha_energies = np.zeros( ( self.N), dtype = np.float64)
        self.sscha_forces = np.zeros((self.N, Nat_sc, 3), dtype = np.float64, order = "F")
        
        self.energies = np.zeros(self.N, dtype = np.float64)
        self.forces = np.zeros( (self.N, Nat_sc, 3), dtype = np.float64, order = "F")
        self.stresses = np.zeros( (self.N, 3, 3), dtype = np.float64, order = "F")
        self.u_disps = np.zeros( (self.N, Nat_sc * 3), dtype = np.float64, order = "F")
        self.xats = np.zeros((self.N, Nat_sc, 3), dtype = np.float64, order = "C")
        for i, s in enumerate(self.structures):
            # Get the displacements
            self.xats[i, :, :] = s.coords
        
        # TODO:
        # Here it is useless to generate the supercell dynamical matrix, 
        # it should be replaced by generating the unit cell structure, 
        # But then the get_energy_forces method should provide the correct implementation.
        new_super_dyn = self.current_dyn.GenerateSupercellDyn(self.current_dyn.GetSupercell())
        self.u_disps[:,:] = np.reshape(self.xats - np.tile(new_super_dyn.structure.coords, (self.N, 1,1)), (self.N, 3 * Nat_sc), order = "C") 

        self.sscha_energies[:], self.sscha_forces[:,:,:] = self.dyn_0.get_energy_forces(None, displacement = self.u_disps)

        
        self.rho = np.ones(self.N, dtype = np.float64)
        self.current_dyn = self.dyn_0.Copy()
        self.current_T = self.T0

        # Setup that both forces and stresses are not computed
        self.stress_computed = np.zeros(self.N, dtype = bool)
        self.force_computed = np.zeros(self.N, dtype = bool)


    def generate(self, N, evenodd = True, project_on_modes = None):
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
            project_on_modes : ndarray(size=(3*nat_sc, nproj)), optional
                If different from None the displacements are projected on the
                given modes.
        """
        
        if evenodd and (N % 2 != 0):
            raise ValueError("Error, evenodd allowed only with an even number of random structures")
            
        self.N = N
        Nat_sc = np.prod(self.supercell) * self.dyn_0.structure.N_atoms
        self.structures = []
        #super_dyn = self.dyn_0.GenerateSupercellDyn(self.supercell)
        super_struct = self.dyn_0.structure.generate_supercell(self.dyn_0.GetSupercell())

        structures = []
        if evenodd:
            structs = self.dyn_0.ExtractRandomStructures(N // 2, self.T0, project_on_vectors = project_on_modes)



            for i, s in enumerate(structs):
                structures.append(s)
                new_s = s.copy()
                # Get the opposite displacement structure
                new_s.coords = super_struct.coords - new_s.get_displacement(super_struct)
                structures.append(new_s)
        else:
            structures = self.dyn_0.ExtractRandomStructures(N, self.T0, project_on_vectors = project_on_modes)

        self.init_from_structures(structures)
        
    def get_unwrapped_ensemble(self, subtract_sscha = True, verbose = True):
        """
        This subroutine gets the displacements, forces and stochastic weights of the ensemble
        unwrapping with the symmetries: for each configuraiton, we add all other configurations equivalent by simmetries
        
        NOTE: it works only if SPGLIB is installed

        Parameter
        ---------
            subtract_sscha : bool
                If true (default), instead of the forces, the method returns the forces subtracted by the
                SCHA forces.
            verbose : bool
                If true, the method prints into stdout the timing.

        Returns
        -------
            u_disps : ndarray(size = (n_configs * n_syms, 3*nat), dtype = np.double)
                The displacements of atomic configurations with respect to the average positions
            forces : ndarray(size = (n_configs * n_syms, 3*nat), dtype = np.double)
                The forces that acts on each configuration (subtracted by the SSCHA if requested)
            weights : ndarray(size = n_configs * n_syms, dytpe = no.double)
                The weights of the configurations
        """

        # First of all, we need to get the symmetries

        # Get the symmetries
        if not __SPGLIB__:
            raise ImportError("Error, get_unwrapped_ensemble mehtod requires spglib")

        # Get the symmetries from spglib
        super_structure = self.current_dyn.structure.generate_supercell(self.supercell)
        spglib_syms = spglib.get_symmetry(super_structure.get_ase_atoms())
        
        # Convert them into the cellconstructor format
        cc_syms = CC.symmetries.GetSymmetriesFromSPGLIB(spglib_syms, False)
        
        print("N syms:", len(cc_syms))
        
        n_syms = len(cc_syms)
        nat_sc = super_structure.N_atoms
        new_N = n_syms * self.N
        
        # Get the IRT atoms
        irts = np.zeros( (n_syms, nat_sc), dtype = int)
        for i in range(n_syms):
            irts[i, :] = CC.symmetries.GetIRT(super_structure, cc_syms[i]) + 1 # Py -> Fortran indexing



        old_udisps = np.zeros( self.u_disps.shape, dtype = np.double)
        old_forces = np.zeros( self.forces.shape, dtype = np.double)
        new_udisps = np.zeros( (new_N, 3 * nat_sc), dtype = np.double)
        new_forces = np.zeros( (new_N, 3 * nat_sc), dtype = np.double)

        # Convert to crystal coordinates
        t1 = time.time()
        for i in range(self.N):
            v = self.u_disps[i, :].reshape((nat_sc, 3))
            old_udisps[i, :] = CC.Methods.cart_to_cryst(super_structure.unit_cell, v).ravel()


            v = self.forces[i, :].reshape((nat_sc, 3))
            if subtract_sscha:
                v -= self.sscha_forces[i, :].reshape((nat_sc, 3))

            old_forces[i, :] = CC.Methods.cart_to_cryst(super_structure.unit_cell, v).ravel()
        t2 = time.time()

        if verbose:
            print("Time to convert everything to crystal coordinates: {} s".format(t2 - t1))

        # Unwrap the ensemble
        new_udisps[:,:] = SCHAMethods.unwrap_ensemble(old_udisps, cc_syms[:3, :3].astype(int), irts, nat_sc, n_syms)
        new_forces[:,:] = SCHAMethods.unwrap_ensemble(old_forces, cc_syms[:3, :3].astype(int), irts, nat_sc, n_syms)


        t3 = time.time()
        if verbose:
            print("Time to unwrap the ensemble: {} s".format(t3 - t2))

        # Convert to cartesian coordinates once more
        v = new_udisps.reshape((new_N, nat_sc, 3))
        new_udisps = CC.Methods.cryst_to_cart(super_structure.unit_cell, v).reshape((new_N, 3 * nat_sc))

        v = new_forces.reshape((new_N, nat_sc, 3))
        new_forces = CC.Methods.cryst_to_cart(super_structure.unit_cell, v).reshape((new_N, 3*nat_sc))

        t4 = time.time()
        if verbose:
            print("Time to convert back to cartesian: {} s".format(t4 - t3))

        weights = np.zeros( (new_N), dtype = np.double)
        for i in range(self.N):
            weights[n_syms * i : n_syms * (i + 1)] = self.rho[i]
        t5 = time.time()

        if verbose:
            print("Time to unwrap the weights: {} s".format(t5 - t4))

            print("    overall time of get_unwrapped_ensemble: {} s".format(t5- t1))

        return new_udisps, new_forces, weights



        
    def _unwrap_symmetries_(self):
        """
        UNWRAP THE ENSEMBLE
        ===================
        
        This function unwraps the current ensemble by introducing the displacements
        (and energies and forces) of the symmetric specular configurations.
        This allows for a simple simmetrization of the odd3 correction.
        
        NOTE: stress tensors are not unwrapped!
        
        NOTE: This works only if spglib is installed
        """
        
        # Get the symmetries
        if not __SPGLIB__:
            raise ImportError("Error, unwrap_symmetries mehtod requires spglib to be importable")
            
        # Get the symmetries from spglib
        super_structure = self.current_dyn.structure.generate_supercell(self.supercell)
        spglib_syms = spglib.get_symmetry(super_structure.get_ase_atoms())
        
        # Convert them into the cellconstructor format
        cc_syms = CC.symmetries.GetSymmetriesFromSPGLIB(spglib_syms, False)
        
        print("N syms:", len(cc_syms))
        
        n_syms = len(cc_syms)
        nat_sc = super_structure.N_atoms
        
        # Get the IRT atoms
        t1 = time.time()
        irts = np.zeros( (n_syms, nat_sc), dtype = int)
        for i in range(n_syms):
            irts[i, :] = CC.symmetries.GetIRT(super_structure, cc_syms[i])
        t2 = time.time()
        
        print ("Time elapsed to compute IRTS:", t2 - t1, "s")
        
        new_N = self.N * n_syms
        u_disps_new = np.zeros( (new_N, 3 * nat_sc), dtype = np.float64, order = "F")
        forces_new = np.zeros( (new_N, nat_sc, 3), dtype = np.float64, order = "F")
        rho_new = np.zeros( (new_N), dtype = np.float64)
        xats_new = np.zeros( (new_N, nat_sc, 3), dtype = np.float64, order = "C")
        energies_new = np.zeros( (new_N), dtype = np.float64)
        sscha_energies_new = np.zeros( (new_N), dtype = np.float64)
        sscha_forces_new = np.zeros( (new_N, nat_sc, 3), dtype = np.float64, order = "F")
        new_structures = []
        
        t1 = time.time()
        for x in range(self.N):
            print ("Config %d" % x)
            for i in range(n_syms):

                index = n_syms * x + i
                
                u_v = self.u_disps[x, :].reshape((nat_sc, 3))
                u_new = CC.symmetries.ApplySymmetryToVector(cc_syms[i], u_v, super_structure.unit_cell, irts[i, :])
                u_disps_new[index, :] = u_new.ravel()
                xats_new[index, :, :] = super_structure.coords + u_new
                
                f_v = self.forces[x, :, :]
                f_new = CC.symmetries.ApplySymmetryToVector(cc_syms[i], f_v, super_structure.unit_cell, irts[i, :])
                forces_new[index, :, :] = f_new
                
                f_v = self.sscha_forces[x, :, :]                
                f_new = CC.symmetries.ApplySymmetryToVector(cc_syms[i], f_v, super_structure.unit_cell, irts[i, :])
                sscha_forces_new[index, :, :] = f_new
                
                rho_new[index] = self.rho[x]
                energies_new[index] = self.energies[x]
                sscha_energies_new[index] = self.sscha_energies[x]
                
                tmp_struct = super_structure.copy()
                tmp_struct.coords = xats_new[index, :, :]
                new_structures.append(tmp_struct)
        t2 = time.time()
        
        print ("Time elapsed unwrapping the ensemble:", t2 - t1, "s")
        
        # Update the ensemble
        self.N = new_N
        self.u_disps = u_disps_new
        self.xats = xats_new
        self.rho = rho_new
        self.structures = new_structures
        self.energies = energies_new
        self.forces = forces_new
        self.sscha_energies = sscha_energies_new
        self.sscha_forces = sscha_forces_new
                
        
        
    def update_weights(self, new_dynamical_matrix, newT, update_q = False):
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
            update_q : bool
                If false the q_vectors are not updated. This is required for some
                methods and application, but not for standard minimization.
                Since it is the most time consuming part, it can be safely avoided. 
        """
        
        self.current_T = newT

            
        
        t1 = time.time()
        # Get the frequencies of the original dynamical matrix
        super_dyn = self.dyn_0.GenerateSupercellDyn(self.supercell)

        w, pols = self.dyn_0.DiagonalizeSupercell()#super_dyn.DyagDinQ(0)
        
        # Exclude translations
        trans_original = CC.Methods.get_translations(pols, super_dyn.structure.get_masses_array()) 
        w = w[~trans_original]

        # Convert from Ry to Ha and in fortran double precision
        w = np.array(w/2, dtype = np.float64)
        
        # Get the a_0
        old_a = SCHAModules.thermodynamic.w_to_a(w, self.T0)
        
        # Now do the same for the new dynamical matrix
        new_super_dyn = new_dynamical_matrix.GenerateSupercellDyn(self.supercell)
        
        w, pols = new_dynamical_matrix.DiagonalizeSupercell()#new_super_dyn.DyagDinQ(0)

        trans_mask = CC.Methods.get_translations(pols, new_super_dyn.structure.get_masses_array()) 

        # Check if the new dynamical matrix satisfies the sum rule
        if (np.sum(trans_mask.astype(int)) != 3) or (np.sum(trans_original.astype(int)) != 3):
            ERR_MSG = """
ERROR WHILE UPDATING THE WEIGHTS
    
Error, one dynamical matrix does not satisfy the acoustic sum rule.
       If this problem arises on a sscha run, 
       it may be due to a gradient that violates the sum rule.
       Please, be sure you are not using a custom gradient function.

DETAILS OF ERROR:
    Number of translatinal modes in the original dyn = {}
    Number of translational modes in the target dyn = {}
    (They should be both 3)
""".format(np.sum(trans_original.astype(int)), np.sum(trans_mask.astype(int)))

            print(ERR_MSG)
            raise ValueError(ERR_MSG)

        w = w[~trans_mask]
        w = np.array(w/2, dtype = np.float64)
        new_a = SCHAModules.thermodynamic.w_to_a(w, newT)
        
        super_structure = new_super_dyn.structure
        Nat_sc = super_structure.N_atoms
        
        # Get the new displacements in the supercell
        t3 = time.time()
        old_disps = np.zeros(np.shape(self.u_disps), dtype = np.double)

        self.u_disps[:,:] = self.xats.reshape((self.N, 3*Nat_sc)) - np.tile(super_structure.coords.ravel(), (self.N,1))
        old_disps[:,:] = self.xats.reshape((self.N, 3*Nat_sc)) - np.tile(super_dyn.structure.coords.ravel(), (self.N,1))

        # for i in xrange(self.N):
        #     self.u_disps[i, :] = (self.xats[i, :, :] - super_structure.coords).reshape( 3*Nat_sc )
            
        #     old_disps[i,:] = (self.xats[i, :, :] - super_dyn.structure.coords).reshape( 3*Nat_sc )
            
        #     # TODO: this method recomputes the displacements, it is useless since we already have them in self.u_disps
        self.sscha_energies[:], self.sscha_forces[:,:,:] = new_super_dyn.get_energy_forces(None, displacement = self.u_disps)

        t4 = time.time()

        if update_q:
            self.current_q = CC.Manipulate.GetQ_vectors(self.structures, new_super_dyn, self.u_disps) 
        
        # Convert the q vectors in the Hartree units
        #old_q = self.q_start * np.sqrt(np.float64(2)) * __A_TO_BOHR__
        #new_q = self.current_q * np.sqrt(np.float64(2)) * __A_TO_BOHR__
        
        
        #t1 = time.time()
        #self.rho = SCHAModules.stochastic.get_gaussian_weight(new_q, old_q, new_a, old_a)
        #t2 = time.time()
        
        if __DEBUG_RHO__:
            print( " ==== [UPDATE RHO DEBUGGING] ==== ")
            print( " INPUT INFO: ")
            np.savetxt("rho_%05d.dat" % self.__debug_index__, self.rho)
            print( " rho saved in ", "rho_%05d.dat" % self.__debug_index__)
        
        
        # Get the covariance matrices of the ensemble
        ups_new = np.real(new_super_dyn.GetUpsilonMatrix(self.current_T))
        ups_old = np.real(super_dyn.GetUpsilonMatrix(self.T0))

        # Get the normalization ratio
        #norm = np.sqrt(np.abs(np.linalg.det(ups_new) / np.linalg.det(ups_old))) 
        norm = np.prod( old_a / new_a)

        t2 = time.time()
        print( "Time elapsed to prepare the rho update:", t2 - t1, "s")
        print ( "(of which to update sscha energies and forces: %.4f s)" % (t4-t3))
        print ( "(of which computing the Upsilon matrix: %.4f s)" % (t2 - t4))
        
        rho_tmp = np.ones( self.N, dtype = np.float64) * norm 
        if __DEBUG_RHO__:
            print("Norm factor:", norm)
        for i in xrange(self.N):
            v_new = self.u_disps[i, :].dot(ups_new.dot(self.u_disps[i, :])) * __A_TO_BOHR__**2
            v_old = old_disps[i, :].dot(ups_old.dot(old_disps[i, :])) * __A_TO_BOHR__**2

            if __DEBUG_RHO__:
                print("CONF {} | displacement = {}".format(i, v_new - v_old))
            rho_tmp[i] *= np.exp(-0.5 * (v_new - v_old) )
        # Lets try to use this one
        self.rho = rho_tmp

        #print("\n".join(["%8d) %16.8f" % (i+1, r) for i, r in enumerate(self.rho)]))
        
        #np.savetxt("upsilon_%05d.dat" % self.__debug_index__, ups_new)
        #np.savetxt("d_upsilon_%05d.dat" % self.__debug_index__, dups)

        
        #print "RHO:", self.rho
        
        #for i in range(self.N):
            # Get the new displacement
            #self.u_disps[i, :] = self.structures[i].get_displacement(new_super_dyn.structure).reshape(3 * new_super_dyn.structure.N_atoms)
            #self.u_disps[i, :] = (self.xats[i, :, :] - super_structure.coords).reshape( 3*Nat_sc )
        t1 = time.time()
        print( "Time elapsed to update weights the sscha energies, forces and displacements:", t1 - t3, "s")
        print( "(of which to update the weights):", t1 - t2, "s")
        self.current_dyn = new_dynamical_matrix.Copy()
        
        
        if __DEBUG_RHO__:
            new_dynamical_matrix.save_qe("ud_%05d" % self.__debug_index__)
            print( " new_dynmat saved in ud_%05d " % self.__debug_index__)
            print( " new_T : ", newT)
            print( " old_T : ", self.T0)
            print( " supercell :", self.supercell)
            self.dyn_0.save_qe("sd_%05d" % self.__debug_index__)
            print( " starting dyn saved in sd_%05d" % self.__debug_index__)
            print( " old_a:", " ".join(["%16.8f" %  x for x in old_a]))
            print( " new_a:", " ".join(["%16.8f" %  x for x in new_a]))
            #np.savetxt("old_q_%05d.dat" %self.__debug_index__, old_q)
            print( " old_q saved in ", "old_q_%05d.dat" %self.__debug_index__)
            #np.savetxt("new_q_%05d.dat" %self.__debug_index__, new_q)
            print( " new_q saved in ", "new_q_%05d.dat" %self.__debug_index__)
            print( " u_disps saved in ", "u_disps_%05.dat" % self.__debug_index__)
            np.savetxt("u_disps_%05d.dat" % self.__debug_index__, self.u_disps)
            
            print( " The last rho in", "rho_last_%05d.dat" % self.__debug_index__)
            np.savetxt("rho_last_%05d.dat" % self.__debug_index__, self.rho)
            print( " The other rho kind saved in", "other_rho_kind_%05d.dat"  % self.__debug_index__)
            np.savetxt("other_rho_kind_%05d.dat"  % self.__debug_index__, rho_tmp)
            print( " The KL according to other rho kind:", np.sum(rho_tmp)**2 / np.sum(rho_tmp**2))
            self.__debug_index__ += 1
            
            
        
        
        
    def get_effective_sample_size(self):
        """
        Get the Kong-Liu effective sample size with the given importance sampling.
        """
        
        sum_rho = np.float64(np.sum(self.rho))
        sum_rho2 = np.float64(np.sum(self.rho**2))
        kl = sum_rho**2 / sum_rho2
        return kl
    
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

        e_energy = np.zeros( self.N, dtype = np.float64)
        e_energy[:] = self.energies[:]
        if subtract_sscha:
            e_energy -= self.sscha_energies[:]

        # Compute the error using the Fortran Module
        value, error = SCHAModules.stochastic.average_error_weight(e_energy, self.rho, "err_yesrho")
        
        if return_error:
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
            - get_error : bool
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
            for i in xrange(nat):
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
        
        free_energy = self.current_dyn.GetHarmonicFreeEnergy(self.current_T)
        
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


    def get_free_energy_interpolating(self, target_supercell, support_dyn_coarse = None, support_dyn_fine = None, error_on_imaginary_frequency = True, return_error = False):
        """
        GET THE FREE ENERGY IN A BIGGER CELL
        ====================================

        This is a trick to interpolate the free energy in the
        infinite volume limit.

        Parameters
        ----------
            target_supercell : list (N, N, N)
               A list of three indices, where N is the dimension
               of the target supercell on which you want to interpolate.
            support_dyn[coarse/fine] : Phonons() Optional
               The harmonic dynamical matrix in the current/target_supercell
               This is optional, it can be used to achieve a better
               interpolation. If provided only the difference between
               the harmonic dyn and the current dyn is interpolated.
            error_on_imaginary_frequency : bool
               If Fase (default True) it will ignore imaginary frequencies
               arising from the interpolation. Otherwise an exception will
               be raised.
            return_error : bool
               As the normal get_free_energy, if this flag is True, the stochastic error is returned.

        Returns
        -------
            free_energy : float
               The free energy in the unit_cell volume [in Ry]. Note.
               This free energy is rescaled on the unit cell volume, 
               it is a different behaviour with respect to get_free_energy.
            error_on free energy : float
               The stochastic error, it is returned only if requested.
        """

        # Check if the support harmonic dyn is of the correct size.
        if not support_dyn_coarse is None:
            assert support_dyn_coarse.GetSupercell() == self.current_dyn.GetSupercell()
            assert support_dyn_fine.GetSupercell() == target_supercell

        
        # Interpolate the dynamical matrix
        new_dyn = self.current_dyn.Interpolate( self.current_dyn.GetSupercell(),
                                                target_supercell,
                                                support_dyn_coarse,
                                                support_dyn_fine)

        # TODO: Allow double interpolation in case of support dyn
        
        # Get the new harmonic free energy
        harm_fe = new_dyn.GetHarmonicFreeEnergy(self.current_T,
                                                not error_on_imaginary_frequency)
        harm_fe /= np.prod(target_supercell)

        # Get the average energy
        av_energy, av_error = self.get_average_energy(subtract_sscha = True, return_error = True)

        av_energy /= np.prod(self.current_dyn.GetSupercell())
        av_error /=  np.prod(self.current_dyn.GetSupercell())

        total_free_energy = harm_fe + av_energy

        if return_error:
            return total_free_energy, av_error
        return total_free_energy
                                                
            
    

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
                                    use_ups_supercell = True, preconditioned = 1,
                                    fast_grad = False, verbose = True):
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

        t1 = time.time()
        supercell_dyn = self.current_dyn.GenerateSupercellDyn(self.supercell)
        
        # Dyagonalize
        w, pols = supercell_dyn.DyagDinQ(0)
        trans = CC.Methods.get_translations(pols, supercell_dyn.structure.get_masses_array())
        ityp = supercell_dyn.structure.get_ityp() + 1 # Py to fortran convertion
        mass = np.array(list(supercell_dyn.structure.masses.values()))
        
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
        for i in xrange(self.N):
            u_disp[i, :, :] = np.reshape(self.u_disps[i,:], (nat, 3)) 
        
        
        # TODO: This may be dangerous
        pols = np.real(pols)
        
        t2 = time.time()
        if verbose:
            print( " [GRADIENT] Time to prepare the gradient calculation:", t2 -t1,"s")
            

        t1 = time.time()
        if fast_grad or not preconditioned:
            grad, grad_err = SCHAModules.get_gradient_supercell(self.rho, u_disp, eforces, w, pols, trans,
                                                                self.current_T, mass, ityp, log_err, self.N,
                                                                nat, 3*nat, len(mass), preconditioned)
        else:
            grad, grad_err = SCHAModules.get_gradient_supercell_new(self.rho, u_disp, eforces, w, pols, trans,
                                                                     self.current_T, mass, ityp, log_err, self.N,
                                                                     nat, 3*nat, len(mass))
        
            
        t2 = time.time()
        if verbose:
            print (" [GRADIENT] Time to call the fortran code:", t2 - t1, "s")
    
        # Perform the fourier transform
        q_grad = CC.Phonons.GetDynQFromFCSupercell(grad, np.array(self.current_dyn.q_tot),
                                                   self.current_dyn.structure, supercell_dyn.structure)
        q_grad_err = CC.Phonons.GetDynQFromFCSupercell(grad_err, np.array(self.current_dyn.q_tot),
                                                       self.current_dyn.structure, supercell_dyn.structure)
        t1 = time.time()
        if verbose:
            print (" [GRADIENT] Time to get back in fourier space:", t1 - t2, "s")
        
        
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
        r"""
        GET THE COVARIANCE STOCASTICALLY
        ================================
        
        This method is for testing, allows to use the ensemble to
        evaluate the covariance matrix stochastically. It should be equal
        to the matrix Upsilon^-1 that is obtained with the GetUpsilonMatrix method
        from the Phonons package.
        
        .. math::
            
            \Upsilon^{-1}_{ab} = \left< u_a u_b\right>
        
        Results
        -------
            cov_mat : 3nat x 3nat, ndarray
                A numpy matrix of the covariance matrix.
        """
        
        # A C style matrix of double precision real values
        cov_mat = np.einsum("i, ij, ik", self.rho, self.u_disps, self.u_disps) / np.sum(self.rho)
        
        return cov_mat
    
    
    def get_stress_tensor(self, offset_stress = None, add_centroid_contrib = False, use_spglib = False):

        """
        GET STRESS TENSOR
        =================
        
        The following subroutine computes the anharmonic stress tensor
        calling the fortran code get_stress_tensor.
        Note that the stress tensor is symmetrized to satisfy the cell constraint.
        
        NOTE: unit of measure is Ry/bohr^3 to match the quantum espresso one
        
        Parameters
        ----------
            offset_stress : 3x3 matrix, optional
                An offset stress to be subtracted to the real stress tensor.
                Usefull if you want to compute just the anharmonic contribution.
            add_centroid_contrib : bool, optional
                If true the contribution of the centroid is added. This is always zero when
                the system is relaxed.
            use_spglib : bool
                If true use the spglib library to perform the symmetrization

        
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
        
        
        # Get frequencies and polarization vectors
        super_dyn = self.current_dyn.GenerateSupercellDyn(self.supercell)
        wr, pols = super_dyn.DyagDinQ(0)
        trans = ~ CC.Methods.get_translations(pols, super_dyn.structure.get_masses_array())
        wr = np.real( wr[trans])
        pols = np.real( pols[:, trans])
        
        nat = super_dyn.structure.N_atoms 
        
        # Volume bohr^3
        volume = np.linalg.det(super_dyn.structure.unit_cell) * __A_TO_BOHR__**3
            
        
        # Get the correctly shaped polarization vectors
        er = np.zeros( (nat, len(wr), 3), dtype = np.float64, order = "F")
        
        for i in xrange(len(wr)):
            for j in xrange(nat):
                er[j, i, 0] = pols[3*j, i] 
                er[j, i, 1] = pols[3*j+1, i] 
                er[j, i, 2] = pols[3*j+2, i] 
                
        # Prepare the displacement in fortran order
        u_disps = np.zeros((self.N, nat, 3), dtype = np.float64, order = "F")
        for i in xrange(self.N):
            u_disps[i,:,:] = np.reshape(self.u_disps[i,:], (nat, 3))
        
        abinit_stress = np.einsum("abc -> cba", self.stresses, order = "F")
        
        stress, err_stress = SCHAModules.get_stress_tensor(volume, self.forces / __A_TO_BOHR__, u_disps * __A_TO_BOHR__, 
                                                           abinit_stress, wr, er, self.current_T, self.rho, "err_yesrho", 
                                                           self.N, nat, len(wr))    
        
        # Correct the stress adding the centroid contribution
        if add_centroid_contrib:
                
            eforces = self.forces - self.sscha_forces
            
            if not np.prod(self.supercell) == 1:
                # Refold the forces in the unit cell
                super_structure = self.current_dyn.structure.generate_supercell(self.supercell)
                itau = super_structure.get_itau(self.current_dyn.structure) - 1 # Fort -> Py
                nat = self.dyn_0.structure.N_atoms
                new_forces = np.zeros((self.N, nat, 3), dtype  =np.float64, order = "C")
                
                # Project in the unit cell the forces
                for i in xrange(nat):
                    #print "%d) ITAU LIST:" % i, itau == i
                    new_forces[:, i, :] = np.sum(eforces[:, itau==i,:], axis = 1) / np.prod(self.supercell)
                    #new_forces[:, i, :] = 
                
                eforces = new_forces
            
            stress_centr = np.zeros( (3,3), dtype = np.float64)
            error_centr = np.zeros( (3,3), dtype = np.float64)
            for i in xrange(0, 3):
                for j in xrange(i, 3):
                    av_array = 0.5 * np.einsum("h, ah", self.current_dyn.structure.coords[:, i],
                                               eforces[:,:,j])
                    av_array += 0.5 * np.einsum("h, ah", self.current_dyn.structure.coords[:, j],
                                                eforces[:,:,i])
                    stress_centr[i,j], error_centr[i,j] = SCHAModules.stochastic.average_error_weight(av_array, self.rho, "err_yesrho")
                    stress_centr[j,i] = stress_centr[i,j]
                    error_centr[j,i] = error_centr[i,j]
            
            
#            f, err_f = self.get_average_forces(True)
#            stress_centroid = 0.5 * np.einsum( "ai,aj", self.current_dyn.structure.coords * __A_TO_BOHR__, f) / volume
#            stress_centroid += np.transpose(stress_centroid)
#            err_stress_centroid = np.einsum( "ai,aj", self.current_dyn.structure.coords**2 , err_f**2)
#            err_stress_centroid = np.sqrt(err_stress_centroid) * __A_TO_BOHR__ / volume
#            err_stress_centroid = np.sqrt( err_stress_centroid**2 + np.transpose(err_stress_centroid**2))
#            divideby = np.ones( (3,3)) * 2
#            divideby[np.eye(3) == 1] = np.sqrt(2)
#            err_stress_centroid /= divideby

            stress += stress_centr
            err_stress = np.sqrt(err_stress**2 + error_centr**2)


        
        # Check the offset
        if not offset_stress is None:
            stress -= offset_stress

        # Symmetrize the stress tensor
        qe_sym = CC.symmetries.QE_Symmetry(self.current_dyn.structure)
        if not use_spglib:
            qe_sym.SetupQPoint()
        else:
            qe_sym.SetupFromSPGLIB()

        qe_sym.ApplySymmetryToMatrix(stress, err_stress)
        
        return stress, err_stress
    
    def get_average_stress(self):
        """
        GET THE AVERAGE STRESS
        ======================
        
        This gets only the ab-initio average of the stress tensor
        
        .. math::
                
            P_{\\alpha\\beta} = \\left<P_{\\apha\\beta}\\right>
        
        """
        stress = np.einsum("abc, a", self.stresses, self.rho) / np.sum(self.rho)
        qe_sym = CC.symmetries.QE_Symmetry(self.current_dyn.structure)
        qe_sym.SetupQPoint()
        qe_sym.ApplySymmetryToMatrix(stress)
        return stress
        
    
#     def get_free_energy_gradient_respect_to_dyn(self):
#         """
#         FREE ENERGY GRADIENT
#         ====================
        
#         Get the free energy gradient respect to the dynamical matrix.
#         The result is in [Ry/bohr^3] as the dynamical matrix are stored
#         in [Ry/bohr^2].
        
#         NOTE: Not working
        
#         .. math::
            
#             \\nabla_\\Phi \\mathcal F = -\\sum_{a\\mu} \\left<\\gamma_\\mu^a q_\\mu\\right>
            
#             \\gamma_\\mu^a = \\frac{e_\\mu^a \\nabla_\\Phi \\ln a_\\mu + \\nabla_\\Phi e_\mu^a}{\\sqrt M_a}(f_a - f^{\\Phi}_a)
            
#             q_\\mu = \\sum_b \\sqrt M_b e_\\mu^b (R_b - \\mathcal R_b)
            
#             \\nabla_\\Phi \\ln a_\\mu = \\frac{1}{2\\omega_\\mu a_\\mu} \\frac{\\partial a_\\mu}{\\partial\\omega_\\mu} \\frac{e_\\mu^a e_\\mu^b}{\\sqrt {M_aM_b}}
            
#             \\nabla_\\Phi e_\mu^c  =\\sum_{\\nu \\neq \\mu} \\frac{e_\\nu^a e_\\mu^b}{\\sqrt {M_aM_b} (\\omega_\\mu^2 - \\omega_\\nu^2)} e_\\nu^c
    
    
#         NOTE: it works only at gamma.
    
    
#         Return
#         ------
#             A 3Nx3N matrix. The gradient of the free energy (To be symmetrized)
            
#         """
#         #K_to_Ry=6.336857346553283e-06
        
#         #T = self.current_T
#         # TODO: TO BE TESTED
        
        
# #        # Get the mass vector
# #        _m_ = np.zeros(self.dyn_0.structure.N_atoms * 3)
# #        for i in range(self.current_dyn.structure.N_atoms):
# #            _m_[ 3*i : 3*i + 3] = self.current_dyn.structure.masses[ self.current_dyn.structure.atoms[i]]
# #        
# #        _m_sqrtinv = 1 / np.sqrt(_m_)
        
#         # Get the frequency and polarization vector of the dynamical matrix
#         w, pols = self.current_dyn.DyagDinQ(0)
        
        
#         # Discard translations and convert in Ha units
#         not_trans = ~CC.Methods.get_translations(pols, self.current_dyn.structure.get_masses_array())
#         w = np.array(w[not_trans] / 2, dtype = np.float64)
#         pols = np.real(pols[:, not_trans])
        
#         #n_modes = len(w)
        
#         # Convert the q vector into Ha units
#         q_new = np.array(self.current_q, dtype = np.float64) * np.sqrt(2) * __A_TO_BOHR__
        
#         # Get the ityp variable 
#         #ityp = self.current_dyn.structure.get_atomic_types()
        
#         # Get the mass and convert in Ha units
#         mass = np.array(self.current_dyn.structure.get_masses_array() * 2,
#                         dtype = np.float64)
        
#         nat = len(mass)
        
#         # Prepare the symmetrization
#         qe_sym = CC.symmetries.QE_Symmetry(self.current_dyn.structure)
#         qe_sym.SetupQPoint(self.current_dyn.q_tot[0])
        
        
# #        # Get the a_mu and its derivatives
# #        a_mu = np.zeros(n_modes, dtype = np.float64)
# #        da_dw = np.zeros(n_modes, dtype = np.float64)
        
#         # Use the fortran subroutines
# #        if T == 0:
# #            a_mu = 1 / np.sqrt(2* w) 
# #            da_dw = -1 /  np.sqrt(8 * w**3)
# #        else:            
# #            beta = 1 / (K_to_Ry*T)
# #            a_mu = 1 / np.sqrt( np.tanh(beta*w / 2) *2* w) 
# #            da_dw = - (w*beta + np.sinh(w*beta)) / (2 * np.sqrt(2) * w**2 * (np.cosh(beta*w) - 1) * np.sqrt(np.cosh(beta*w / 2) / (np.sinh(beta*w/2) * w)))
# #            
# #    

#         # Print the sscha forces converted
#         print ("SCHA forces:")
#         for i in xrange(self.N):
#             for j in xrange(self.current_dyn.structure.N_atoms):
#                 print ("Conf\t%d\tAtom\t%d\t" % (i, j), self.sscha_forces[i, j, :]/ (__A_TO_BOHR__))
                
                
#         # Convert the forces in Ha / bohr units and in the same type as fortran
#         e_forces = np.array( self.forces - self.sscha_forces, dtype = np.float64, order = "F") / (2 * __A_TO_BOHR__)
        
#         # Get df_da
#         df_da = SCHAModules.anharmonic.get_df_da_nonav(w, w, self.current_T, pols,
#                                                        e_forces,
#                                                        q_new, mass, "stat_schappp")
#         #print np.shape(e_forces)
#         # Now get the rest of the derivative
        
#         df_dfc = np.zeros( np.shape(self.current_dyn.dynmats[0]), dtype = np.float64)
#         err_df_dfc = np.zeros( np.shape(self.current_dyn.dynmats[0]), dtype = np.float64)
        
#         # Just to do something good
#         da_dcr_mat = np.zeros( (nat * 3, nat * 3, len(w)), dtype = np.float64)
        
#         for x_i in xrange(self.current_dyn.structure.N_atoms * 3):
#             for y_i in xrange(x_i, self.current_dyn.structure.N_atoms * 3):
#                 da_dcr, de_dcr = SCHAModules.anharmonic.get_da_dcr_and_de_dcr(w, pols, self.current_T,
#                                                                               mass, x_i+1, y_i+1)
                
#                 print ("(%d, %d): DA_DCR = " % (x_i+1, y_i+1), da_dcr)
#                 da_dcr_mat[x_i, y_i, :] = da_dcr
#                 da_dcr_mat[y_i, x_i, :] = da_dcr
                

#                 df_dc, delta_df_dc = SCHAModules.anharmonic.get_df_dcoeff_av_new(df_da, da_dcr, e_forces,
#                                                                                  q_new, mass, de_dcr, 
#                                                                                  self.rho, 1, "err_yesrho")
#                 # Fill the matrix
#                 df_dfc[x_i, y_i] = df_dc
#                 df_dfc[y_i, x_i] = df_dc
                
#                 err_df_dfc[x_i, y_i] = delta_df_dc
#                 err_df_dfc[y_i, x_i] = delta_df_dc
        
#         # Get the generator
#         ghr = np.zeros( (3*nat, 3*nat), dtype = np.float64, order = "F")
#         ghr[0,0] = 1
#         # Apply the sum rule
#         qe_sym.ImposeSumRule(ghr)
#         # Apply symmetries
#         qe_sym.SymmetrizeDynQ(ghr, self.current_dyn.q_tot[0])
#         ghr /= np.sqrt(np.trace(ghr.dot(ghr)))
#         print ("Generator:")
#         print (ghr)

#         print ("dA/dGhr = ", np.einsum("ijk, ij", da_dcr_mat, ghr) )       
        
#         # Force the symmetrization
#         qe_sym.ImposeSumRule(df_dfc)
#         qe_sym.ImposeSumRule(err_df_dfc)
#         qe_sym.SymmetrizeDynQ(df_dfc, self.current_dyn.q_tot[0])
#         qe_sym.SymmetrizeDynQ(err_df_dfc, self.current_dyn.q_tot[0])
        
#         # Convert from [Ha/bohr] in [Ry/bohr]
#         df_dfc *= 2
#         err_df_dfc *=  2
        

# #        # Prepare the w as a matrix
# #        _w_ = np.tile(w, (n_modes, 1))
# #        # 1 / (w_mu^2 - w_nu^2)
# #        one_over_omegamunu = 1 / (_w_**2 - _w_.transpose()**2)
# #        #one_over_omegamunu *= 1 - np.eye(n_modes) # Remove the therms for mu equal to nu
# #        one_over_omegamunu[ (_w_ - _w_.transpose()) < __EPSILON__] = 0
# #        
# #        #print "freqs:", w
# #        #print "w", _w_
# #        #print "one_over_omega:", one_over_omegamunu
# #                                        
# #        # Get the derivative of the lna_mu respect to the dynamical matrix
# #        # Inner product
# #        d_lna_d_dyn = np.einsum("i, ai, bi, ci, a, b, c->abic", da_dw/(2 * w * a_mu), pols, pols, pols, _m_sqrtinv, _m_sqrtinv, _m_sqrtinv)
# #        
# #        # Get the derivative respect to the polarization vector
# #        d_pol_d_dyn = np.einsum("ai,bj,ci,ji,a,b,c->abjc", pols, pols, pols, one_over_omegamunu, _m_sqrtinv, _m_sqrtinv, _m_sqrtinv)
# #        
# #        #print "d_lna:", d_lna_d_dyn
# #        #print "d_pol:", d_pol_d_dyn
# #        
# #        pre_sum = d_lna_d_dyn + d_pol_d_dyn
# #        
# #        # Get the q vector
# #        d_F_d_dyn = np.zeros(np.shape(self.current_dyn.dynmats[0]))
# #        for i in range(self.N):
# #            # Get the displacements of the structure
# #            u_disp = self.structures[i].get_displacement(self.current_dyn.structure).reshape(3 * self.current_dyn.structure.N_atoms)
# #            
# #            # Get the forces on the configuration
# #            delta_f = (self.forces[i,:,:] - self.sscha_forces[i,:,:]).reshape(3 * self.current_dyn.structure.N_atoms)
# #            
# #            # Get the q vector
# #            q = np.einsum("i, ij, i", np.sqrt(_m_), pols, u_disp)
# #            
# #            # Get gamma matrix
# #            gamma = np.einsum("abcd, d", pre_sum, delta_f)
# #            
# #            #print "%d) delta_f = " % (i+1), delta_f
# #            #print "%d) q = " % (i+1), q
# #            #print "%d) gamma = " % (i+1), gamma
# #            
# #            # Contract the gamma matrix and multiply it for the weight
# #            partial_gradient = - np.einsum("abc, c", gamma, q)
# #            d_F_d_dyn += partial_gradient * self.rho[i]
# #            
# #            #print "conf %d | weight %.4e | partial gradient:" % (i, self.rho[i]), partial_gradient
# #            
# #            
# #        # Normalization
# #        d_F_d_dyn /= np.sum(self.rho)
#         #print "Grad:"
#         #for i in range(np.shape(d_F_d_dyn)[0]):
#         #    print " ".join(["%7.2e" % x for x in list(d_F_d_dyn[i,:])])
        
#         #TODO: apply symmetries
            
#         return df_dfc, err_df_dfc


    def get_d3_muspace(self):
        r"""
        GET V3 IN MODE SPACE
        ====================

        This subroutine gets the d3 directly in the space of the modes.

        ..math::

            D^{(3)}_{abc} = \sum_{xyz} \frac{\Phi^{(3)}_{xyz} e_a^x e_b^y e_c^z}{\sqrt{m_x m_y m_z}}


        """

        # Be shure to have the correct units
        self.convert_units(UNITS_DEFAULT)

        supersturct = self.current_dyn.structure.generate_supercell(self.supercell)

        # Convert from A to Bohr the space 
        u_disps = self.u_disps * __A_TO_BOHR__
        n_rand, n_modes = np.shape(u_disps)
        forces = (self.forces - self.sscha_forces).reshape(self.N, n_modes)  / __A_TO_BOHR__ 

        Ups = self.current_dyn.GetUpsilonMatrix(self.current_T)
        v_disp = u_disps.dot(Ups)

        # pass in the polarization space
        w, pols = self.current_dyn.DiagonalizeSupercell()

        # Discard translations
        trans = CC.Methods.get_translations(pols, supersturct.get_masses_array())
        pols = pols[:, ~trans]

        m = np.tile(supersturct.get_masses_array(), (3,1)).T.ravel()

        pol_vec = np.einsum("ab, a->ab", pols, 1 / np.sqrt(m))

        v_mode = v_disp.dot(pol_vec)
        f_mode = forces.dot(pol_vec)

        # Now compute the d3 as <vvf>
        N_eff = np.sum(self.rho)
        f_mode = np.einsum("ia, i->ia", f_mode, self.rho)
        d3_noperm = np.einsum("ia,ib,ic->abc", v_mode, v_mode, f_mode)
        d3_noperm /= -N_eff # there is a minus

        # Apply the permuatations
        d3 = d3_noperm.copy()
        d3 += np.einsum("abc->acb", d3_noperm)
        d3 += np.einsum("abc->bac", d3_noperm)
        d3 += np.einsum("abc->bca", d3_noperm)
        d3 += np.einsum("abc->cab", d3_noperm)
        d3 += np.einsum("abc->cba", d3_noperm)
        d3 /= 6

        # TODO: symmetrize

        return d3


    def get_v3_realspace(self):
        """
        This is a testing function that computes the V3 matrix in real space:
            
        ..math::
            
            \\Phi^{(3)}_{xyz} = - \sum_{pq} \\Upsilon_{xp}\\Upsilon_{yq} \\left<u_pu_q f_z\\\right>
        """
        
        nat_sc = self.structures[0].N_atoms
        Ups = self.current_dyn.GetUpsilonMatrix(self.current_T)
        f2 = np.reshape(self.forces - self.sscha_forces, (self.N, 3 * nat_sc))
        
        # Get the average <uuf>
        t1 = time.time()
        N_eff = np.sum(self.rho)
        uuf = np.einsum("ix,iy,iz,i", self.u_disps, self.u_disps, f2, self.rho)
        uuf /= N_eff
        t2 = time.time()
        
        print("Time elapsed to compute <uuf>:", t2-t1, "s")
        
        # Get the v3
        v3 = - np.einsum("xp,yq,pqz->xyz", Ups, Ups, uuf) * __A_TO_BOHR__
        # Symmetrize
        #v3 = np.einsum("xyz,xzy,yxz,yzx,zxy,zyx->xyz", v3, v3, v3, v3, v3, v3) / 6
        v3 -=  np.einsum("xp,yq,pqz->xzy", Ups, Ups, uuf) * __A_TO_BOHR__
        v3 -=  np.einsum("xp,yq,pqz->zyx", Ups, Ups, uuf) * __A_TO_BOHR__
        v3 -=  np.einsum("xp,yq,pqz->yxz", Ups, Ups, uuf)* __A_TO_BOHR__
        v3 -=  np.einsum("xp,yq,pqz->zxy", Ups, Ups, uuf)* __A_TO_BOHR__
        v3 -=  np.einsum("xp,yq,pqz->yzx", Ups, Ups, uuf)* __A_TO_BOHR__
        v3 /= 6
        t3 = time.time()
        print("Time elapsed to compute v3:", t3-t1, "s")
        return v3

    def get_odd_realspace(self):
        """
        This is a testing function to compute the odd3 correction 
        using the real space v3 (similar to the raffaello first implementation)
        """
        # Get the dynamical matrix in the supercell
        super_dyn = self.current_dyn.GenerateSupercellDyn(self.supercell)
        w_sc, pols_sc = super_dyn.DyagDinQ(0)
        
        # Remove translations
        no_trans_mask = ~CC.Methods.get_translations(pols_sc, super_dyn.structure.get_masses_array())
        w_sc = w_sc[no_trans_mask]
        pols_sc = pols_sc[:, no_trans_mask]
        
        # Get phi3
        phi3 = self.get_v3_realspace()
        
        # Get Gmunu
        Gmunu = super_dyn.GetGmunu(self.current_T)
        
        # Divide the polarization vectors by the mass
        nat_sc = super_dyn.structure.N_atoms
        epol = pols_sc.copy()
        m = super_dyn.structure.get_masses_array()
        for x in range(nat_sc):
            epol[3*x : 3*x + 3, :] = pols_sc[3*x : 3*x + 3, :] / np.sqrt(m[x])
        
        first_part = np.einsum("xab,ij,ai,bj->xij", phi3, Gmunu, epol, epol)
        second_part = np.einsum("ci,dj,cdy->ijy", epol, epol, phi3)
        odd_correction = np.einsum("xij, ijy->xy", first_part, second_part)
        
        fakedyn = super_dyn.Copy()
        fakedyn.dynmats[0] = odd_correction
        fakedyn.save_qe("odd_new")
        return super_dyn.dynmats[0] + odd_correction
        

    def get_v3_qspace(self, q, k):
        r"""
        GET THE PHONON-PHONON SCATTERING ELEMENT
        ========================================

        This subroutine computes the 3-body phonon-phonon scatternig within the sscha.
        It evaluates the vertex where the q phonon splits in a q+k and -k phonon:
                    
                  /---> q + k
           q ____/
                 \
                  \---> -k  

        
        This computes v3 on the fly in real space.
        
        .. math ::

            V^3_{abc} (q, -q-k, k)

        Where :math:`a`, :math:`b` and :math:`c` are the atomic indices in the unit cell.

        Parameters
        ----------
            q : ndarray(size = 3, dtype = float)
                The q vector for the v3 compuation V3(q, -q-k, k).
            k : ndarray(size = 3, dtype = float)
                The k vector for the v3 computation V3(q, -q-k, k).
        
        Returns
        -------
            v3 : ndarray( size = (3*nat, 3*nat, 3*nat), dtype = np.complex128)
                The 3-rank tensor vertext V3(q, -q-k, k) of the phonon-phonon scattering
        """

        # Define the q vectors
        q1 = -q -k
        q2 = k

        superdyn = self.current_dyn.GenerateSupercellDyn(self.supercell)
        superstruc = superdyn.structure
        ups_mat = np.real(superdyn.GetUpsilonMatrix(self.current_T))

        # Get Upsilon dot u
        vs = self.u_disps.dot(ups_mat) * __A_TO_BOHR__
        
        # Get the corrispondance between unit cell and super cell
        itau = superdyn.structure.get_itau(self.current_dyn.structure) - 1
        nat_sc = superdyn.structure.N_atoms
        nat = self.current_dyn.structure.N_atoms
        struct = self.current_dyn.structure

        D3 = np.zeros( (3*nat, 3*nat, 3*nat), dtype = np.complex128)
        N_eff = np.sum(self.rho)
        fc = np.zeros((3,3,3), dtype = np.float64)
        for i in range(nat_sc):
            i_uc = itau[i]

            # The forces and displacement along this atom
            v_i = vs[:, 3*i:3*i+3]
            f_i = self.forces[:, i, :] - self.sscha_forces[:, i, :]
            f_i /= __A_TO_BOHR__
            for j in range(nat_sc):
                j_uc = itau[j]

                R1 = superstruc.coords[i, :] - struct.coords[i_uc, :]
                R1 -= superstruc.coords[j, :] - struct.coords[j_uc, :]
                q1dotR = q1.dot(R1)

                # Forces and displacement along this atom
                v_j = vs[:, 3*j:3*j+3]
                f_j = self.forces[:, j, :] - self.sscha_forces[:, j, :] 
                f_j /= __A_TO_BOHR__
                

                for k in range(nat_sc):
                    k_uc = itau[k]

                    R2 = superstruc.coords[i, :] - struct.coords[i_uc, :]
                    R2 -= superstruc.coords[k, :] - struct.coords[k_uc, :]
                    q2dotR = q2.dot(R2)
                        
                    # Forces and displacement along this atom
                    v_k = vs[:, 3*k:3*k+3]
                    f_k = self.forces[:, k, :] - self.sscha_forces[:, k, :] 
                    f_k /= __A_TO_BOHR__

                    fc[:,:] = np.einsum("ia,ib,ic,i", v_i, v_j, f_k, self.rho)
                    fc += np.einsum("ia,ib,ic,i", v_i, f_j, v_k, self.rho)
                    fc += np.einsum("ia,ib,ic,i", f_i, v_j, v_k, self.rho)
                    fc /= 3*N_eff
                    D3[3*i_uc: 3*i_uc+3, 3*j_uc: 3*j_uc+3, 3*k_uc : 3*k_uc+3] += fc * np.exp(-1j* q1dotR - 1j*q2dotR)

                    
        return D3


    def get_dynamical_bubble(self, q, w, smearing = 1e-5):
        r"""
        GET THE DYNAMICAL BUBBLE SELF ENERGY
        ====================================

        This function returns the dynamical bubble self-energy:

        .. math::

            \Sigma_{af}(q, w) = \sum_{q'q''}\sum_{bc,\mu\nu} D^{(3)}_{abc} \left(-\frac 1 2 \chi_{\mu\nu}(\omega, q', q'')\right) \frac{e_\nu^b e_\mu^c e_\nu^d e_\mu^e}{\sqrt{M_bM_cM_dM_e}} D^{(3)}_{def}

        
        NOTE: The integral in the q space is performed over the mesh grid given by the supercell.


        Parameters
        ----------
            q : vector
                The q vector to compute the dynamical self energy
            w : float or array
                The frequency(ies) to compute the dynamical self-energy

        Results
        -------
            Sigma : ndarray(size = (3*nat, 3*nat), dtype = np.complex128)
                The dynamical self energy. Note it could be a list of Sigma
                if the provided frequency is an array
        """


        # Perform the summation over the allowed q points
        q_list = self.current_dyn.q_tot

        bg = CC.Methods.get_reciprocal_vectors(self.current_dyn.structure.unit_cell)
        m = self.current_dyn.structure.get_masses_array()

        nat = self.current_dyn.structure.N_atoms

        # Extend m to 3*nat
        # This is a numpy hack: tile creates a replica matrix of m (3xN_nat)
        # .T: makes a transposition to N_nat x 3 and ravel convert it in a 1d array
        m = np.tile(m, (3,1)).T.ravel()
        minvsqrt = 1 / np.sqrt(m)

        # Check how many frequencies has been provided
        N_w = 1
        try:
            N_w = len(w)
        except:
            pass

        # Initialize the bubble self energy
        sigma = np.zeros((3*nat, 3*nat), dtype = np.complex128)
        if N_w > 1:
            sigmas = [sigma.copy() for x in range(N_w)]
        for ik, k in enumerate(q_list):
            k1 = -q -k

            # Get the v3 
            d3_1 = self.get_v3_qspace(k1, k)
            print ("Sum of v3:", np.sum(d3_1))

            # Get the phonon-propagator
            if N_w > 1:
                bubbles = []
                for i in range(N_w):
                    bubbles.append( self.current_dyn.get_phonon_propagator(w[i], self.current_T, k1, k, smearing))
            else:
                bubble = self.current_dyn.get_phonon_propagator(w, self.current_T, k1, k, smearing)

            # Get the index of k1 to extract the polarization vectors
            k1_dists = [CC.Methods.get_min_dist_into_cell(bg, k1, x) for x in q_list]
            k1_i = np.argmin(k1_dists)

            wk, polk = self.current_dyn.DyagDinQ(ik)
            wk1, polk1 = self.current_dyn.DyagDinQ(k1_i)

            # Get |e><e| / sqrt(m_a m_b)
            e_mat = np.einsum("am,bn, a, b->abmn", polk1, polk, minvsqrt, minvsqrt)

            # Convert the d3 in the mu basis
            d3_mubasis = np.einsum("abc, bcmn -> amn", d3_1, e_mat)

            # Compute the bubble
            if N_w > 1:
                for i in range(N_w):
                    sigmas[i] -= np.einsum("amn, mn, bmn->ab", d3_mubasis, bubbles[i], np.conj(d3_mubasis)) / 2
            else:
                sigma -= np.einsum("amn, mn, bmn->ab", d3_mubasis, bubble, np.conj(d3_mubasis)) / 2

        if N_w > 1:
            return sigmas
        return sigma


    def get_free_energy_hessian(self, include_v4 = False, get_full_hessian = True, verbose = False, \
        use_symmetries = True):
        """
        GET THE FREE ENERGY ODD CORRECTION
        ==================================

        This subroutines computes the odd correction
        to the free energy hessian using the fortran subroutines, as describe in the
        Bianco paper ...

        The calculation is performed in the supercell

        Parameters
        ----------
            include_v4 : bool
                If True we include the fourth order force constant matrix.
                This requires a lot of memory
            get_full_hessian : bool
                If True the full hessian matrix is returned, if false, only the correction to
                the SSCHA dynamical matrix is returned.
            verbose : bool
                If True a lot things are written in output.
                This is usefull for debugging purpouses.
            use_symmetries : bool
                If true, the d3 and d4 are symmetrized in real space.
                It requires that spglib is installed to detect symmetries in the supercell correctly.

        Returns
        -------
            phi_sc : Phonons()
                The dynamical matrix of the free energy hessian in (Ry/bohr^2)
        """
        # For now the v4 is not implemented
        #     if include_v4:
        #         ERROR_MSG = """
        # Error, the v4 computation has not yet been implemented.
        # """
        #         raise NotImplementedError(ERROR_MSG)

        # Convert anything into the Ha units
        # This is needed for the Fortran subroutines
        self.convert_units(UNITS_HARTREE)

        # Get the dynamical matrix in the supercell
        dyn_supercell = self.current_dyn.GenerateSupercellDyn(self.supercell)
        w, pols = dyn_supercell.DyagDinQ(0)
        a = SCHAModules.thermodynamic.w_to_a(w, self.current_T)
        

        n_modes = len(w)
        nat_sc = int(np.shape(pols)[0] / 3)

        # Get the polarization vectors in the correct format
        new_pol = np.zeros( (nat_sc, n_modes, 3), dtype = np.double)
        for i in range(nat_sc):
            for j in range(n_modes):
                new_pol[i, j, :] = pols[3*i : 3*(i+1), j]
        

        # Get the translational modes
        trans = CC.Methods.get_translations(pols, dyn_supercell.structure.get_masses_array())


        # Get the atomic types
        ityp = dyn_supercell.structure.get_ityp() + 1 #Py to Fortran indexing
        n_typ = len(self.current_dyn.structure.masses)

        amass = np.zeros(n_typ, dtype = np.double)

        for at_type in self.current_dyn.structure.masses:
            index = ityp[self.current_dyn.structure.atoms.index(at_type)] - 1
            amass[index] = self.current_dyn.structure.masses[at_type]

        # Get the forces and conver in the correct units
        f = (self.forces - self.sscha_forces)# * Bohr 
        u = self.u_disps.reshape((self.N, nat_sc, 3), order = "C") #/ Bohr

        log_err = "err_yesrho"

        # Lets call the Fortran subroutine to compute the v3
        if verbose:
            print ("Going into d3")
        d3 = SCHAModules.get_v3(a, new_pol, trans, amass, ityp,
                                f, u, self.rho, log_err)
        if verbose:
            print("Outside d3")
        

        # Symmetrize the d3
        if use_symmetries:
            if verbose:
                print("Symmetrizing the d3")
                np.save("d3_realspace_nosym.npy", d3)
            qe_sym = CC.symmetries.QE_Symmetry(dyn_supercell.structure)
            qe_sym.SetupFromSPGLIB()
            qe_sym.ApplySymmetryToTensor3(d3)

        if verbose:
            print("Saving the third order force constants as d3_realspace_sym.npy [Ha units]")
            np.save("d3_realspace_sym.npy", d3)

        # Check if the v4 must be included
        if include_v4:
            print("Computing the v4, this requires some time...")
            t1 = time.time()
            d4 = SCHAModules.get_v4(a, new_pol, trans, amass, ityp, \
                f, u, self.rho, log_err)
            t2 = time.time()
            print("Time elapsed to compute the v4: {} s".format(t2-t1))

            # Symmetrize the v4 
            if use_symmetries:
                qe_sym = CC.symmetries.QE_Symmetry(dyn_supercell.structure)
                qe_sym.SetupFromSPGLIB()
                qe_sym.ApplySymmetryToTensor4(d4)
            
            if verbose: print("Inside odd straight")
            phi_sc_odd = SCHAModules.get_odd_straight_with_v4(a, w, new_pol, trans, \
                amass, ityp, self.current_T, d3, d4)
            if verbose : print("Outside odd straight")
        else:
            # Only v3
            # Get the odd correction (In Ha/bohr^2)
            if verbose: 
                print ("Inside odd straight")
                print (" A = ", a)
                print (" W = ", w)
                print (" TRANS = ", trans)
                print (" AMASS = ", amass)
                print (" ITYP = ", ityp)
                print (" T = ", self.current_T)
            phi_sc_odd = SCHAModules.get_odd_straight(a, w, new_pol, trans, amass, ityp, 
                                                    self.current_T, d3)
                                        
            if verbose:
                print ("Outside odd straight.")
                print ("Saving the odd correction (Ha) as phi_odd.npy")
                np.save("phi_odd.npy", phi_sc_odd)

                # Try to save this matrix
                dyn_supercell.dynmats[0] = phi_sc_odd 
                dyn_supercell.save_qe("SupercellOddDynHa")
        


        # Lets fourier transform
        dynq_odd = CC.Phonons.GetDynQFromFCSupercell(phi_sc_odd, np.array(self.current_dyn.q_tot), 
                                                     self.current_dyn.structure, dyn_supercell.structure)
        
        
        # Convert back the ensemble in Default units
        self.convert_units(UNITS_DEFAULT)
        dynq_odd *= 2 # Ha/bohr^2 -> Ry/bohr^2

        # Generate the Phonon structure by including the odd correction
        dyn_hessian = self.current_dyn.Copy()
        for iq in range(len(self.current_dyn.q_tot)):
            if get_full_hessian:
                dyn_hessian.dynmats[iq] = self.current_dyn.dynmats[iq] + dynq_odd[iq, :, :] 
            else:
                dyn_hessian.dynmats[iq] = dynq_odd[iq, :, :] 

        
        return dyn_hessian






    def get_odd_correction(self, include_v4 = False, store_v3 = True, 
            store_v4 = True, progress = False, frequencies = 0, smearing = 5e-5, v4_conv_thr = 1e-2,
            return_only_correction = False, save_all = False, load_d3 = None, use_omp = False):
        """
        RAFFAELLO BIANCO'S ODD CORRECTION TESTING VERSION
        =================================================

        NOTE: THIS IS DEPRECATED! USE get_free_energy_hessian instead!
        
        This function returns the odd correction to the dynamical matrix. 
        It can be used to get the free energy curvature, to study structure stability.
        It contains also the static phonons (w -> 0). 
        
        .. math ::
            
            \\Phi^(F) = \\Phi^{(2)} + \\Phi^{(3)} \\Lambda (I - \\Phi^{(4)})^{-1} \\Lambda \\Phi^{(3)}
        
        where :math:`\\Phi^{(2)}` is the standard dynamical matrix, while :math:`\\Phi^{(n)}` are
        the high order force constant matrices.
        
        See "Phys. Rev. B 96, 014111" for further details (consider also to cite it if
        you use this method).
        
        Parameters
        ----------
            include_v4 : bool
                If false (default) the :math:`\\Phi^{(4)} matrix is not computed and it is
                considered as the null matrix. This speed up a lot the calculation,
                and :math:`\\Phi^{(4)}` is found to give a neglectable contribution in many
                systems.
            store_v3 : bool
                If true (default) the d3 matrix is stored in memory. This is much
                cheaper from a computational point of view, but in large system can
                require a lot of memory.
            store_v4 : bool
                If true (default) the v4 matrix is stored in memory. This is very fast, but
                it can drain the memory.
            progress : bool
                If true (default false), shows a progress bar while computing 
            frequencies : double or list
                if different from zero (default = 0) it computes the dynamical correction.
                If it is a list it will compute the correction for each frequency in the list.
                In this case the returned quantity will be a list of self-energies.
            smearing : double
                This is the smearing applied only in the case of dynamical correction
            return_only_correction: bool, optional
                If true it returns only the correction behond the sscha matrix.
            save_all : bool, default false
                If true the self-energy in the supercell is saved for each w value. 
                It is usefull if you want to analyze the results in a new run.
            load_d3 : ndarray((n_modes, n_modes, n_modes), dtype = np.float64)
                The d3 matrix. If None (default) it is recomputed.
            use_omp : bool
                If true the OpenMP parallelization is used to perform some calculation.
                If you enable this flag be carefull to the scaling with the number of processors.
                It can have a negative impact. You can select the number of threads
                by using the OMP_NUM_THREADS variable before running the python command.
        
        Results
        -------
            new_dyn : ndarray(3xnat_sc, 3xnat_sc, dtype = float64)
                The real dynamical matrix. Note, if return_only_correction = True, it
                will be only the bubble (+ bubble chain) terms. 
        """
        # Check if freqiencies is an array
        # And setup if it is a static or dynamic calcuation
        is_dynamic = False
        try:
            N_w = len(frequencies)
            is_dynamic = True
        except:
            N_w = 1
            if np.abs(frequencies) < __EPSILON__:
                is_dynamic = False
        

        # Setup the type of the calculation
        __TYPE__ = np.float64
        if is_dynamic:
            __TYPE__ = np.complex128
        
        # Get the dynamical matrix in the supercell
        super_dyn = self.current_dyn.GenerateSupercellDyn(self.supercell)
        w_sc, pols_sc = super_dyn.DyagDinQ(0)
        
        # Remove translations
        no_trans_mask = ~CC.Methods.get_translations(pols_sc, super_dyn.structure.get_masses_array())
        w_sc = w_sc[no_trans_mask]
        pols_sc = pols_sc[:, no_trans_mask]
        
        # Check if the polarization vectors have an imaginary part
        imag_part = np.sqrt(np.sum( np.imag(pols_sc)**2))
        if imag_part > __EPSILON__:
            print("Imaginary part of the polarization vectors:", imag_part)
            raise ValueError("Error the imaginary part of the polarization vector exceeded the threshold.")
            
        pols_sc = np.real(pols_sc)
        
        nat_sc = super_dyn.structure.N_atoms
        
        # Get the masses array of the correct shape
        _m_ = np.zeros( (3*nat_sc), dtype = __TYPE__)
        _m_[:] = np.tile(super_dyn.structure.get_masses_array(), (3, 1)).T.ravel()
        
        n_modes_sc = len(w_sc)
        
        # Get the forces of the correct shape
        _f_ = (self.forces - self.sscha_forces).reshape((self.N, 3*nat_sc))
        
        # Get the X and Y arrays
        X = np.zeros( (n_modes_sc, self.N), dtype = __TYPE__)
        Y = np.zeros( (n_modes_sc, self.N), dtype = __TYPE__)

        # Get the Upsilon eigenvalue
        w_ups = 2 * w_sc        
        n_bos = 0 * w_sc
        dn_dw = 0 * w_sc
        if self.current_T > __EPSILON__:
            n_bos = 1 / (np.exp(w_sc * __RyToK__ / self.current_T) - 1)
            w_ups /=  1 + 2*n_bos
            dn_dw = - n_bos**2 * np.exp(w_sc * __RyToK__ / self.current_T) * __RyToK__ / self.current_T
            
        X[:,:] = np.einsum("ab,ca,a,b -> bc", pols_sc, self.u_disps, np.sqrt(_m_), w_ups) * __A_TO_BOHR__
        Y[:,:] = np.einsum("ab,ca,a,c -> bc", pols_sc, _f_, 1 / np.sqrt(_m_), self.rho) / __A_TO_BOHR__
        N_eff = np.sum(self.rho)
        
        if save_all:
            # Save the X and Y to file.
            # This is for debugging.
            np.save("X.npy", X)
            np.save("Y.npy", Y)

        # Get lambda matrix
        if N_w > 1:
            Lambdas = []

        # Get the propagator at all the given frequencies
        for i_lambda in range(N_w):
            if N_w > 1:
                frequency = frequencies[i_lambda]
            else:
                frequency = frequencies
        
            Lambda_G = np.zeros( (n_modes_sc, n_modes_sc), dtype = __TYPE__)
            for mu in range(n_modes_sc):
                w_mu = w_sc[mu]
                n_mu = n_bos[mu]
                
                for nu in range(n_modes_sc):
                    w_nu = w_sc[nu]
                    n_nu = n_bos[nu]
                    if abs(w_mu - w_nu) < __EPSILON__ and not is_dynamic:
                        Lambda_G[mu, nu] = (2*n_nu + 1)/(2*w_nu) - dn_dw[nu]
                    else:
                        if not is_dynamic:
                            Lambda_G[mu, nu] = (n_mu + n_nu + 1)/(w_mu + w_nu) - (n_mu - n_nu)/(w_mu - w_nu)
                        else:
                            Lambda_G[mu, nu] =(w_mu + w_nu)*(n_mu + n_nu + 1)/((w_mu + w_nu)**2 - (frequency + 1j*smearing)**2) - \
                                (w_mu - w_nu)*(n_mu - n_nu)/((w_mu - w_nu)**2 - (frequency + 1j*smearing)**2)

                    Lambda_G[mu, nu] /= - 4*w_mu*w_nu
            if N_w > 1:
                Lambdas.append(Lambda_G)
                
        #print ("Consistent check on Lambda:", np.sqrt(np.sum( (Lambda_G - Lambda_G.T)**2)))
        
        # Compute the odd3 correction

        if N_w > 1:
            odd_corrs = []
        #odd_corr = np.zeros( (n_modes_sc, n_modes_sc), dtype = __TYPE__)
        if store_v3:
            if progress:
                print("Computing v3...")
            d3 = np.zeros( (n_modes_sc, n_modes_sc, n_modes_sc), dtype = np.float64, order = "C")
            if load_d3 is None:
                if not use_omp:
                    # Compute the d3 using python
                    d3 = np.einsum("ai,bi,ci", X, X, Y)
                    d3 += np.einsum("ai,bi,ci", X, Y, X)
                    d3 += np.einsum("ai,bi,ci", Y, X, X)
                    d3 /= - 3 * N_eff
                else:
                    # Lets call the C code with openMP support
                    # to compute the d3 faster

                    if not _SSCHA_ODD_:
                        raise ImportError("Error, sscha_HP_odd is required to use openmp, thid is deprecated. Pleas use get_free_energy_hessian instead.")

                    sscha_HP_odd.GetV3(X, Y, n_modes_sc, self.N, d3)
                    d3 *= self.N / N_eff
            else:
                d3 = load_d3
            #d3 /= N_eff

            print (type(d3[0,0,0]))
            print ("N_eff:", N_eff)

            #odd_corr[:,:] = np.einsum("abc,bc,de,def -> af", d3, Lambda_G, Lambda_G, d3)
            if not include_v4:
                if progress:
                    print("Computing the bubble...")
                for i_freq in range(N_w):
                    if N_w > 1:
                        Lambda_G = Lambdas[i_freq]
                    odd_corr = np.einsum("abc,bc,bcd -> ad", d3, Lambda_G, d3)
                    if N_w > 1: 
                        odd_corrs.append(odd_corr)

            if include_v4:

                if store_v4:
                    d4 = np.einsum("ai,bi,ci,di", X, X, X, Y)
                    d4 += np.einsum("ai,bi,ci,di", X, X, Y, X)
                    d4 += np.einsum("ai,bi,ci,di", X, Y, X, X)
                    d4 += np.einsum("ai,bi,ci,di", Y, X, X, X)
                    d4 /= -4 * N_eff

                    print ("d3 = ", np.sqrt(np.sum(d3**2)))
                    print ("d3 * Lambda = ", np.sqrt(np.sum(np.einsum("ab,abc", Lambda_G, d3)**2)))
                    print ("d4 = ", np.sqrt(np.sum(d4**2)))
                    print ("Lambda = ", np.sqrt(np.sum(Lambda_G**2)))
                    print ("d4*Lambda =", np.sqrt(np.sum(np.einsum("ab, abxy", Lambda_G, d4)**2)))
                    # Reshape Lambda d4

                    # Get the dynamical response
                    for i_freq in range(N_w):
                        if N_w > 1:
                            Lambda_G = Lambdas[i_freq]
                        Ld4 = np.einsum("ab, abxy->abxy", Lambda_G, d4).reshape((n_modes_sc**2, n_modes_sc**2))
                        print ("Ld4 spectrum:", np.linalg.eigvals(Ld4))

                        # Get the solution
                        Inv_mat = np.linalg.inv(np.eye(n_modes_sc**2) - Ld4)
                        Inv_mat_4d = Inv_mat.reshape((n_modes_sc, n_modes_sc, n_modes_sc, n_modes_sc))
                        New_Prop = np.einsum("abcd, cd->abcd", Inv_mat_4d, Lambda_G)
                        odd_corr = np.einsum("xab, abcd, cdy", d3, New_Prop, d3)
                        if N_w > 1:
                            odd_corrs.append(odd_corr)
                    """ 
                    # Perform the cycle for the geometric sum
                    old_odd = odd_corr.copy()
                    new_corr = np.zeros( np.shape(old_odd), dtype = __TYPE__)
                    running = True
                    new_d3 = d3.copy()

                    # Setup the biconjugate method
                    R = np.einsum("xyab, ab, abz", d4, Lambda_G, d3)
                    Rbar = R.copy()
                    P = R.copy()
                    Pbar = R.copy()
                    
                    while running:
                        Ap = P - np.einsum("xyab, ab, abz", d4, Lambda_G, P)
                        ATpbar = Pbar - np.einsum("abxy,xy, abz->xyz", d4, Lambda_G, Pbar)

                        alpha = np.einsum("abc, abc->c", Rbar, R) / np.einsum("abc,abc->c", Pbar, Ap)
                        new_d3 += np.einsum("c, abc->abc", alpha, P)
                        Rnew = R - np.einsum("c, abc->abc", alpha, Ap)
                        Rbarnew = Rbar - np.einsum("c, abc->abc", alpha, ATpbar)
                        beta = np.einsum("abc,abc->c", Rbarnew, Rnew) / np.einsum("abc, abc->c", Rbar, R)
                        P = R + np.einsum("c, abc->abc", beta, P)
                        Pbar = Rbar + np.einsum("c, abc->abc", beta, Pbar)
                        R = Rnew
                        Rbar = Rbarnew

                        rest_mod =np.sqrt(np.einsum("abc,abc", R, R))
                        if progress:
                            print("")
                            print ("Rest: ", rest_mod)

                        if rest_mod < v4_conv_thr:
                            running = False   """
                else:
                    raise NotImplementedError("Not yet implemented")
                    # We can do exactly the same as before
                    # But performing the sum over I explicitely
                    odd_dim = np.sum(odd_corr**2)
                    old_odd = odd_corr.copy()
                    new_corr = np.zeros(np.shape(old_odd), dtype = __TYPE__)
                    running = True
                    new_d3 = d3.copy()


                    while running:
                        for x in xrange( n_modes_sc):
                            aux1 = new_d3[x, :, :] * Lambda_G
                            new_d3[x, :, :] = np.einsum("ab,ai,bi,ci,di->cd", aux1, X, X, X, Y)
                            new_d3[x, :, :] += np.einsum("ab,ai,bi,ci,di->cd", aux1, X, X, Y, X)
                            new_d3[x, :, :] += np.einsum("ab,ai,bi,ci,di->cd", aux1, X, Y, X, X)
                            new_d3[x, :, :] += np.einsum("ab,ai,bi,ci,di->cd", aux1, Y, X, X, X)
                            new_d3[x, :, :] /= -4*N_eff
                            aux2 = new_d3[x, :, :] * Lambda_G
                            for y in xrange(x, n_modes_sc):
                                new_corr[x,y] = np.sum(aux2 * d3[y, :, :])
                                new_corr[y,x] = new_corr[x,y]
                            
                            if progress:
                                sys.stdout.write("\r Computing v4 correction ... %d / %d" % (x+1, n_modes_sc))
                        
                        # Add the correction
                        odd_corr += new_corr
                        ratio = np.sum(new_corr**2) / np.sum(odd_corr**2)
                        
                        if progress:
                            print("")
                            print("Ratio: %.4e | Converged after: %.4e" % (ratio, v4_conv_thr))
                            print("")
                        if np.sqrt(ratio) < v4_conv_thr:
                            running = False  


                    # Do not store d4
                    #raise NotImplementedError("Error, v4 not stored not implemented")
            
#            for b in range(len(w_sc)):
#                for c in range(len(w_sc)):
#                    odd_corr[:,:] += np.outer(d3[:, b, c], d3[b, c, :]) * Lambda_G[b,c]
        else:
            if progress:
                print ("D3 correction computation:")
                
            
            for a in range(n_modes_sc):
                d3_a = np.einsum("i,bi,ci", X[a,:], X, Y)
                d3_a += np.einsum("i,bi,ci", X[a,:], Y, X)
                d3_a += np.einsum("i,bi,ci", Y[a,:], X, X)
                d3_a /= - 3 * N_eff
                for b in range(a, n_modes_sc):
                    d3_b = d3_a
                    if a != b:
                        d3_b = np.einsum("i,bi,ci", X[a,:], X, Y)
                        d3_b += np.einsum("i,bi,ci", X[a,:], Y, X)
                        d3_b += np.einsum("i,bi,ci", Y[a,:], X, X)
                        d3_b /= - 3 * N_eff
                    
                    odd_corr[a,b] = np.sum(d3_a * d3_b * Lambda_G)
                    odd_corr[b,a] = odd_corr[a,b]
                
                if progress:
                    sys.stdout.write("\r%2d %%" % (a * 100 / n_modes_sc))
                    sys.stdout.flush()

                if include_v4:
                    raise NotImplementedError("Error, v4 without v3 stored not yet implemented.")

                    
            
#                
#            
#            odd_corr[:,:] += np.einsum("ai,bi,ci,bc,bi,ci,di->ad", X, X, Y, Lambda_G, X, X, Y)
#            print ("1")
#            odd_corr[:,:] += np.einsum("ai,bi,ci,bc,bi,ci,di->ad", X, X, Y, Lambda_G, X, Y, X)
#            print ("2")
#            odd_corr[:,:] += np.einsum("ai,bi,ci,bc,bi,ci,di->ad", X, X, Y, Lambda_G, Y, X, X)
#            print ("3...")
#            odd_corr[:,:] += np.einsum("ai,bi,ci,bc,bi,ci,di->ad", X, Y, X, Lambda_G, X, X, Y)
#            odd_corr[:,:] += np.einsum("ai,bi,ci,bc,bi,ci,di->ad", X, Y, X, Lambda_G, X, Y, X)
#            odd_corr[:,:] += np.einsum("ai,bi,ci,bc,bi,ci,di->ad", X, Y, X, Lambda_G, Y, X, X)
#            odd_corr[:,:] += np.einsum("ai,bi,ci,bc,bi,ci,di->ad", Y, X, X, Lambda_G, X, X, Y)
#            odd_corr[:,:] += np.einsum("ai,bi,ci,bc,bi,ci,di->ad", Y, X, X, Lambda_G, X, Y, X)
#            odd_corr[:,:] += np.einsum("ai,bi,ci,bc,bi,ci,di->ad", Y, X, X, Lambda_G, Y, X, X)
#                
#                
#            for i in range(self.N):
#                odd_corr[:,:] += np.einsum("a,b,c,bc,b,c,d->ad", X[:,i], X[:,i], Y[:,i], Lambda_G, X[:,i], X[:,i], Y[:,i])
#                odd_corr[:,:] += np.einsum("a,b,c,bc,b,c,d->ad", X[:,i], X[:,i], Y[:,i], Lambda_G, X[:,i], Y[:,i], X[:,i])
#                odd_corr[:,:] += np.einsum("a,b,c,bc,b,c,d->ad", X[:,i], X[:,i], Y[:,i], Lambda_G, Y[:,i], X[:,i], X[:,i])
#                odd_corr[:,:] += np.einsum("a,b,c,bc,b,c,d->ad", X[:,i], Y[:,i], X[:,i], Lambda_G, X[:,i], X[:,i], Y[:,i])
#                odd_corr[:,:] += np.einsum("a,b,c,bc,b,c,d->ad", X[:,i], Y[:,i], X[:,i], Lambda_G, X[:,i], Y[:,i], X[:,i])
#                odd_corr[:,:] += np.einsum("a,b,c,bc,b,c,d->ad", X[:,i], Y[:,i], X[:,i], Lambda_G, Y[:,i], X[:,i], X[:,i])
#                odd_corr[:,:] += np.einsum("a,b,c,bc,b,c,d->ad", Y[:,i], X[:,i], X[:,i], Lambda_G, X[:,i], X[:,i], Y[:,i])
#                odd_corr[:,:] += np.einsum("a,b,c,bc,b,c,d->ad", Y[:,i], X[:,i], X[:,i], Lambda_G, X[:,i], Y[:,i], X[:,i])
#                odd_corr[:,:] += np.einsum("a,b,c,bc,b,c,d->ad", Y[:,i], X[:,i], X[:,i], Lambda_G, Y[:,i], X[:,i], X[:,i])
#      
#                if progress:
#                    if i % 10 == 0:
#                        sys.stdout.write("\r%2d %%" % (i * 100 / self.N))
#                        sys.stdout.flush()
#            
#            odd_corr /= (3 * self.N)**2
            if progress:
                print ()
            #raise NotImplementedError("Error, the store_v3 = .false. has still not been implemented")
        
        # Get back in the cartesian basis
        total_self_energies = []
        total_self_energy = np.zeros((3*nat_sc, 3*nat_sc), dtype = __TYPE__)
        for i_freq in range(N_w):
            if N_w > 1:
                odd_corr = odd_corrs[i_freq]
            cart_odd_corr = np.einsum("ai,bj,ij -> ab", pols_sc, pols_sc, odd_corr)
            cart_odd_corr *= np.outer(np.sqrt(_m_),  np.sqrt(_m_))


            total_self_energy = cart_odd_corr 
            if not return_only_correction:
                if is_dynamic:
                    total_self_energy += super_dyn.dynmats[0]
                else:
                    total_self_energy += np.real(super_dyn.dynmats[0])
            if N_w >1:
                total_self_energies.append(total_self_energy)
            
            if save_all:
                np.save("odd_corr_%d.npy" % i_freq, cart_odd_corr)
        
        # Return the new dynamical matrix in the supercell
        if N_w > 1:
            return total_self_energies
        return total_self_energy


    def compute_ensemble(self, calculator, compute_stress = True, stress_numerical = False,
                         cluster = None):
        """
        GET ENERGY AND FORCES
        =====================
        
        This is the generic function to compute forces and stresses.
        It can be used both with clusters, and with simple ase calculators
        
        Paramters
        ---------
            calculator:
                The ase calculator
            compute_stress: bool
                If true compute the stress
            stress_numerical : bool
                Compute the stress tensor with finite difference, 
                this is not possible with clusters
            cluster: Cluster, optional
                The cluster in which to send the calculation.
                If None the calculation is performed on the same computer of
                the sscha code.
        """

        
        # Check if the calculator is a cluster
        is_cluster = False
        if not cluster is None:
            is_cluster = True
        
        # Check consistency
        if stress_numerical and is_cluster:
            raise ValueError("Error, stress_numerical is not implemented with clusters")
    
        # Check if not all the calculation needs to be done
        n_calcs = np.sum( self.force_computed.astype(int))
        computing_ensemble = self

        if compute_stress:
            self.has_stress = True

        # Check wheter compute the whole ensemble, or just a small part
        should_i_merge = False 
        if n_calcs != self.N:
            should_i_merge = True
            computing_ensemble = self.get_noncomputed()
            self.remove_noncomputed() 

        if is_cluster:
            cluster.compute_ensemble(computing_ensemble, calculator, compute_stress)
        else:
            computing_ensemble.get_energy_forces(calculator, compute_stress, stress_numerical)
        
        if should_i_merge:
            # Remove the noncomputed ensemble from here, and merge 
            self.merge(computing_ensemble)


    def merge(self, other):
        """
        MERGE TWO ENSEMBLES
        ===================

        This function will merge two ensembles together. 
        
        Parameters
        ----------
            other : Ensemble()
                Another ensemble to be merge with. It must be generated by the same dynamical matrix
                as this one, otherwise wired things will happen.
        """

        self.N += other.N 
        self.forces = np.concatenate( (self.forces, other.forces), axis = 0)
        self.stresses = np.concatenate( (self.stresses, other.stresses), axis = 0)
        self.structures += other.structures
        self.u_disps = np.concatenate((self.u_disps, other.u_disps), axis = 0)
        self.xats = np.concatenate((self.xats, other.xats), axis = 0)
        self.energies = np.concatenate( (self.energies, other.energies))

        self.stress_computed = np.concatenate( (self.stress_computed, other.stress_computed))
        self.force_computed = np.concatenate( (self.force_computed, other.force_computed))


        self.sscha_forces = np.concatenate( (self.sscha_forces, other.sscha_forces), axis = 0)
        self.sscha_energies = np.concatenate( (self.sscha_energies, other.sscha_energies))

        self.rho = np.concatenate( (self.rho, other.rho))

        # Now update everything
        self.update_weights(self.current_dyn, self.current_T)


    def split(self, split_mask):
        """
        SPLIT THE ENSEMBLE
        ==================

        This method will return an ensemble with only the configurations matched by the split_mask array.
        NOTE: The original ensemble will remain untouched.

        Parameters
        ----------
            split_mask : ndarray(size = self.N, dtype = bool)
                A mask array. It must be of the same size of the number of configurations, 
                and contain a True or False if you want that the corresponding configuration to be included in the
                splitted ensemble
        
        Results
        -------
            splitted_ensemble : Ensemble()
                An ensemble tath will contain only the configurations in the split mask.
        """

        structs = [self.structures[x] for x in np.arange(len(split_mask))[split_mask]]

        N = np.sum(split_mask.astype(int))
        ens = Ensemble(self.dyn_0, self.T0, self.dyn_0.GetSupercell())
        ens.init_from_structures(structs) 
        ens.force_computed[:] = self.force_computed[split_mask]
        ens.stress_computed[:] = self.stress_computed[split_mask]
        ens.energies[:] = self.energies[split_mask]
        ens.forces[:, :, :] = self.forces[split_mask, :, :]
        ens.has_stress = self.has_stress
        if self.has_stress:
            ens.stresses[:, :, :] = self.stresses[split_mask, :, :]

        ens.update_weights(self.current_dyn, self.current_T)

        return ens


    def remove_noncomputed(self):
        """
        Removed all the incomplete calculation from the ensemble.
        It may be used to run a minimization even if the ensemble was not completely calculated.
        """

        good_mask = self.force_computed
        if self.has_stress:
            good_mask = good_mask & self.stress_computed

        self.N = np.sum( good_mask.astype(int))
        self.forces = self.forces[good_mask, :, :]    
        self.sscha_forces = self.sscha_forces[good_mask, :, :]    
        self.stresses = self.stresses[good_mask, :, :]      
        self.energies = self.energies[good_mask]
        self.sscha_energies = self.sscha_energies[good_mask]
        self.xats = self.xats[good_mask, :, :]
        self.u_disps = self.u_disps[good_mask, :]

        self.structures = [self.structures[x] for x in np.arange(len(good_mask))[good_mask]]

        self.rho = self.rho[good_mask]

        # Check everything and update the weights
        self.update_weights(self.current_dyn, self.current_T)

    def get_noncomputed(self):
        """
        Get another ensemble with only the non computed configurations.
        This may be used to resubmit only the non computed values
        """

        non_mask = ~self.force_computed
        if self.has_stress:
            non_mask = non_mask & (~self.stress_computed)

        return self.split(non_mask)

    def get_energy_forces(self, ase_calculator, compute_stress = True, stress_numerical = False, skip_computed = False):
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
            stress_numerical : bool
                If the calculator does not support stress, it can be computed numerically
                by doing finite differences.
            skip_computed : bool
                If true the configurations already computed will be skipped. 
                Usefull if the calculation crashed for some reason.
            
        """
        
        # Setup the calculator for each structure
        if __MPI__:
            comm = MPI.COMM_WORLD
            size = comm.Get_size()
            rank = comm.Get_rank()
            
            # Broad cast to all the structures
            structures = comm.bcast(self.structures, root = 0)            
            nat3 = comm.bcast(self.current_dyn.structure.N_atoms* 3* np.prod(self.supercell), root = 0)
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
            nat3 = self.current_dyn.structure.N_atoms* 3 * np.prod(self.supercell)
            N_rand = self.N
            
        # Only for the master
        
        # Prepare the energy, forces and stress array
        energies = np.zeros( N_rand // size, dtype = np.float64)
        forces = np.zeros( ( N_rand // size) * nat3 , dtype = np.float64)
        if compute_stress:
            stress = np.zeros( (N_rand // size) * 9, dtype = np.float64)

        if rank == 0:
            total_forces = np.zeros( N_rand * nat3, dtype = np.float64)
            total_stress = np.zeros( N_rand * 9, dtype = np.float64)
        else:
            total_forces = np.empty( N_rand * nat3, dtype = np.float64)
            total_stress = np.empty( N_rand * 9, dtype = np.float64)
            

        # If an MPI istance is running, split the calculation
        for i0 in xrange(N_rand // size):
            i = i0 + size * rank

            # Avoid performing this calculation if already done
            if skip_computed:
                if self.force_computed[i]:
                    if compute_stress:
                        if self.stress_computed[i]:
                            continue
                    else:
                        continue
            
            
            struct = structures[i]
            atms = struct.get_ase_atoms()
            
            # Setup the ASE calculator
            atms.set_calculator(ase_calculator)


            # Print the status
            if rank == 0:
                print ("Computing configuration %d / %d" % (i0+1, N_rand / size))
            
            # Avoid for errors
            run = True
            count_fails = 0
            while run:
                try:
                    energy = atms.get_total_energy() / Rydberg # eV => Ry
                    run = False
                except:
                    print ("Rerun the job %d" % i)
                    count_fails += 1
                    if count_fails >= 5:
                        run = False
                        sys.stderr.write("Error in the ASE calculator for more than 5 times\n")
                        raise
            
            # Get energy, forces (and stress)
            energy = atms.get_total_energy() / Rydberg # eV => Ry
            forces_ = atms.get_forces() / Rydberg # eV / A => Ry / A
            if compute_stress:
                if not stress_numerical:
                    stress[9*i0 : 9*i0 + 9] = -atms.get_stress(False).reshape(9) * Bohr**3 / Rydberg  # ev/A^3 => Ry/bohr
                else:
                    stress[9*i0 : 9*i0 + 9] = -ase_calculator.calculate_numerical_stress(atms, voigt = False).ravel()* Bohr**3 / Rydberg 
            
            # Copy into the ensemble array
            energies[i0] = energy
            forces[nat3*i0 : nat3*i0 + nat3] = forces_.reshape( nat3 )
            
            

            
        
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
        self.forces[:, :, :] = np.reshape(total_forces, (N_rand, self.current_dyn.structure.N_atoms*np.prod(self.supercell), 3), order = "C")
        self.force_computed[:] = True

        if compute_stress:
            self.stresses[:,:,:] = np.reshape(total_stress, (N_rand, 3, 3), order = "C")
            self.has_stress = True
            self.stress_computed[:] = True
        else:
            self.has_stress = False
            
            
        
        
