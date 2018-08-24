# -*- coding: utf-8 -*-

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


"""
This module contains usefull subroutines to work with
"""

import cellconstructor as CC
import scipy, scipy.linalg
import numpy as np

# ------------------ GENERATORS ---------------------------
# Here we work with the old sscha generators, to enable them
# to represent any given matrix

class Generators:
    """
    This is a class that allows for the interface with the old
    Fortran sscha code, that is interely based on the generators.
    """
    
    def __init__(self):
        """
        Initialize the generator class
        """
        
        # The number of q points represented by the generators
        self.nq = 0
        self.nat = 0
        self.wyck_gen = None
        self.wyck_ncoeff = 0
        
        self.dyn_ncoeff = 0
        self.dyn_gen = []
    
    def LoadFromFileFC(self, filename, natoms, nqirr):
        """
        This subroutine loads the list of generators from a file.
        The FC means that the generators are expected to be on the force constant
        matrix, not the wyckoff positions.
        
        Parameters
        ----------
            filename : string
                The path to the file that contains all the generators
            natoms : int
                The number of atoms in the current structure
            nqirr : int
                The number of irreducible q points
        """
        
        f = open(filename, "r")
        flines = [l.strip() for l in f.readlines()]
        
        self.dyn_gen = []
        self.nat = natoms
        
        # Read how many generator are for this particular q point
        n_gen = int(flines[0])
        current_i = 0
        
        ghrs = []
        fc = np.zeros( (3*natoms, 3*natoms), dtype = np.complex128) 
        new_gen = True
        
        na = 0
        nb = 0
        index  = 0
        for i, line in enumerate(flines):
            
            if new_gen:
                # Append the generator
                if i != 0:
                    ghrs.append(fc)
                
                if current_i == n_gen:
                    n_gen = int(line)
                    current_i = 0
                    self.dyn_gen.append( np.array(ghrs))
                    self.dyn_ncoeff.append(len(ghrs))
                    continue
            
                
                fc = np.zeros( (3*natoms, 3*natoms), dtype = np.complex128)
                new_gen = False
                current_i += 1
                continue
            
            # Polish the line
            line = line.replace(",", "")
            line = line.replace("(", "")
            line = line.replace(")", "")
            
            line_list = [np.float64(x) for x in line.split()]
            
            # Select the atomic indexs
            if len(line_list) == 2:
                na = line_list[0] -1
                nb = line_list[1] -1
                index = 0
                continue
            
            fc[3 * na + index, 3 * nb] = line_list[0] + 1j*line_list[1]
            fc[3 * na + index, 3 * nb + 1] = line_list[2] + 1j*line_list[3]
            fc[3 * na + index, 3 * nb + 2] = line_list[4] + 1j*line_list[5]
            
            if (na+1  == natoms ) and (nb+1 == natoms):
                new_gen = True
            
                
        # Append also the last generators
        self.dyn_gen.append( np.array(ghrs))
        self.dyn_ncoeff.append(len(ghrs))
        
        
    def ProjectDyn(self, fc, iq = 0):
        """
        Project the force constant matrix in the
        basis of the generators
            
        
        Parameters
        ----------
            fc : ndarray (3n x 3n)
                The force constant matrix to be projected on the generator subspace
            iq : the index of the irreducible q point.
                
        Results
        -------
            ndarray :
                The coefficients of the generators that decompose the number.
        """
        
        res = np.einsum("ijk, kj", self.dyn_gen[iq], fc)
        return res                
    
    
    def GetDynFromCoeff(self, coeffs, iq=0):
        """
        This subroutine generate the dynamical matrix starting from the generator coefficients
        
        
        Parameters
        ----------
            coeffs : ndarray
                The coefficients that represent the dynamical matrix. 
                Must be of the correct dimension.
            iq : int
                The index of the q point
                
                
        Result
        ------
            ndarray 3N x 3N
                The fc generated by the coefficients.
        """
        
        # Check if the coeff are of the correct length
        if len(coeffs) != self.dyn_ncoeff[iq]:
            raise ValueError("Error, the number of coeff %d does not match the number of generator %d. (iq=%d)" % (len(coeffs), self.dyn_coeff[iq], iq))
        
        fc = np.einsum("ijk, i", self.dyn_gen[iq], coeffs)
        return fc
    
    def Generate(self, dyn, qe_sym = None):
        """
        GENERATE THE GENERATORS
        =======================
        
        
        The following subroutine generate the generators for the given dynamical
        matrix and the given symmetries.
        
        NOTE: this subroutine must be test for supercells, in particular complex generator should be added
        NOTE: This must be tested in general
        Parameters
        ----------
            dyn : CC.Phonons.Phonons()
                The dynamical matrix represented by the generators.
            qe_sym : CC.symmetries.QE_Symmetry()
                If given, the selected symmetries will be used to generate
                the generators. Otherwise symmetries will be generated from the 
                dynamical matrix using the default parameters.
        """
        
        # Check if the symmetries must be initialize
        if qe_sym is None:
            qe_sym = CC.symmetries.QE_Symmetry(dyn.structure)
            
        
        # Get the number of irreducible q points from the matrix
        self.nq = dyn.nqirr
        self.nat = dyn.structure.N_atoms
        
        # Initialize the symmetries at q = 0
        qe_sym.SetupQPoint()
        
        # Prepare the wyckoff basis
        tmp_wyck_gen = np.zeros((3 * self.nat, self.nat, 3), dtype = np.float64)
        
        for i in range( 3 * self.nat):
            x = i % 3
            n = i / 3
            tmp_wyck_gen[i, n, x] = 1
            
            # Symmetrize the vector
            qe_sym.SymmetrizeVector(tmp_wyck_gen[i, :, :])
        
        # Apply the gram-schmidt
        new_gen = tmp_wyck_gen.reshape((3 * self.nat, 3 * self.nat)).transpose()
        new_gen = scipy.linalg.orth(new_gen).transpose()
        
        # Get the number of wyckoff coefficients
        self.wyck_ncoeff = new_gen.shape()[0]
        
        # Reshape the array and get the coefficients
        self.wyck_gen = new_gen.reshape((self.wyck_ncoeff, self.nat, 3))
        
        r = np.arange(3 * self.nat)
        
        self.dyn_ncoeff = np.zeros(self.nq, dtype = int)
        self.dyn_gen = []
        
        # Cycle for each irreducible q point of the matrix
        for iq in range(self.nq):
            q = dyn.q_stars[iq][0]
            # Setup the symmetries for this q point
            qe_sym.SetupQPoint(q)
            
            gh = []
            
            for i in range(self.nat * 3):
                for j in range(i, self.nat * 3):
                    # Take the generator
                    fc = np.zeros((3 * self.nat, 3 * self.nat), dtype = np.complex128)
                    fc[i, j] = 1
                    
                    # Apply the symmetry
                    qe_sym.SymmetrizeDynQ(q, fc)
                    
                    # Check if the generator has already be generated
                    is_new = True
                    for k in range(i+1):
                        mask = fc[k, :] != 0
                        first_value = r[mask]
                        if len(first_value):
                            if k == i:
                                if first_value[0] < j:
                                    is_new = False
                                    break
                        else:
                            is_new = False
                            break
                    
                    # If the generator is new
                    if is_new:
                        qe_sym.ImposeSumRule(fc, "simple")
                        
                        # Check if the sum rule makes this generator desappearing
                        if np.sum ((fc != 0).as_type(int)) != 0:
                            gh.append(fc / np.sqrt(np.trace(fc.dot(fc))))
        
            dim = len(gh)
        
            # Prepare the gram-shmidt
            gh = np.array(gh, dtype = np.complex128)
        
            gh_new = np.reshape((dim, 9 * self.nat**2)).transpose()
            gh_new = scipy.linalg.orth(gh_new).transpose()
        
            self.dyn_ncoeff = np.shape(gh_new)[0]
        
            self.dyn_gen.append(np.reshape(gh_new, (self.dyn_ncoeff, 3*self.nat, 3*self.nat)))
            
                            
                    
                    
                    
                
            