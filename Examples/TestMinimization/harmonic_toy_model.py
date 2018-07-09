# -*- coding: utf-8 -*-

"""
This example shows a prototype of a minimization
on an harmonic system.

The harmonic dynamical matrix is loaded from file, an
harmonic toy model is generated on this matrix, 
then the matrix is slightly modified and the minimization is
started with the new matrix.

In the end we check if the final matrix coincides with the harmonic one.
"""

import cellconstructor as CC
import cellconstructor.Structure
import cellconstructor.Phonons

import sscha

# Load the harmonic dynamical matrix
harm_dyn = CC.Phonons.Phonons("dyn_harmonic", full_name = True)

# Load the initial sscha matrix
start_dyn = CC.Phonons.Phonons("start_dyn", full_name = True)

# Temperature of the simulation
T = 100 #K

# The number of configurations
N = 1000

# The number of minimization steps
M = 100

# Generate the ensemble
ensemble = sscha.Ensemble.Ensemble(start_dyn, T)
ensemble.generate(N)


# Get the forces and energies with the harmonic dynamical matrices
for i in range(N):
    energy, force = harm_dyn.get_energy_forces(ensemble.structures[i])
    ensemble.energies[i] = energy
    ensemble.forces[i, :,:] = force


# Setup the minimization
minim = sscha.SchaMinimizer.SSCHA_Minimizer(ensemble)

# Perform M minimization steps
for i in range(M):
    minim.minimization_step()
    
    # Print a set of info about the step
    print " ---------------- "
    print " Step :", i
    print " Free energy :", minim.get_free_energy()
    print ""