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
import sscha.Ensemble
import sscha.SchaMinimizer

import numpy as np

# Rydberg to cm-1 conversion factor
RyToCm  = 109691.40235

# Load the harmonic dynamical matrix
harm_dyn = CC.Phonons.Phonons("dyn_harmonic", full_name = True)
cell = harm_dyn.structure.unit_cell * 1.01
new_dyn = harm_dyn.GetStrainMatrix(cell, T = 0)
harm_dyn.dynmats[0] = new_dyn.dynmats[0]

# Load the initial sscha matrix
start_dyn = CC.Phonons.Phonons("start_dyn", full_name = True)

w_harm, p = harm_dyn.DyagDinQ(0)
w_sscha, p = start_dyn.DyagDinQ(0)

print "Starting freqs | exact freqs"
print "\n".join(["\t".join( (str(w_sscha[i]), str(w_harm[i]))) for i in range(len(w_harm))])

# Temperature of the simulation
T = 0 #K

# The number of configurations
N = 10

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
minim.min_step_dyn = 1e-5
minim.min_step_struc = 1e-5

# Perform M minimization steps
for i in range(M):
    
    # Print a set of info about the step
    print " ---------------- "
    print " Step :", i
    print " Free energy :", minim.get_free_energy()
    print " Weights :", minim.ensemble.rho
    w, pols = minim.dyn.DyagDinQ(0)
    print " Freqs :", "\t".join(["%.3f cm-1" % (freq * RyToCm) for freq in w])
    if i != 0:
        print " |gc| :", np.sqrt(np.sum(np.diag(np.dot( minim.prev_grad, minim.prev_grad))))
    print ""
    minim.minimization_step()