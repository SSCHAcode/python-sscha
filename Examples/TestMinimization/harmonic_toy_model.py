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
import cellconstructor.Methods

import sscha
import sscha.Ensemble
import sscha.SchaMinimizer

import spglib

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['figure.dpi'] = 200
#plt.ion()

# Rydberg to cm-1 conversion factor
RyToCm  = 109691.40235

# Load the harmonic dynamical matrix
harm_dyn = CC.Phonons.Phonons("dyn_harmonic", full_name = True)
#cell = harm_dyn.structure.unit_cell * 1.01
#new_dyn = harm_dyn.GetStrainMatrix(cell, T = 0)
#harm_dyn.dynmats[0] = new_dyn.dynmats[0]

# Load the initial sscha matrix
#start_dyn = CC.Phonons.Phonons("start_dyn", full_name = True)
start_dyn = harm_dyn.Copy()

w_harm, p = harm_dyn.DyagDinQ(0)
w_sscha, p = start_dyn.DyagDinQ(0)

print "Starting freqs | exact freqs"
print "\n".join(["\t".join( (str(w_sscha[i]), str(w_harm[i]))) for i in range(len(w_harm))])

# Temperature of the simulation
T = 0 #K

# The number of configurations
N = 1600

# The number of minimization steps
M = 5

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

# Test gradient
MOVE_STEP = 3e-3

print "----- MOVING DYN ------- "
print "# Step, free energy IS, error, free energy STOC, error, gradient IS, gradient STOC, Effective sample size"
steps = np.zeros(M)
grad_IS = np.zeros(M)
free_IS = np.zeros(M)
free_IS_err = np.zeros(M)
grad_STOC = np.zeros(M)
free_STOC = np.zeros(M)
free_STOC_err = np.zeros(M)
kong_liu = np.zeros(M)
av_energy_STOC= np.zeros((M,2))
av_energy_IS = np.zeros((M,2))
sscha_energy = np.zeros(M)

frequencies = np.zeros( (M, 6))

for i in range(M):
    minim.update()
    
    
    free_IS[i], free_IS_err[i] = minim.ensemble.get_free_energy(return_error = True)
    g = minim.ensemble.get_free_energy_gradient_respect_to_dyn()
    grad_IS[i] = .5 * (g[5,0] + g[0,5])
    av_energy_IS[i,0], av_energy_IS[i, 1]  = minim.ensemble.get_average_energy(return_error = True, subtract_sscha = False)
    
    # Generate a new ensemble
    test_ensemble = sscha.Ensemble.Ensemble(minim.dyn, T)
    test_ensemble.generate(N)
    print np.sum(minim.ensemble.rho), np.sum(test_ensemble.rho)
    
        
    for k in range(N):
        energy, force = harm_dyn.get_energy_forces(test_ensemble.structures[k])
        test_ensemble.energies[k] = energy
        test_ensemble.forces[k, :,:] = force
                            
    test_ensemble.update_weights(minim.dyn, T)
    print test_ensemble.energies - test_ensemble.sscha_energies
                       
    free_STOC[i], free_STOC_err[i] = test_ensemble.get_free_energy(return_error = True)
    av_energy_STOC[i,0], av_energy_STOC[i, 1]  = test_ensemble.get_average_energy(return_error = True, subtract_sscha= False)

    g = test_ensemble.get_free_energy_gradient_respect_to_dyn()
    grad_STOC[i] = .5* (g[5,0] + g[0,5])
    
    steps[i] = i * MOVE_STEP
    kong_liu[i] = minim.ensemble.get_effective_sample_size() / float(N)
    
    
    # Save the ensemble
    test_ensemble.save("data_dir", population = i)
    test_ensemble.dyn_0.save_qe("dyn_pop%d_" % i)
    
    
    print "HARM:"
    w, p = harm_dyn.DyagDinQ(0)
    print w
    print "Current:"
    w,p = test_ensemble.current_dyn.DyagDinQ(0)
    print w
    frequencies[i, :] = w.copy()
    
    # Get the sscha contribution to the energy
    sscha_energy[i] = np.sum(w[~ CC.Methods.get_translations(p)] / 2)
         
    print "%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e" % (steps[i], free_IS[i], free_IS_err[i], free_STOC[i], free_STOC_err[i], grad_IS[i], grad_STOC[i], kong_liu[i])
    
    minim.dyn.dynmats[0][5,0] += MOVE_STEP / np.sqrt(2)
    minim.dyn.dynmats[0][0,5] += MOVE_STEP / np.sqrt(2)
    
    # Apply the sum rule 
    minim.dyn.ApplySumRule()
    
    # Get the symmetries from the matrix
    print "Symmetry CLASS:", spglib.get_spacegroup(minim.dyn.structure.get_ase_atoms())
    sym_spglib = spglib.get_symmetry(minim.dyn.structure.get_ase_atoms())
    symmetries = CC.Methods.GetSymmetriesFromSPGLIB(sym_spglib)
    print "Symmetry counts:", len(symmetries)
    
    # Symmetrize the matrix
    minim.dyn.ForceSymmetries(symmetries)
    
# Plot the Free energy
plt.figure()
plt.errorbar(steps, free_IS, yerr = free_IS_err, label="IS")
plt.errorbar(steps, free_STOC, yerr = free_STOC_err, label = "STOC")
plt.xlabel("Step")
plt.ylabel("Free energy")
plt.title("Free energy test")
plt.legend()
plt.tight_layout()

plt.figure()
plt.errorbar(steps, av_energy_IS[:,0], yerr = av_energy_IS[:,1], label = "IS")
plt.errorbar(steps, av_energy_STOC[:,0], yerr = av_energy_STOC[:,1], label = "STOC")
plt.plot(steps, sscha_energy[:], label="SSCHA energy")
plt.legend()
plt.xlabel("Step")
plt.ylabel("Energy")
plt.title("Total energy comparison")
plt.tight_layout()


# Plot the gradient
plt.figure()
plt.plot(steps, grad_IS, label="IS")
plt.plot(steps, grad_STOC, label="STOC")
plt.plot(steps[:-1] + MOVE_STEP * .5, np.diff(free_IS) / MOVE_STEP, label="IS finite difference")
plt.plot(steps[:-1] + MOVE_STEP * .5, np.diff(free_STOC) / MOVE_STEP, label="STOC finite difference")
plt.xlabel("Step")
plt.ylabel("gradient along the direction")
plt.title("Gradient test")
plt.legend()
plt.tight_layout()

plt.figure()
plt.plot(steps, kong_liu)
plt.title("Kong-Liu effective sample size")
plt.xlabel("Step")
plt.ylabel("$N_{eff} / N$")
plt.tight_layout()
plt.show()
    
    

## Perform M minimization steps
#for i in range(M):
#    
#    # Print a set of info about the step
#    print " ---------------- "
#    print " Step :", i
#    print " Free energy :", minim.get_free_energy()
#    print " Weights :", minim.ensemble.rho
#    w, pols = minim.dyn.DyagDinQ(0)
#    print " Freqs :", "\t".join(["%.3f cm-1" % (freq * RyToCm) for freq in w])
#    if i != 0:
#        print " |gc| :", np.sqrt(np.sum(np.diag(np.dot( minim.prev_grad, minim.prev_grad))))
#    print ""
#    minim.minimization_step()