 # -*- coding: utf-8 -*-

"""
TEST GRADIENT

This example generate an input harmonic matrix that uses as a toy model
to compute forces and energies.
The dynamical matrix is moved in one of the allowed direction by symmetry.
The ensemble is generated and the Free energy and its gradient is computed.
The gradient is computed also by generating a new ensemble to test the
importance sampling.

This example shows how to manipulate the sscha methods in detail.
"""


import cellconstructor as CC
import cellconstructor.Structure
import cellconstructor.Phonons
import cellconstructor.Methods
import cellconstructor.symmetries

import sscha
import sscha.Ensemble
import sscha.SchaMinimizer


import numpy as np
import matplotlib.pyplot as plt





plt.rcParams['figure.dpi'] = 120

# Rydberg to cm-1 and meV conversion factor
RyToCm  = 109691.40235
RyTomev = 13605.698066


# Load the harmonic dynamical matrix
harm_dyn = CC.Phonons.Phonons("dyn_harmonic", full_name = True)

# Prepare the sscha matrix by copying the harmonic one
start_dyn = harm_dyn.Copy()

# Temperature of the simulation
T = 0 #K

# The number of configurations
N = 4000

# The number of step
M = 20

# STEP SIZE
MOVE_STEP = 2e-3

# Generate the initial ensemble
ensemble = sscha.Ensemble.Ensemble(start_dyn, T)
ensemble.generate(N)


# Get energies and forces for the Harmonic Toy Model.
for i in range(N):
    energy, force = harm_dyn.get_energy_forces(ensemble.structures[i])
    ensemble.energies[i] = energy
    ensemble.forces[i, :,:] = force

# Prepare the symmetries at Gamma
qe_sym = CC.symmetries.QE_Symmetry(start_dyn.structure)
qe_sym.SetupQPoint(np.array( [0,0,0]), verbose = True)

# Extract the random move direction
nat = harm_dyn.structure.N_atoms
rand_vect = np.random.normal(size = (3 * nat, 3 * nat))
rand_vect /= np.sqrt( np.sum(rand_vect**2) )

# Impose that the random move must satisfy symmetries and translations
qe_sym.ImposeSumRule(rand_vect)
qe_sym.SymmetrizeDynQ(rand_vect, np.array([0,0,0]))


# Prepare the output to store information of the SSCHA free energy and gradient
print "----- MOVING DYN ------- "
print "# Step, free energy IS, error, free energy STOC, error, gradient IS, gradient STOC, Effective sample size"
steps = np.zeros(M)
grad_IS = np.zeros(M)
grad_IS_err = np.zeros(M)
free_IS = np.zeros(M)
free_IS_err = np.zeros(M)
grad_STOC = np.zeros(M)
grad_STOC_err = np.zeros(M)
free_STOC = np.zeros(M)
free_STOC_err = np.zeros(M)
kong_liu = np.zeros(M)

# Move the dynamical matrix
for i in range(M):
    
    # Compute free energy and its error
    free_IS[i], free_IS_err[i] = ensemble.get_free_energy(return_error = True)
    
    # Get the gradient of the free energy and its error
    g, g_err = ensemble.get_free_energy_gradient_respect_to_dyn()
    
    # Project the gradient along the moving direction
    grad_IS[i] = np.sum ( g * rand_vect )
    g_err = np.sum ( g_err * rand_vect )
        
    # Prepare a new ensemble to test also the importance sampling
    test_ensemble = sscha.Ensemble.Ensemble(ensemble.current_dyn.Copy(), T)
    test_ensemble.generate(N)
        
    # Use the harmonic toy model to compute energies and forces for this new ensemble
    for k in range(N):
        energy, force = harm_dyn.get_energy_forces(test_ensemble.structures[k])
        test_ensemble.energies[k] = energy
        test_ensemble.forces[k, :,:] = force
                            
    
    # Save the ensemble in the old sscha format
    #test_ensemble.save("data_dir", population = i+1)
    #test_ensemble.dyn_0.save_qe("dyn_pop%d_" % (i+1))
    
    #test_ensemble.update_weights(minim.dyn, T)

    # Compute the same quantities with the new ensemble
    free_STOC[i], free_STOC_err[i] = test_ensemble.get_free_energy(return_error = True)
    g, g_new_err = test_ensemble.get_free_energy_gradient_respect_to_dyn()
    
    # Project the gradient on the moving direction
    grad_STOC[i] = np.sum(g * rand_vect)
    g_new_err = np.sum(g_new_err * rand_vect)
    
    # Compute the Kong-Liu effective sample size
    kong_liu[i] = ensemble.get_effective_sample_size() / float(N)
    
    # Perform the step for the next calculation
    new_dyn = start_dyn.Copy()
    new_dyn.dynmats[0] += i * MOVE_STEP * rand_vect
    
    # Save the step 
    steps[i] = i * MOVE_STEP
    
    # Update the weights of the ensemble
    ensemble.update_weights(new_dyn, T)
    
    # Print on stdout all the info
    print "%5d%16.8e%16.8e%16.8e%16.8e%16.8e%16.8e%16.8e%16.8e%16.8e" % (steps[i] , 
                                                                         free_IS[i] , 
                                                                         free_IS_err[i],
                                                                         free_STOC[i], 
                                                                         free_STOC_err[i], 
                                                                         grad_IS[i], 
                                                                         grad_IS_err[i],
                                                                         grad_STOC[i], 
                                                                         grad_STOC_err[i], 
                                                                         kong_liu[i])



print ""
print "Data saved on 'python_sscha_minim.dat'"
# Save all the data stored
save_mat = np.transpose([steps,
                         free_IS*RyTomev, free_IS_err*RyTomev, 
                         free_STOC*RyTomev, free_STOC_err*RyTomev,
                         grad_IS, grad_STOC, kong_liu])
np.savetxt("python_sscha_minim.dat", save_mat, 
           header = "Steps|freeIS|err|freeSTOC|err|av_enIS|err|av_enSTOC|err|Escha|gradIS|gradSTOC|KL" )
                         
print ""
print "Plotting results"
# Plot the Free energy
plt.figure()
plt.errorbar(steps, free_IS * RyTomev, yerr = free_IS_err * RyTomev, label="IS")
plt.errorbar(steps, free_STOC * RyTomev, yerr = free_STOC_err * RyTomev, label = "STOC")
plt.xlabel("Step")
plt.ylabel("Free energy")
plt.title("Free energy test")
plt.legend()
plt.tight_layout()

# Plot the gradient (compare with the numerical difference of Free energy)
plt.figure()
plt.errorbar(steps, grad_IS, yerr=grad_IS_err, label="IS")
plt.errorbar(steps, grad_STOC, yerr=grad_STOC_err, label="STOC")
plt.plot(steps[:-1] + MOVE_STEP * .5, np.diff(free_IS) / MOVE_STEP, label="IS finite difference")
plt.plot(steps[:-1] + MOVE_STEP * .5, np.diff(free_STOC) / MOVE_STEP, label="STOC finite difference")
plt.xlabel("Step")
plt.ylabel("gradient along the direction")
plt.title("Gradient test")
plt.legend()
plt.tight_layout()

# Plot the Kong-Liu effective sample size
plt.figure()
plt.plot(steps, kong_liu)
plt.title("Kong-Liu effective sample size")
plt.xlabel("Step")
plt.ylabel("$N_{eff} / N$")
plt.tight_layout()
plt.show()