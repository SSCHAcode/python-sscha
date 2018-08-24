 # -*- coding: utf-8 -*-

"""
TEST MINIMIZATION

This example generate an input harmonic matrix that uses as a toy model
to compute forces and energies.
Then the dynamical matrix is moved from its equilibrium position and a SSCHA
minimization is started.

The minimum should be reached when the dynamical matrix is equal to the
harmonic one.

This example shows a trivial minimization and the use of an harmonic
Toy Model.
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


# Load the harmonic dynamical matrix
harm_dyn = CC.Phonons.Phonons("dyn_harmonic", full_name = True)

# Prepare the sscha matrix by copying the harmonic one
start_dyn = harm_dyn.Copy()

# Temperature of the simulation
T = 0 #K

# The number of configurations
N = 4000

# The distance to move the initial dynamical matrix respect to the minimum
DISTANCE = 1e-2

LAMBDA_A = 0.70486658940500413e-1

# If this flag is true, the ensemble is loaded from an existing one
LOAD_ENSEMBLE = True

# Prepare the starting dynamcal matrix by summing a random vector
# That respect all the symmetries of the given structure

# Prepare the symmetries at Gamma
qe_sym = CC.symmetries.QE_Symmetry(start_dyn.structure)
qe_sym.SetupQPoint(np.array( [0,0,0]), verbose = True)



if not LOAD_ENSEMBLE:
    # Extract the random move direction
    nat = harm_dyn.structure.N_atoms
    rand_vect = np.random.normal(size = (3 * nat, 3 * nat))
    rand_vect /= np.sqrt( np.sum(rand_vect**2) )
    
    # Impose that the random move must satisfy symmetries and translations
    qe_sym.ImposeSumRule(rand_vect)
    qe_sym.SymmetrizeDynQ(rand_vect, np.array([0,0,0]))
    
    # Move the starting dynamical matrix respect to the origin
    start_dyn.dynmats[0] += rand_vect * DISTANCE
    
    
    # Prepare the initial ensemble
    ensemble = sscha.Ensemble.Ensemble(start_dyn, T)
    ensemble.generate(N)
    
    # Get energies and forces for the Harmonic Toy Model.
    for i in range(N):
        energy, force = harm_dyn.get_energy_forces(ensemble.structures[i])
        ensemble.energies[i] = energy
        ensemble.forces[i, :,:] = force
        
        
    # Save the ensemble to do the SSCHA test
    ensemble.save("data_dir", population = 1)
    ensemble.dyn_0.save_qe("dyn")
else :
    # Load the ensemble
    start_dyn = CC.Phonons.Phonons("dyn", nqirr = 1)
    ensemble = sscha.Ensemble.Ensemble(start_dyn, T)
    ensemble.load("data_dir", population = 1, N = N)
    

# Start the minimization
minim = sscha.SchaMinimizer.SSCHA_Minimizer(ensemble)
minim.min_step_dyn = LAMBDA_A

# Define the custom function to print both the gradient and the dynamical 
# matrix projected in the generator
fc = np.zeros((6,6), dtype = np.float64)
fc[0,0] = 1

# Apply the sum rule and the symmetries
qe_sym.ImposeSumRule(fc)
qe_sym.SymmetrizeDynQ(fc, np.array([0,0,0]))

# Orthonormalize
fc /= np.sqrt(np.trace(fc.dot(fc)))

# Custom function
def print_dyn_grad_projection(sscha_minim):
    # Project the dyn
    dyn_coeff = np.trace( fc.dot(sscha_minim.dyn.dynmats[0]) )
    gc_coeff = np.trace( fc.dot(sscha_minim.prev_grad))
    
    print "Dyn_coeff [Ha/bohr^2]: ", dyn_coeff / 2
    print "GC_coeff [Ha/bohr^3]: ", gc_coeff / 2

# Initialize the minimization
minim.init()

# Cycle over the ensebles
running = True
while running:
    print "Running the minimization..."
    minim.run(2, print_dyn_grad_projection)
    
    # If the minimization is converged, exit.
    if minim.is_converged():
        running = False
        print "The minimization is converged"
    else:
        # Extract the new ensemble
        print "We are out of the statistical sampling."
        print "Generating a new ensemble..."
        new_ensemble = sscha.Ensemble.Ensemble(minim.dyn, T)
        new_ensemble.generate(N)
        
        break
        
        # Compute energies and forces with the harmonic Toy Model
        for i in range(N):
            energy, force = harm_dyn.get_energy_forces(ensemble.structures[i])
            new_ensemble.energies[i] = energy
            new_ensemble.forces[i, :,:] = force
        
        # Update the minimization ensemble and start again
        minim.ensemble = new_ensemble
        
    # Only one ensemble
    running = False
                

# The minimization is done.
# Now plot the result
print "Minimization done."
print "Plotting the results"
minim.plot_results("minim_python_data.dat")

# Save the final dynamical matrix
minim.dyn.save_qe("dyn_converged")