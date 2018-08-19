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

# The distance to move the initial dynamical matrix respect to the minimum
DISTANCE = 1e-2

# Prepare the starting dynamcal matrix by summing a random vector
# That respect all the symmetries of the given structure

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


# Start the minimization
minim = sscha.SCHAMinimizer.SSCHA_Minimizer(ensemble)
minim.run()




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