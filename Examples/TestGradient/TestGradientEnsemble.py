# -*- coding: utf-8 -*-

"""
This script test the gradient usnig the ensemble_data_test
"""
import time
import numpy as np
import matplotlib.pyplot as plt

import cellconstructor as CC
import cellconstructor.Phonons

import sscha, sscha.Ensemble
import sscha.SchaMinimizer

DATADIR="../ensemble_data_test"
POPULATION = 2
N = 1000
T0 = 0
T = 0
STEP_SIZE = 2e-4
M_STEPS = 10
EQ_ENERGY = -144.40680397

# Seed the generators to use always the same random matrix
np.random.seed(0)

# Load the dynamical matrix
dyn0 = CC.Phonons.Phonons(DATADIR + "/dyn", nqirr = 1)

# Load the ensemble
print "Loading the ensembe..."
t = time.time()
ens = sscha.Ensemble.Ensemble(dyn0.Copy(), T0)
ens.load(DATADIR, POPULATION, N)
t_new = time.time()

print "Time ellapsed (seconds):", t_new - t

# Update the esnemble at the correct temperature if required
if T != T0:
    ens.update_weights(dyn0, T)

# Choose the random direction for the dynamical matrix
print "Initializing symmetries..."
qe_sym = CC.symmetries.QE_Symmetry(dyn0.structure)
qe_sym.SetupQPoint()

# Extract the random matrix
rand_mat = np.random.uniform(size = np.shape(dyn0.dynmats[0]))

#Apply the symmetry to the random matrix
qe_sym.SymmetrizeDynQ(rand_mat, dyn0.q_tot[0])
qe_sym.ImposeSumRule(rand_mat)

gc = []
gc_old = []
free = []
kl = []

print "Start the cycle:"
t = time.time()
for ka in range(M_STEPS):
    #g = ens.get_preconditioned_gradient()
    
    # Get also the old good gradient
    g_old = ens.get_fc_from_self_consistency()
    
    # Symmetrize the gradient
    #qe_sym.SymmetrizeFCQ(g, dyn0.q_stars)
    qe_sym.SymmetrizeDynQ(g_old, dyn0.q_tot[0])
    qe_sym.ImposeSumRule(g_old)
    
    # Multiply dot the lambda matrix
    #g[0, :,:] = sscha.SchaMinimizer.ApplyLambdaTensor(dyn0, g[0,:,:], T)
    g_old = sscha.SchaMinimizer.ApplyLambdaTensor(dyn0, g_old, T)
    
    # Symmetrize the gradient
    #qe_sym.SymmetrizeFCQ(g, dyn0.q_stars)
    qe_sym.SymmetrizeDynQ(g_old, dyn0.q_tot[0])
    qe_sym.ImposeSumRule(g_old)
    
    # Project into the direction
    #gc.append( np.einsum("ab, ba", g[0,:,:], rand_mat))
    gc_old.append( np.einsum("ab, ba", g_old, rand_mat))
    kl.append(ens.get_effective_sample_size())
    
    # Get the free energy
    free.append(ens.get_free_energy() - EQ_ENERGY)
    
    # Update the ensemble
    dyn0.dynmats[0] += rand_mat * STEP_SIZE
    
    ens.update_weights(dyn0, T)
    
    print " # ------ STEP ------ "
    print "  ka = ", ka
    print "  Free energy = ", free[-1] * sscha.SchaMinimizer.__RyTomev__
    print "  gc = ", gc_old[-1]  * sscha.SchaMinimizer.__RyTomev__
    print "  kl = ", kl[-1]
    t_new = time.time()
    print "  time = ", t_new - t, " s"
    t = t_new
    print ""
    
# Convert in numpy array
free = np.array(free)
#gc = np.array(free)
gc_old = np.array(gc_old)
kl = np.array(kl)

# Plot the result
print " Plotting the results..."
x_axis = np.arange(M_STEPS) * STEP_SIZE

plt.figure()
plt.title("Effective sample size")
plt.xlabel("Step size")
plt.ylabel("Kong Liu")
plt.plot(x_axis, kl)
plt.tight_layout()

plt.figure()
plt.title("Free energy")
plt.xlabel("Step size")
plt.ylabel("Free energy [meV]")
plt.plot(x_axis, free* sscha.SchaMinimizer.__RyTomev__)
plt.tight_layout()



plt.figure()
plt.title("Gradient")
#plt.plot(x_axis, gc, label="Real gradient")
plt.plot(x_axis, -np.gradient(free) / STEP_SIZE, label="Finite difference")
plt.plot(x_axis, gc_old, label="Gradient")
plt.xlabel("Step size")
plt.ylabel("Gradient")
plt.legend()
plt.tight_layout()

print "Done"
plt.show()
