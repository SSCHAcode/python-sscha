# -*- coding: utf-8 -*-

"""
This program attempt the minimzation of only the wyckoff positions.
This is an explicit minimization, whitout passing through the
SchaMinimizer. It can be compared with a pure minimization using the old
sscha code.
"""
import numpy as np
import matplotlib.pyplot as plt

import cellconstructor as CC
import cellconstructor.Phonons
import sscha, sscha.Ensemble, sscha.SchaMinimizer

# Where are the data
DATA_DIR="../ensemble_data_test"
DYNPATH="../ensemble_data_test/dyn"
POPULATION=2
N_RANDOM=2000
T=0
EQ_ENERGY = -144.40680397
# number of steps
MAX_KA = 50
STEP_SIZE = 1

USE_PRECONDITIONING = True

# Start witht the interactive mode
plt.ion()


# Load the ensenble and the dynamical matrix
dynmat = CC.Phonons.Phonons(DYNPATH)
ens = sscha.Ensemble.Ensemble(dynmat, T)

nat = dynmat.structure.N_atoms

print "Loading the ensemble..."
ens.load(DATA_DIR, POPULATION, N_RANDOM)

print "Initializing symmetries..."
qe_sym = CC.symmetries.QE_Symmetry(dynmat.structure)
qe_sym.SetupQPoint( verbose = True ) # Get the symmetry in Gamma

prev_grad = None

# Get the preconditioner
precond = sscha.SchaMinimizer.GetStructPrecond(dynmat)
print "Best step size = ", 1 / np.max(np.real(np.linalg.eigvals(dynmat.dynmats[0]))) * CC.Phonons.BOHR_TO_ANGSTROM**2
    
# Perform the minimization
__fe__ = []
__gw__ = []
__kl__ = []

for ka in range(MAX_KA):
    # Compute the gradient
    grad = ens.get_average_forces(False)
    
    # Symmetrize the gradient 
    qe_sym.SymmetrizeVector(grad)
    
    
    # Print the free energy
    fe, fe_err = ens.get_free_energy(True)
    fe -= EQ_ENERGY
    gw = np.sqrt( np.sum(grad**2) )
    kl= ens.get_effective_sample_size()
    
    __fe__.append(fe)
    __gw__.append(gw)
    __kl__.append(kl)
    
    
    if USE_PRECONDITIONING:
        #STEP_SIZE = 2e-2
        grad = precond.dot(grad.reshape(3*nat)).reshape((nat, 3))
#    
#    if not prev_grad is None:
#        # Use the line minimization
#        y0 = np.sum(prev_grad * prev_grad)
#        y1 = np.sum(grad * prev_grad) 
#        
#        STEP_SIZE *= y0 / (y0 - y1) 
#        print "STEP SIZE: ", STEP_SIZE
#    prev_grad = grad
    
    # Apply the gradient in the structure
    dynmat.structure.coords += grad * STEP_SIZE
    
    # Upsate the ensemble
    ens.update_weights(dynmat, T)    
    
    
    print ""
    print " ---------- NEW MINIMIZATION STEP -----------"
    print "ka = ", ka
    print "Fe = ", fe * sscha.SchaMinimizer.__RyTomev__, " +- ", fe_err  * sscha.SchaMinimizer.__RyTomev__, " meV"
    print "gw = ", gw
    print "KL = ", kl

__fe__ = np.array(__fe__)


plt.figure()
plt.title("Free energy")
plt.ylabel("Free energy [meV]")
plt.xlabel("Step")
plt.plot(__fe__ * sscha.SchaMinimizer.__RyTomev__)
plt.tight_layout()

plt.figure()
plt.title("Gradient")
plt.ylabel("gw")
plt.xlabel("Step")
plt.plot(__gw__)
plt.tight_layout()

plt.figure()
plt.title("Kong-Liu sample size")
plt.ylabel("KL")
plt.xlabel("Step")
plt.plot(__kl__)
plt.tight_layout()

plt.show()