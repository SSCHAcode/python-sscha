"""
Test the multiplication for the lambda matrix made by python and fortran.

They should match
"""
import numpy as np
import time

import cellconstructor as CC
import cellconstructor.Phonons
import sscha, sscha.Ensemble, sscha.SchaMinimizer
import SCHAModules

FILDYN="../ensemble_data_test/dyn"
T = 0


# Load the dynamical matrix
dyn = CC.Phonons.Phonons(FILDYN)

# Extract the force constant matrix
fc = dyn.dynmats[0]

t1 = time.time()
new_py = sscha.SchaMinimizer.ApplyLambdaTensor(dyn, fc, T)
t2 = time.time()
print "Time for the python Lambda application:", t2 - t1

t1 = time.time()
# Do the same with the fortran matrix
w, pols = dyn.DyagDinQ(0)
w /= 2 # Ha convertion
trans = CC.Methods.get_translations(pols, dyn.structure.get_masses_array())
mass = np.array(dyn.structure.masses.values())
mass *= 2
ityp = dyn.structure.get_ityp() +1 # Fortran index conversion

fc /= 2 # Ha conversion

new_f = SCHAModules.multiply_lambda_tensor(w, pols, trans.astype(np.intc), mass,
                                           ityp, T, fc, 0)
t2 = time.time()
print "Time for the fortran Lambda application:", t2 - t1


np.savetxt("LambdaPhi_py.dat", new_py, header="Lambda . Phi in python")
np.savetxt("LambdaPhi_f.dat", new_f, header="Lambda . Phi in fortran")

diff = new_py - new_f
print "Distance between them:", np.sqrt(np.einsum("ab,ba", diff, diff))




