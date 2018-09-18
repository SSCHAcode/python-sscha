# -*- coding: utf-8 -*-
import numpy as np
import time

import cellconstructor as CC
import cellconstructor.Phonons

import sscha, sscha.Ensemble
import sscha.SchaMinimizer

import SCHAModules


"""
This example compares the result of the fortran vs python in computing the gradient
"""

N = 1000

# Load the dynamical matrix and the ensemble
dynmat = CC.Phonons.Phonons("../ensemble_data_test/dyn")
ens = sscha.Ensemble.Ensemble(dynmat, 0)

print "Loading the ensemble..."
t1 = time.time()
ens.load("../ensemble_data_test", 2, N)
t2 = time.time()
print "Time elapsed to load the ensemble ", t2 -t1

nat = dynmat.structure.N_atoms



# Get the gradient from python
t1 = time.time()
grad_py, err_grad_py = ens.get_fc_from_self_consistency(True, True)
grad_py = sscha.SchaMinimizer.ApplyLambdaTensor(dynmat, grad_py, 0)
t2 = time.time()
print "Time to compute the python gradient:", t2 - t1

u_disp_f = np.zeros((N, nat,3), dtype = np.float64, order = "F")
eforces = np.zeros((N, nat,3), dtype = np.float64, order = "F")

eforces[:,:,:] = ens.forces - ens.sscha_forces
fforces = np.zeros((N, 3*nat), dtype = np.float64, order = "F")
for i in range(N):
    u_disp_f[i, :,:] = np.reshape(ens.u_disps[i,:], (nat, 3))
    fforces[i,:] = np.reshape(ens.forces[i, :,:] - ens.sscha_forces[i,:,:], (nat*3))


w, pols = dynmat.DyagDinQ(0)
pols = np.real(pols)
trans = CC.Methods.get_translations(pols, dynmat.structure.get_masses_array())

# Convert in the correct units
w /= 2
eforces /= 2 
rho = np.array(ens.rho, dtype = np.float64)
ityp = dynmat.structure.get_ityp() + 1
mass = np.array(dynmat.structure.masses.values())
mass *= 2
print "ITYP:", ityp
print "MASS:", mass



# Get the gradient from fortran
t1 = time.time()
grad_f, err_grad_f = SCHAModules.get_gradient_supercell(rho, u_disp_f, eforces, w, pols, trans.astype(np.intc), 
                                                        0, mass, ityp, "err_yesrho",
                                                        N, nat, 3*nat, preconditioned = 0)
#grad_f *= -2
t2 = time.time()

print "Elapsed time to compute the fortran gradient: ", t2 -t1

# Lets try also the manual computation
uf2 = np.einsum("i, ij, ik", ens.rho, ens.u_disps , fforces) / np.sum(ens.rho)
np.savetxt("uf_py.dat", uf2, header="UF computed with python")
ups2 = dynmat.GetUpsilonMatrix(T = 0.0)
new_grad = ups2.dot(uf2)
new_grad = 0.5*(new_grad + np.transpose(new_grad))


new_dyn = dynmat.Copy()
new_dyn.dynmats[0][:,:] = grad_py
new_dyn.save_qe("grad_py")
new_dyn.dynmats[0][:,:] = grad_f
new_dyn.save_qe("grad_f")
new_dyn.dynmats[0][:,:] = new_grad
new_dyn.save_qe("grad_new")

# Difference between the two
diff = grad_py - grad_f
print "Difference between the two:", np.sqrt(np.einsum("ab, ba", diff, diff))

diff = new_grad - grad_py
print "Difference between python two:", np.sqrt(np.einsum("ab, ba", diff, diff))

diff = new_grad - grad_f
print "Difference between fortran two:", np.sqrt(np.einsum("ab, ba", diff, diff))
 
