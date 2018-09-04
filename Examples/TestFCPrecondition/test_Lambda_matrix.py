# -*- coding: utf-8 -*-

"""
This is a test to check if the Lambda matrix (the hessian of the free energy
with respect to the force-constant matrix) is computed correctly by
the methods.

First of all, we check that the subsequent application of its inversion
is the identity.
"""

import numpy as np
import cellconstructor as CC
import cellconstructor.Phonons

import sscha, sscha.SchaMinimizer

DYNMAT = "../ensemble_data_test/dyn"
ENSEMBLE = "../ensemble_data_test/"
N = 1000

# Load the dynamical matrix
dyn = CC.Phonons.Phonons(DYNMAT)

print "Loading the ensemble..."
ens = sscha.Ensemble.Ensemble(dyn, 0)
ens.load(ENSEMBLE, 2, 1000)

print "Computing the preconditioned gradient..."
# Get the gradient
gc = ens.get_fc_from_self_consistency(True, False)

# Fix the translations
qe_sym = CC.symmetries.QE_Symmetry(dyn.structure)
qe_sym.SetupQPoint()

# Symmetrize the gradient
qe_sym.SymmetrizeDynQ(gc, np.array([0,0,0]))
qe_sym.ImposeSumRule(gc)

#fc = dyn.dynmats[0]

#qe_sym.ImposeSumRule(fc)

# Apply the lambda matrix to dyn
print "Applying the precondition to dyn..."

lambda_dyn = sscha.SchaMinimizer.ApplyFCPrecond(dyn, gc)

print "Applying the Lambda matrix..."
new_dyn = sscha.SchaMinimizer.ApplyLambdaTensor(dyn, lambda_dyn)

print ""

# We want to compare it with the application of
# the Identity minus the translations
w, pols = dyn.DyagDinQ(0)
m = dyn.structure.get_masses_array()
trans = ~ CC.Methods.get_translations(pols, m)
w = w[trans]
pols = pols[:, trans]
nat = dyn.structure.N_atoms
pols1 = pols.copy()
pols2 = pols.copy()
for i in range(nat):
    pols1[3*i : 3*i + 3, :] /= np.sqrt(m[i])
    pols2[3*i: 3*i + 3, :] *= np.sqrt(m[i])
    
print "Computing the matrix with the correct transformation..."
new_dyn2 = np.einsum("ah, bk, ch, dk, cd", pols1, pols1, pols2, pols2, gc)

print ""


print "Distance between the matrices:", np.sqrt(np.sum( (np.real(new_dyn - new_dyn2))**2))

print "Distance before and after the double preondition:", np.sqrt(np.sum( (new_dyn - gc)**2))

# Impose symmetry on the new matrix
#qe_sym.SymmetrizeDynQ(new_dyn, np.array([0,0,0]))
qe_sym.ImposeSumRule(new_dyn)
print "Distance after the application of the sum rule:", np.sqrt(np.sum( (new_dyn - gc)**2))
