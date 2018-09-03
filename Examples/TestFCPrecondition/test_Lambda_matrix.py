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

# Load the dynamical matrix
dyn = CC.Phonons.Phonons(DYNMAT)

# Fix the translations
qe_sym = CC.symmetries.QE_Symmetry(dyn.structure)
qe_sym.SetupQPoint()

fc = dyn.dynmats[0]

qe_sym.ImposeSumRule(fc)

# Apply the lambda matrix to dyn
print "Applying the precondition to dyn..."

lambda_dyn = sscha.SchaMinimizer.ApplyFCPrecond(dyn, fc)

print "Applying the Lambda matrix..."
new_dyn = sscha.SchaMinimizer.ApplyLambdaTensor(dyn, lambda_dyn)

print ""

qe_sym.ImposeSumRule(new_dyn)
qe_sym.ImposeSumRule(fc)

w0 = np.linalg.eigvals(new_dyn)
w1 = np.linalg.eigvals(fc)

print np.sort(np.abs(w0 - w1))

# Check if the original matrix is equal to the new one
dist = np.sum( (w0 -w1)**2 )
dist = np.sqrt(dist)

print "Distance between the two matrices (should be 0):", dist