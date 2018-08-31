"""
This script test whether the wyckoff preconditioner is actually the inverse of the
Dynamical matrix or not
"""
import numpy as np

import cellconstructor as CC
import cellconstructor.Phonons

import sscha, sscha.SchaMinimizer

# Load the dynamical matrix
dyn = CC.Phonons.Phonons("dyn")

# Get the preconditioner
prec = sscha.SchaMinimizer.GetStructPrecond(dyn)

# Multiply the two matrix
identity = prec.dot(np.real(dyn.dynmats[0]))

# Check whether the eigenvalues of the identity are close to 1
eigvals = np.linalg.eigvals(identity)
eigvals = eigvals[ eigvals > 0.1 ]
print np.sqrt(np.sum( (eigvals- 1)**2))


