# -*- coding: utf-8 -*-

"""
Here we test the BFGS optimizer
using a custom stress tensorm function
"""
import numpy as np
import sscha, sscha.Optimizer


def free_energy(uc):
    """
    Free energy sample function 
    """
    return np.trace((uc.dot(np.diag([1,10,100])) - eye(3)).dot(uc.dot(np.diag([1,10,100])) - eye(3)))
def stress(uc):
    """
    Compute the stress tensor of the given free energy
    """
    return 2 * (uc.dot(np.diag([1,10,100])) - eye(3))


# Initialize a random unit cell
uc = np.random.uniform(size=(3,3))

N_STEPS = 30
opt = sscha.Optimizer.BFGS_UC()
for i in range(N_STEPS):
    print "---------- STEP %d ----------" % i
    print "UC:"
    print uc
    print "free energy:", free_energy(uc)
    print "Stress:"
    print stress(uc)
    
    # Step
    opt.UpdateCell(uc, stress(uc))
    
    