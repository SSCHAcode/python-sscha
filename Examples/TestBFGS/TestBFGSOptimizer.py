# -*- coding: utf-8 -*-

"""
Here we test the BFGS optimizer
using a custom stress tensorm function
"""
import numpy as np
import sscha, sscha.Optimizer

I = np.eye(3)
def free_energy(uc):
    """
    Free energy sample function 
    """
    return np.trace( (uc - I).dot( np.transpose(uc - I)))

def stress(uc):
    """
    Compute the stress tensor of the given free energy
    """
    Omega = np.linalg.det(uc)
    return - 2 * uc.dot(uc - I) / Omega


# Initialize a random unit cell
uc = np.random.uniform(size=(3,3))

N_STEPS = 10
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
    
    