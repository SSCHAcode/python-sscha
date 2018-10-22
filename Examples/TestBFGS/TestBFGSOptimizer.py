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

    gradient = -2 *np.transpose(uc).dot(uc - I) / Omega

    return gradient


# Initialize a random unit cell
#uc = np.random.uniform(size=(3,3))
d = np.array([5,10,20], dtype = np.float64)
uc = np.diag( d)
#v1 = uc[0,:] + uc[1,:]
#v2 = uc[0,:] - uc[1,:]
#uc[0,:] = v1
#uc[1,:] = v2

N_STEPS = 5
opt = sscha.Optimizer.BFGS_UC()
#opt = sscha.Optimizer.UC_OPTIMIZER()
opt.alpha = 1
for i in range(N_STEPS):
    print "---------- STEP %d ----------" % i
    print "UC:"
    print uc
    print "free energy:", free_energy(uc)
    print "Stress:"
    print stress(uc)
    print "ALPHA:", opt.alpha

    print "TEST GRAD:"
    print "Simple:"
    print 2 * (uc -I)
    print "Complex:"
    g = 2 * (uc - I)
    new_g = 0.5 * (g  + np.transpose(np.linalg.inv(uc)).dot(np.transpose(g).dot(uc)))
    print new_g 
    # Step
    opt.UpdateCell(uc, stress(uc))
    
    
