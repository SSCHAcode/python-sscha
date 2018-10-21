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

    grad1 = 2 * (uc - I)

    du_deps = np.zeros( (3,3,3,3))

    for i in range(3):
        for j in range(3):
            for x in range(3):
                for y in range(3):
                    du_deps[i,j,x,y] = 0.5 * (uc[i,y] * I[x,j] + uc[i,x] * I[y,j])
    
    dF_deps = np.einsum("ab, abcd", grad1, du_deps)

    return - dF_deps / Omega


# Initialize a random unit cell
uc = np.random.uniform(size=(3,3))

N_STEPS = 20
#opt = sscha.Optimizer.BFGS_UC()
opt = sscha.Optimizer.UC_OPTIMIZER()
opt.alpha = 0.1
for i in range(N_STEPS):
    print "---------- STEP %d ----------" % i
    print "UC:"
    print uc
    print "free energy:", free_energy(uc)
    print "Stress:"
    print stress(uc)
    print "ALPHA:", opt.alpha
    
    # Step
    opt.UpdateCell(uc, stress(uc))
    
    