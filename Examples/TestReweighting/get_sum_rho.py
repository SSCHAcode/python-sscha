# Cross compatibility between python 2 and python 3
from __future__ import print_function
from __future__ import division

import numpy as np

import sscha, sscha.Ensemble
import cellconstructor as CC
import cellconstructor.Phonons

N_RANDOM = 10
POPULATION = 2

TESTED_VALUE = 1.7445855558823957 


dyn = CC.Phonons.Phonons("super_final")
dyn_start = CC.Phonons.Phonons("../ensemble_data_test/dyn")

ens = sscha.Ensemble.Ensemble(dyn_start, 0)
ens.load("../ensemble_data_test/", POPULATION, N_RANDOM)

# Store the old displacements
old_disp = ens.u_disps.copy()
ens.update_weights(dyn, 100)
new_disp = ens.u_disps.copy()

# print("OLD Disp | NEW Disp")
# print("\n".join("{:8d} -> {}".format(i, old_disp[i,:] - new_disp[i,:]) for i in range (N_RANDOM)))

# print()
# print("Real disp:")
# print(dyn_start.structure.coords - ens.current_dyn.structure.coords)


print( "The sum over rho:", np.sum(ens.rho) / ens.N)
print( "The tested value is:", TESTED_VALUE)

