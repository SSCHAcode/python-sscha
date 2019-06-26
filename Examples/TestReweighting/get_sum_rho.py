# Cross compatibility between python 2 and python 3
from __future__ import print_function
from __future__ import division

import numpy as np

import sscha, sscha.Ensemble
import cellconstructor as CC
import cellconstructor.Phonons

N_RANDOM = 10
POPULATION = 2

dyn = CC.Phonons.Phonons("final_dyn")
dyn_start = CC.Phonons.Phonons("../ensemble_data_test/dyn")

ens = sscha.Ensemble.Ensemble(dyn_start, 0)
ens.load("../ensemble_data_test/", POPULATION, N_RANDOM)
ens.update_weights(dyn, 0)

print( "The sum over rho:", np.sum(ens.rho) / ens.N)

