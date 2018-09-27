# -*- coding: utf-8 -*-

import cellconstructor as CC
import cellconstructor.Phonons

import sscha, sscha.Ensemble

T = 0
NQIRR = 3
SUPERCELL = (2,1,2)
N_RANDOM = 1000

DATA_DIR = "data"
POPULATION = 1

# Load the dynamical matrix
dyn = CC.Phonons.Phonons("dyn", NQIRR)

# Generate the ensemble
print "Generating the ensemble ..."
ens = sscha.Ensemble.Ensemble(dyn, T, SUPERCELL)
ens.generate(N_RANDOM)
print "Done."

# Save the ensemble
ens.save(DATA_DIR, POPULATION)
