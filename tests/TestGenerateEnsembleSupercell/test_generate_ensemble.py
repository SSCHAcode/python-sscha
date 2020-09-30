# -*- coding: utf-8 -*-
from __future__ import print_function
from __future__ import division

import sys, os
import shutil
import cellconstructor as CC
import cellconstructor.Phonons

import sscha, sscha.Ensemble

def test_generate_ensemble_sup():
    total_path = os.path.dirname(os.path.abspath(__file__))
    os.chdir(total_path)

    T = 0
    NQIRR = 3
    SUPERCELL = (2,1,2)
    N_RANDOM = 10

    DATA_DIR = "_data_tmp_"
    POPULATION = 1

    # Load the dynamical matrix
    dyn = CC.Phonons.Phonons("dyn", NQIRR)

    # Generate the ensemble
    ens = sscha.Ensemble.Ensemble(dyn, T, SUPERCELL)
    ens.generate(N_RANDOM)

    # Save the ensemble
    if not os.path.isdir(DATA_DIR):
        os.mkdir(DATA_DIR)
    ens.save(DATA_DIR, POPULATION)
    shutil.rmtree(DATA_DIR, ignore_errors= True)

if __name__ == "__main__":
    test_generate_ensemble_sup()
