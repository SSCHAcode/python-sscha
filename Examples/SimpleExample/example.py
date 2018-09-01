# -*- coding: utf-8 -*-

"""
This example shows a working minimization of a unit cell of ICE XI,
a complex molecular crystal with a lot of free parameters.
Here we show how to perform the minimization using the ensemble and force
calculation present in the data_dir.
"""

import cellconstructor as CC
import cellconstructor.Phonons

import sscha
import sscha.Ensemble
import sscha.SchaMinimizer

# Load the dynamical matrix
dyn1 = CC.Phonons.Phonons("dyn")

# Load the ensemble
ens = sscha.Ensemble.Ensemble(dyn1, 0)
ens.load("data", 2, 10000) # This was a second population (2)


# Setup the minimizer
minim = sscha.SchaMinimizer.SSCHA_Minimizer(ens)
minim.min_step_struc = 1e-4
minim.min_step_dyn = 0
minim.meaningful_factor = 1e-4
minim.eq_energy = -144.40680397
minim.precond_wyck = False

# Run the minimization
minim.init()
minim.run()
