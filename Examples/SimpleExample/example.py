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

# Input info
DYNPATH="../ensemble_data_test/dyn"
DATADIR="../ensemble_data_test"
POPULATION=2
NRANDOM=1000
T=0

# Load the dynamical matrix
dyn1 = CC.Phonons.Phonons(DYNPATH)

# Load the ensemble
ens = sscha.Ensemble.Ensemble(dyn1, T)
ens.load(DATADIR, POPULATION, NRANDOM, verbose = True) # This was a second population (2)

# Setup the minimizer
minim = sscha.SchaMinimizer.SSCHA_Minimizer(ens)
minim.min_step_struc = 0
minim.min_step_dyn = 1
minim.meaningful_factor = 1e-7
minim.eq_energy = -144.40680397
minim.root_representation = "root4"
minim.precond_dyn = False
#minim.precond_wyck = False
minim.max_ka = 40
#minim.fake_precond = False
#minim.precond_dyn = False

minim.gradi_op = "all"

# Run the minimization
minim.init()
minim.run()
minim.finalize()        
# Plot the results
minim.plot_results()
