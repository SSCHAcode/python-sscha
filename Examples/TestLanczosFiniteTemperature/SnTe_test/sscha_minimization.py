from __future__ import print_function
from __future__ import division

import cellconstructor as CC
import cellconstructor.Phonons

import spglib

# Import the modules of the force field
import fforces as ff
import fforces.Calculator

# Import the modules to run the sscha
import sscha, sscha.Ensemble, sscha.SchaMinimizer
import sscha.Relax, sscha.Utilities

import numpy as np
import sys,os

# Load the dynamical matrix for the force field
ff_dyn = CC.Phonons.Phonons("ffield_dynq", 3)

# Setup the forcefield with the correct parameters
ff_calculator = ff.Calculator.ToyModelCalculator(ff_dyn)
ff_calculator.type_cal = "pbtex"
ff_calculator.p3 = 0.036475
ff_calculator.p4 = -0.022
ff_calculator.p4x = -0.014

T = 250

# Load the original dynamical matrix
sscha_dyn = ff_dyn.Copy()
sscha_dyn.ForcePositiveDefinite()
sscha_dyn.Symmetrize()

# Setup the minimization
print("Generate the ensemble...")
ens = sscha.Ensemble.Ensemble(sscha_dyn, 250, sscha_dyn.GetSupercell())
minim = sscha.SchaMinimizer.SSCHA_Minimizer(ens)
minim.min_step_dyn = 0.01
minim.root_representation = "root2"
#minim.precond_dyn = False
minim.min_step_struc = 0.1

# Perform the automatic relaxation
print("Prepare relaxation...")
relax = sscha.Relax.SSCHA(minim, ff_calculator, N_configs = 10000, max_pop = 3)
print("Start the relaxation...")
relax.relax()

print("New relaxation")
relax.max_pop = 12
relax.minim.min_step_dyn = 0.8
relax.relax()

print("Saving results in ensemble (population 1) and SnTe_sscha")
relax.minim.dyn.save_qe("SnTe_sscha")
relax.minim.ensemble.save_bin("ensemble", 1)

relax.minim.plot_results()


