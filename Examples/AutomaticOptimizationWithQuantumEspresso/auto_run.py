# -*- coding: utf-8 -*-

"""
This code set up a quantum espresso ASE calculator
to be used as an engine for force and energy calculation.
In this way all the minimization is automatized by the code
"""

import cellconstructor as CC
import cellconstructor.Phonons

import sscha, sscha.Ensemble
import sscha.SchaMinimizer

# The Quantum ESPRESSO calculator
from ase.calculators.espresso import Espresso

# Adjust the pseudo you want to use
# by default ase will look for them in $HOME/espresso/pseudo
pseudo = {"H": "H.pbe-rrkjus_psl.0.1.UPF",
          "O": "O.pbe-n-rrkjus_psl.0.1.UPF"}


# Setup some info on the calculator (the wavefunction and density cutoff)
# Note they are super below convergence
# Good values would be (45, 45*8) for the cutoffs, and (3,3,2) for the k points
esp_info = {"ecutwfc" : 25, "ecutrho" : 25*4, "disk_io" : "none"}
K_POINTS = (1,1,1) # Only gamma calculation

# Setup the quantum espresso command for run (parallel on 4 processors here)
cmd = "/usr/bin/mpirun -np 2 $HOME/Downloads/QuantumESPRESSO/qe-6.3/bin/pw.x -in PREFIX.pwi > PREFIX.pwo"


# Setup the calculator
calc = Espresso(pseudopotentials = pseudo, tstress = True ,
                tprnfor = True, kpts = K_POINTS, input_data = esp_info,
                command = cmd)




#Now load the initial dynamical matrix
dyn_start = CC.Phonons.Phonons("../ensemble_data_test/dyn")
T = 0 # Temperature of the calculation
N_RAND = 500 # number of random configurations

# Setup the ensemble
ens = sscha.Ensemble.Ensemble(dyn_start, T)

# Generate the ensemble
ens.generate(N_RAND)

# Get forces and energy using the quantum espresso calculator
ens.get_energy_forces(calc, True)



# Setup the minimizer
minim = sscha.SchaMinimizer.SSCHA_Minimizer(ens, root_representation = "root4")
minim.precond_dyn = False
minim.min_step_dyn = 0.5

# Run the minimization
minim.init()
minim.run()
minim.finalize(2)
minim.plot_results()



