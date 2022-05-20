# Import the sscha code
import sscha, sscha.Ensemble, sscha.SchaMinimizer, sscha.Relax, sscha.Utilities

# Import the cellconstructor library to manage phonons
import cellconstructor as CC, cellconstructor.Phonons
import cellconstructor.Structure, cellconstructor.calculators

# Import the force field of Gold
import ase, ase.calculators
from ase.calculators.emt import EMT

# Import numerical and general pourpouse libraries
import numpy as np, matplotlib.pyplot as plt
import sys, os


"""
You need first to run the
get_gold_free_energy.py

Here we use NPT simulation to compute the gold thermal expansion.
"""

# Define the temperature range (in K)
T_START = 300
T_END = 1000
DT = 50

N_CONFIGS = 50
MAX_ITERATIONS = 10

# Import the gold force field
calculator = EMT()

# Import the starting dynamical matrix (final result of get_gold_free_energy.py)
dyn = CC.Phonons.Phonons("sscha_T300_dyn", nqirr = 13)

# Create the directory on which to store the output
DIRECTORY = "thermal_expansion"
if not os.path.exists(DIRECTORY):
    os.makedirs("thermal_expansion")

# We cycle over several temperatures
t = T_START


volumes = []
temperatures = []
while t <= T_END:
    # Change the temperature
    ensemble = sscha.Ensemble.Ensemble(dyn, t)
    minim = sscha.SchaMinimizer.SSCHA_Minimizer(ensemble)
    minim.set_minimization_step(0.1)

    relax = sscha.Relax.SSCHA(minim, calculator, N_configs = N_CONFIGS,
                              max_pop = MAX_ITERATIONS)

    # Setup the I/O
    ioinfo = sscha.Utilities.IOInfo()
    ioinfo.SetupSaving( os.path.join(DIRECTORY, "minim_t{}".format(t)))
    relax.setup_custom_functions( custom_function_post = ioinfo.CFP_SaveAll)


    # Run the NPT simulation
    relax.vc_relax(target_press = 0)

    # Save the volume and temperature
    volumes.append(relax.minim.dyn.structure.get_volume())
    temperatures.append(t)

    # Start the next simulation from the results
    relax.minim.dyn.save_qe( os.path.join(DIRECTORY, "sscha_T{}_dyn".format(t)))
    dyn = relax.minim.dyn
    relax.minim.finalize()
    
    # Update the temperature
    t += DT

# Save the thermal expansion
np.savetxt(os.path.join(DIRECTORY, "thermal_expansion.dat"),
           np.transpose([temperatures, volumes]),
           header = "Temperature [K]; Volume [A^3]")

    




