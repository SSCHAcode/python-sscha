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
Here we load the primitive cell of Gold from a cif file.
And we use CellConstructor to compute phonons from finite differences.
The phonons are computed on a q-mesh 4x4x4
"""

gold_structure = CC.Structure.Structure()
gold_structure.read_generic_file("Au.cif")

# Get the force field for gold
calculator = EMT()

# Relax the gold structure (useless since for symmetries it is already relaxed)
relax = CC.calculators.Relax(gold_structure, calculator)
gold_structure_relaxed = relax.static_relax()

# Compute the harmonic phonons
# NOTE: if the code is run with mpirun, the calculation goes in parallel
gold_harmonic_dyn = CC.Phonons.compute_phonons_finite_displacements(gold_structure_relaxed, calculator, supercell = (4,4,4))

# Impose the symmetries and 
# save the dynamical matrix in the quantum espresso format
gold_harmonic_dyn.Symmetrize()
gold_harmonic_dyn.save_qe("harmonic_dyn")


# If the dynamical matrix has imaginary frequencies, remove them
gold_harmonic_dyn.ForcePositiveDefinite()

"""
gold_harmonic_dyn is ready to start the SSCHA calculation.

Now let us initialize the ensemble, and the calculation at 300 K.
We will run a NVT calculation, using 100 configurations at each step
"""

TEMPERATURE = 300
N_CONFIGS = 50
MAX_ITERATIONS = 20

# Initialize the random ionic ensemble
ensemble = sscha.Ensemble.Ensemble(gold_harmonic_dyn, TEMPERATURE)

# Initialize the free energy minimizer
minim = sscha.SchaMinimizer.SSCHA_Minimizer(ensemble)
minim.set_minimization_step(0.01) 

# Initialize the NVT simulation
relax = sscha.Relax.SSCHA(minim, calculator, N_configs = N_CONFIGS,
                          max_pop = MAX_ITERATIONS)

# Define the I/O operations
# To save info about the free energy minimization after each step
ioinfo = sscha.Utilities.IOInfo()
ioinfo.SetupSaving("minim_info")
relax.setup_custom_functions(custom_function_post = ioinfo.CFP_SaveAll)


# Run the NVT simulation (save the stress to compute the pressure)
relax.relax(get_stress = True)

# If instead you want to run a NPT simulation, use
# The target pressure is given in GPa.
#relax.vc_relax(target_press = 0)

# You can also run a mixed simulation (NVT) but with variable lattice parameters
#relax.vc_relax(fix_volume = True)

# Now we can save the final dynamical matrix
# And print in stdout the info about the minimization
relax.minim.finalize()
relax.minim.dyn.save_qe("sscha_T{}_dyn".format(TEMPERATURE))










