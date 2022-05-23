# Import the sscha code
import sscha, sscha.Ensemble, sscha.SchaMinimizer, sscha.Relax, sscha.Utilities

# Import the cellconstructor library to manage phonons
import cellconstructor as CC, cellconstructor.Phonons
import cellconstructor.Structure, cellconstructor.calculators

# Import the DFT calculator
import cellconstructor.calculators

# Import numerical and general pourpouse libraries
import numpy as np, matplotlib.pyplot as plt
import sys, os


# Initialize the DFT (Quantum Espresso) calculator for gold
# The input data is a dictionary that encodes the pw.x input file namelist
input_data = {
    'control' : {
        # Avoid writing wavefunctions on the disk
        'disk_io' : 'None',
        # Where to find the pseudopotential
        'pseudo_dir' : '.'
    },
    'system' : {
        # Specify the basis set cutoffs
        'ecutwfc' : 45,   # Cutoff for wavefunction
        'ecutrho' : 45*4, # Cutoff for the density
        # Information about smearing (it is a metal)
        'occupations' : 'smearing',
        'smearing' : 'mv',
        'degauss' : 0.03
    },
    'electrons' : {
        'conv_thr' : 1e-8
    }
}

# the pseudopotential for each chemical element
# In this case just Gold
pseudopotentials = {'Au' : 'Au_ONCV_PBE-1.0.oncvpsp.upf'}

# the kpoints mesh and the offset
kpts = (1,1,1)
koffset = (1,1,1)

command = 'mpirun -np 4 pw.x -npool 1 -i PREFIX.pwi > PREFIX.pwo'

# Prepare the quantum espresso calculator
calculator = CC.calculators.Espresso(input_data,
                                     pseudopotentials,
                                     command = command,
                                     kpts = kpts,
                                     koffset = koffset)



TEMPERATURE = 300
N_CONFIGS = 50
MAX_ITERATIONS = 20
START_DYN = 'start_dyn'
NQIRR = 13

# Let us load the starting dynamical matrix
gold_dyn = CC.Phonons.Phonons(START_DYN, NQIRR)

# Initialize the random ionic ensemble
ensemble = sscha.Ensemble.Ensemble(gold_dyn, TEMPERATURE)

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










