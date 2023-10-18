"""Example of an actual AiiDA-powered SSCHA run."""
import os
import numpy as np

from cellconstructor.Phonons import Phonons
from sscha.aiida_ensemble import AiiDAEnsemble
from sscha.SchaMinimizer import SSCHA_Minimizer
from sscha.Relax import SSCHA
from sscha.Utilities import IOInfo

from aiida import load_profile

from aiida_quantumespresso.common.types import ElectronicType

load_profile()

# =========== DEFINITION OF INPUTS =============== #
np.random.seed(0)
# Define the temperature range (in K)
t_start = 300
t_end = 800
delta_t = 50
target_pressure = 0
number_q_irreducible = 4
kpoints_distance = 0.4 # kpoints mesh
number_of_configurations = 100
max_iterations = 15
data_directory = './minimization'
pw_code_label = 'pw_7.2@hlrn-gottingen'
dyn_filepath = './harmonic_dyn/dynamical-matrix-'
directory_output = "./thermal_expansion" 

aiida_inputs = dict(
    pw_code_label=pw_code_label,
    protocol='precise',
    overrides={
        'pseudo_family': 'SSSP/1.3/PBEsol/precision',
        'kpoints_distance': 0.4,
        'pw':{
            'parameters':{
                'SYSTEM':{
                    # ...
                },
                'ELECTRONS':{
                    'mixing_beta': 0.7
                },
            },
            'settings':{
                'cmdline': ['-nk', '8']
            },
        }
    },
    options={
      'resources':{
          'num_machines': 1,
          'num_mpiprocs_per_machine': 96,
          'num_cores_per_mpiproc': 1
      },
    #   'account': 'myproject',
      'max_wallclock_seconds': int(1*60*60),
    },
    electronic_type=ElectronicType.INSULATOR,
    # group_label='Group/Where/To/Add/PwNodes', # You need to create it first
)
# =========== DEFINITION OF INPUTS =============== #

dyn = Phonons(dyn_filepath, number_q_irreducible)
dyn.ForcePositiveDefinite()
dyn.Symmetrize()

if not os.path.exists(directory_output):
    os.makedirs("./thermal_expansion")

# We cycle over several temperatures
t = t_start
volumes = [] 
temperatures = [] 

while t <= t_end:
    # Change the temperature
    ensemble = AiiDAEnsemble(dyn, t)
    minim = SSCHA_Minimizer(ensemble)
    minim.set_minimization_step(0.1)
    relax = SSCHA(
        minimizer=minim,
        aiida_inputs=aiida_inputs,
        N_configs=number_of_configurations,
        max_pop=max_iterations,
        save_ensemble=True,
    )
    # Setup the I/O
    
    ioinfo = IOInfo()
    ioinfo.SetupSaving( os.path.join(directory_output, "minim_t{}".format(t))) 
    relax.setup_custom_functions( custom_function_post = ioinfo.CFP_SaveAll)
    # Run the NPT simulation
    relax.vc_relax(
        target_press = target_pressure, 
        restart_from_ens = False,
        ensemble_loc = f'ensembles_T{t}'
    )

    # Save the volume and temperature
    volumes.append(relax.minim.dyn.structure.get_volume())
    temperatures.append(t)
    # Start the next simulation from the converged value at this temperature
    relax.minim.dyn.save_qe(os.path.join(directory_output, f"sscha_T{t}_dyn"))
    dyn = relax.minim.dyn
    # Print in standard output
    relax.minim.finalize()
    # Update the temperature
    t += delta_t
    # Save the thermal expansion
    np.savetxt(os.path.join(directory_output, "thermal_expansion.dat"),
                np.transpose([temperatures, volumes]),
                header = "Temperature [K]; Volume [A^3]")