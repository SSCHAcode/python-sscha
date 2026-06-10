"""Test for an actual FLARE-powered ensemble computation."""
import os
import numpy as np
from ase.io import read

from cellconstructor.Phonons import Phonons, compute_phonons_finite_displacements
from cellconstructor.Structure import Structure
from sscha.Ensemble import Ensemble
from sscha.SchaMinimizer import SSCHA_Minimizer
from sscha.Relax import SSCHA
from sscha.Utilities import IOInfo

from flare.bffs.sgp.calculator import SGP_Calculator


def main():
    """Run with FLARE + SSCHA @ NPT."""
    # =========== GENERAL INPUTS =============== #
    np.random.seed(0)
    structure_filename        = 'Si.pwi'
    number_of_configurations  = 50
    max_iterations            = 10
    temperature_i             = 300
    temperature_f             = 300
    temperature_step          = 10
    pressure                  = 0
    meaningful_factor         = 0.01
    kong_liu_ratio            = 0.5
    minimization_step         = 0.1
    supercell                 = [2,2,2]
    restart_from_previous_dyn = True
    restart_from_ens          = False
    
    atoms = read(structure_filename)

    # =========== FLARE MODEL =============== #
    flare_calc, _ = SGP_Calculator.from_file('./otf_run_flare.json')
    
    # =========== DYNAMICAL MATRIX =============== #
    dyn = Phonons("Si-dynamical-matrix", nqirr=3)
    
    # (*) If you have a pre-existing model, you can compute the dynamical matrix with it
    # structure = Structure()
    # structure.generate_from_ase_atoms(atoms)
    # dyn = compute_phonons_finite_displacements(structure, flare_calc, supercell=supercell)
    
    dyn.Symmetrize()
    dyn.ForcePositiveDefinite()

    if not os.path.exists("./thermal_expansion"):
        os.makedirs("./thermal_expansion")

    # We cycle over several temperatures
    temperature = temperature_i
    volumes = [] 
    temperatures = [] 

    while temperature <= temperature_f:
        ensemble = Ensemble(dyn, temperature)
        
        minim = SSCHA_Minimizer(ensemble)
        minim.set_minimization_step(minimization_step)
        minim.kong_liu_ratio = kong_liu_ratio # default 0.5
        minim.meaningful_factor = meaningful_factor

        relax = SSCHA(
            minimizer=minim,
            ase_calculator=flare_calc,
            N_configs=number_of_configurations,
            max_pop=max_iterations,
            save_ensemble=True,
        )

        ioinfo = IOInfo()
        ioinfo.SetupSaving(f'./thermal_expansion/minim_t{temperature}')
        relax.setup_custom_functions( custom_function_post = ioinfo.CFP_SaveAll)
        
        # Run the NVT simulation
        relax.vc_relax(
            target_press = pressure,
            restart_from_ens = restart_from_ens,
            ensemble_loc = f'./ensembles_P{pressure}_T{temperature}_flare',
        )

        # Print in standard output
        relax.minim.finalize()
        
        # Save the volume and temperature
        volumes.append(relax.minim.dyn.structure.get_volume())
        temperatures.append(temperature)
        relax.minim.dyn.save_qe(f"./thermal_expansion/sscha_T{temperature}_dyn") 
        
        # Start the next simulation from the converged value at this temperature
        if restart_from_previous_dyn:
            dyn = relax.minim.dyn
        
        # Increase temperature
        temperature += temperature_step
        
        # Save thermal expension
        np.savetxt(
            "./thermal_expansion/thermal_expansion.dat",
            np.transpose([temperatures, volumes]),
            header = "Temperature [K]; Volume [A^3]",
        )


if __name__ == '__main__':
    main()

