import cellconstructor as CC, cellconstructor.Phonons
import sscha, sscha.Ensemble
import sscha.SchaMinimizer, sscha.Relax
import sscha.Utilities

import ase, ase.calculators, ase.calculators.emt
import numpy as np
import os, sys

import pytest


def test_gold(use_julia=True, verbose=False):
    supercell=(2,2,2)
    
    # Change directory to the current script
    total_path = os.path.dirname(os.path.abspath(__file__))
    os.chdir(total_path)
    
    cif_file="../test_release/Au_mp-81_primitive.cif"
    temperature = 300

    # Load the structure
    struct = CC.Structure.Structure()
    struct.read_generic_file(cif_file)
    gold_dyn = CC.Phonons.compute_phonons_finite_displacements(struct, 
        ase.calculators.emt.EMT(), supercell=supercell)

    # Setup a sscha
    gold_dyn.Symmetrize()
    gold_dyn.ForcePositiveDefinite()

    ensemble = sscha.Ensemble.Ensemble(gold_dyn, temperature)
    minim = sscha.SchaMinimizer.SSCHA_Minimizer(ensemble)
    minim.set_minimization_step(0.001)
    minim.meaningful_factor = 1e-5
    minim.use_julia = use_julia

    relax = sscha.Relax.SSCHA(minim,
            ase_calculator=ase.calculators.emt.EMT(),
            N_configs = 50,
            max_pop = 3)
    relax.save_ensemble=False

    cfp=None
    if verbose:
        ioinfo = sscha.Utilities.IOInfo()
        ioinfo.SetupSaving("minim")
        cfp = ioinfo.CFP_SaveAll

    relax.setup_custom_functions(
            custom_function_post=cfp)

    relax.vc_relax(target_press=0)
    
    if verbose:
        relax.minim.dyn.save_qe("final_dyn")

if __name__ == "__main__":
    use_julia=True
    if len(sys.argv) > 1:
        if sys.argv[1].lower() != "julia":
            use_julia=False
                
    test_gold(verbose=True, use_julia=use_julia)
