import cellconstructor as CC, cellconstructor.Phonons
import sscha, sscha.Ensemble
import sscha.SchaMinimizer, sscha.Relax
import sscha.Utilities

import ase, ase.calculators, ase.calculators.emt
import numpy as np
import os, sys

import pytest


@pytest.mark.release
def test_gold(use_julia=True, verbose=False):
    supercell=(4,4,4)
    
    # Change directory to the current script
    total_path = os.path.dirname(os.path.abspath(__file__))
    os.chdir(total_path)
    
    cif_file="Au_mp-81_primitive.cif"
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
            N_configs = 500,
            max_pop = 50)
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

    # Load the final dynamical matrix and compare
    final_dyn = CC.Phonons.Phonons("dyn_good_gold_sscha_300_", 13)


    gold_dyn = relax.minim.dyn
    for i in range(len(final_dyn.q_tot)):
        delta = final_dyn.dynmats[i] - gold_dyn.dynmats[i]
        if verbose:
            print("Delta q = ", final_dyn.q_tot[i])
            print("Delta max = ", np.max(np.abs(delta)))

        assert delta < 0.005, "Error in the dynamical matrix at q = " + str(final_dyn.q_tot[i])

    delta = final_dyn.structure.unit_cell[0,0] - gold_dyn.structure.unit_cell[0,0]
    assert np.abs(delta) < 0.0005, "Error in the lattice parameter"

    if verbose:
        print("Delta cell = ", np.abs(delta))


if __name__ == "__main__":
    use_julia=False
    if len(sys.argv) > 1:
        use_julia = sys.argv[1] == "julia"
                
    test_gold(verbose=True, use_julia=use_julia)
