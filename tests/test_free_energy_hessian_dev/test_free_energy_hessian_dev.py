import numpy as np
import pytest
import cellconstructor as CC
import cellconstructor.Phonons

import sscha.Ensemble

@pytest.mark.release
def test_free_energy_hessian_dev():
    N_pop = 1
    Temperature = 0
    dyn_sscha_final = CC.Phonons.Phonons("./sscha/dyn_end_population1_", nqirr=3)
    dyn_sscha_final.Symmetrize()

    ensemble = sscha.Ensemble.Ensemble(dyn_sscha_final, T0=Temperature, supercell = dyn_sscha_final.GetSupercell())
    ensemble.load_bin("./sscha/data_ensemble", population_id = N_pop)
    ensemble.update_weights(dyn_sscha_final, Temperature)
    ensemble.current_dyn.Symmetrize() # The dev method requires symmetrized polvecs.

    dyn_hessian = CC.Phonons.Phonons("hessianv4_", nqirr=3)
    dyn_hessian_dev = ensemble.get_free_energy_hessian_dev(include_v4 = True)

    dyn_hessian_dev.Symmetrize()

    np.testing.assert_allclose(
        dyn_hessian_dev.dynmats, 
        dyn_hessian.dynmats, 
        atol=5e-5,
        err_msg="Hessian matrix from dev implementation does not match reference!"
    )

    dyn_hessian.save_qe("hessianv4_new_")

if __name__ == "__main__":
    test_free_energy_hessian_dev()


