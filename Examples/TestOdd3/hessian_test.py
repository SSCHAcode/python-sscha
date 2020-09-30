import cellconstructor as CC
import cellconstructor.Phonons

import sscha, sscha.Ensemble


# Load the dynamical matrix (sscha)
# nqirr = 1, fildyn = "../ensemble_data_test/dyn"
dyn_sscha = CC.Phonons.Phonons("../ensemble_data_test/dyn", 1)

# Load the ensemble (T = 100 K)
ensemble = sscha.Ensemble.Ensemble(dyn_sscha, 100, supercell = (1,1,1))
# data_dir = "../ensemble_data_test"
# population = 2
# n_random = 1000
ensemble.load("../ensemble_data_test", 2, 1000)


# Now compute the free energy hessian
# Using the raffaello code
# You can include also fourth order by setting include_v4 to true.
hessian_dyn = ensemble.get_free_energy_hessian(include_v4 = False)

# Save the free energy hessian in qe format
hessian_dyn.save_qe("free_energy_hessian")
