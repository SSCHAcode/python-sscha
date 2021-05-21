import pytest
import sys, os
import cellconstructor as CC, cellconstructor.Phonons

import sscha, sscha.Ensemble, sscha.Parallel
from sscha.Parallel import pprint as print
import numpy as np
__SKIP_TEST__ = False
try:
    import fforces as ff
    import fforces.Calculator
    import mpi4py
except:
    __SKIP_TEST__ = True
    
@pytest.mark.skipif(__SKIP_TEST__, reason = "mpi4py and F3Calc required.")
def test_compute_ensemble_parallel(verbose = False):
    total_path = os.path.dirname(os.path.abspath(__file__))
    os.chdir(total_path)
    np.random.seed(0)

    dyn = CC.Phonons.Phonons("dyn_eff", 3)
    ens = sscha.Ensemble.Ensemble(dyn, 250)
    ens.generate(2)

    # Setup a simple harmonic calculator
    calc = ff.Calculator.ToyModelCalculator(dyn, type_cal = "harmx")
    ens.compute_ensemble(calc, compute_stress = False)

    comm = mpi4py.MPI.COMM_WORLD
    if verbose and comm.Get_rank() == 1:# sscha.Parallel.am_i_the_master():
        print("Computed forces:")
        print(ens.forces)
        print("SCHA forces:")
        print(ens.sscha_forces)
        np.savetxt("cf.dat", ens.forces.reshape( (ens.N, 3 * ens.structures[0].N_atoms)))
        np.savetxt("sf.dat", ens.sscha_forces.reshape( (ens.N, 3 * ens.structures[0].N_atoms)))

    comm.barrier()
    
    assert np.max(np.abs( ens.forces - ens.sscha_forces)) < 1e-8, "Rank: {} is not correct.".format(comm.Get_rank())
    
    
if __name__ == "__main__":
    test_compute_ensemble_parallel(verbose = True)
