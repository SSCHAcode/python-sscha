import cellconstructor as CC
import cellconstructor.Phonons

import sscha, sscha.Ensemble, sscha.SchaMinimizer

import numpy as np
import sys, os


def test_root_step_identity():
    total_path = os.path.dirname(os.path.abspath(__file__))
    os.chdir(total_path)
    
    dyn = CC.Phonons.Phonons("dyn_mono_10x10x1_full", 14)
    dyn.Symmetrize()

    # Prepare the virtual root step
    nat = dyn.structure.N_atoms
    
    dyn_q  = np.zeros((len(dyn.q_tot), 3*nat, 3*nat), dtype = np.complex128)
    grad_q = np.zeros((len(dyn.q_tot), 3*nat, 3*nat), dtype = np.complex128)

    # Complete the dynq
    for iq in range(len(dyn.q_tot)):
        dyn_q[iq, :, :] = dyn.dynmats[iq]

    new_dyn = sscha.SchaMinimizer.PerformRootStep(dyn_q, grad_q, step_size = 0,
                                                  root_representation = "root2")

    disp = np.max(np.abs(new_dyn - dyn_q))
    assert disp < 1e-12, "Error on the root step, margin = {}".format(disp)


if __name__ == "__main__":
    test_root_step_identity()

    
