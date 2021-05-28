from __future__ import print_function
from __future__ import division

import sys, os

import numpy as np
import cellconstructor as CC
import cellconstructor.Phonons

import sscha, sscha.Ensemble, sscha.Utilities
import sscha.SchaMinimizer


MODE_LOCKED = 3
def test_lock_manyq():
    total_path = os.path.dirname(os.path.abspath(__file__))
    os.chdir(total_path)

    # Here we test the lock modes on many q points

    dyn = CC.Phonons.Phonons("data_graphite/dyn_gen_pop1_", 13)
    ens = sscha.Ensemble.Ensemble(dyn, 0, dyn.GetSupercell())
    ens.load_bin("data_graphite", 1)

    # Setup the saving
    IO_freq = sscha.Utilities.IOInfo()
    IO_freq.SetupSaving("frequencies_graphite.dat") # The frequencies will be saved at each step of the minimization

    minim = sscha.SchaMinimizer.SSCHA_Minimizer(ens)
    minim.min_step_dyn = 4e-2
    minim.min_step_struc = 4e-2
    minim.meaningful_factor = 1e-6 
    minim.max_ka = 3

    minim.init()
    minim.print_info()

    # save the frequencies for each q
    ws = []
    for i in range(len(dyn.q_tot)):
        w, p = dyn.DyagDinQ(i)
        ws.append(w)

    # Setup mode locking (the lowest optical mode)
    mode_cons = sscha.Utilities.ModeProjection(dyn)
    mode_cons.SetupFreeModes(MODE_LOCKED, MODE_LOCKED + 1)

    minim.run(custom_function_gradient = mode_cons.CFG_ProjectOnModes,
              custom_function_post = IO_freq.CFP_SaveFrequencies)

    minim.finalize()
    minim.dyn.save_qe("data_graphite/final_dyn")

    # Check if the third mode frequency is the only one that changed for each q point
    for i in range(len(dyn.q_tot)):
        new_w, p = minim.dyn.DyagDinQ(i)

        delta_w = np.abs(new_w - ws[i]) * CC.Units.RY_TO_CM

        for j in range(len(new_w)):
            if j == MODE_LOCKED:
                continue
            
            assert delta_w[j] < 1e-4, "Error, mode {} of q = {} changed (locked = {})".format(j, i, MODE_LOCKED)
            


if __name__ == "__main__":
    test_lock_manyq()
