    # -*- coding: utf-8 -*-
from __future__ import print_function
from __future__ import division

import sys, os
import numpy as np
import cellconstructor as CC
import cellconstructor.Phonons

import sscha, sscha.Ensemble
import sscha.SchaMinimizer

"""
This test makes a simple relaxation of the sample ensemble
provided within this distribution
"""

def test_update_weights(verbose = False):
    total_path = os.path.dirname(os.path.abspath(__file__))
    os.chdir(total_path)

    rho_final = np.loadtxt("final_rho.dat")
    DATA_PATH = "../../Examples/ensemble_data_test/"

    dyn_start = CC.Phonons.Phonons(os.path.join(DATA_PATH, "dyn"))
    dyn_target = CC.Phonons.Phonons(os.path.join(DATA_PATH, "dyn1_population2"), full_name = True)

    # Perform the minimization
    ens = sscha.Ensemble.Ensemble(dyn_start, 0, dyn_start.GetSupercell())
    ens.load(DATA_PATH, 2, 1000)

    #Update the ensemble
    ens.update_weights_fourier(dyn_target, 0)
    #ens.update_weights(dyn_target, 0)

    delta_rho = np.max(np.abs(ens.rho - rho_final))
    EPS = 1e-7
    
    if verbose:
        print("Maximum difference on rho:", delta_rho)

        if delta_rho > EPS:
            print("Python RHO:")
            print(ens.rho)
            
    assert delta_rho < 2e-5

def test_simple_relax(verbose = False, use_julia = True):
    total_path = os.path.dirname(os.path.abspath(__file__))
    os.chdir(total_path)

    DATA_PATH = "../../Examples/ensemble_data_test/"

    dyn_start = CC.Phonons.Phonons(os.path.join(DATA_PATH, "dyn"))

    # We do not use the dyn1_population2 matrix (the fortran one) because it seems it has problem with the wyckoff minimization.
    # It is printing the structure slightly displaced
    dyn_target = CC.Phonons.Phonons(os.path.join(DATA_PATH, "dyn1_population2"), full_name = True)

    cfp = None
    cfg = None
    if verbose:
        def save_all(minim):
            ka = len(minim.__fe__)

            minim_data = {}
            minim_data["avforce"] = minim.ensemble.get_fourier_forces(False)
            minim_data["rho"] = minim.ensemble.rho
            minim_data["disp_q"] = minim.ensemble.u_disps_qspace.ravel()
            minim_data["disp_r"] = minim.ensemble.u_disps.ravel()
            minim.dyn.save_qe("dyn_relax_ka{}_".format(ka))

            fname = "minim_data_ka{}.npz".format(ka)
            np.savez(fname, **minim_data)
    
        

        grad_structs = []

        def save_gradient(grad_dyn, grad_struct, minim):
            grad_structs.append(grad_struct.ravel())

            np.savetxt("grad_struct.dat", 
                np.concatenate([
                    np.transpose([np.arange(len(grad_structs))]),
                    np.array(grad_structs)
                ], axis=1), header = "step; gradient structure")

        cfp = save_all
        cfg = save_gradient

    # Perform the minimization
    ens = sscha.Ensemble.Ensemble(dyn_start, 0, dyn_start.GetSupercell())
    ens.load(DATA_PATH, 2, 1000)

    minim = sscha.SchaMinimizer.SSCHA_Minimizer(ens)
    minim.use_julia = use_julia
    minim.minim_struct = True
    minim.min_step_dyn = 0.5
    minim.min_step_struc = 0.5
    minim.meaningful_factor = 1e-10
    
    # Go in the working directory to avoid missing up files
    if verbose:
        if minim.use_julia:
            if not os.path.isdir("julia"):
                os.mkdir("julia")
            os.chdir("julia")
        else:
            if not os.path.isdir("standard"):
                os.mkdir("standard")
            os.chdir("standard")

    minim.init()
    minim.run(custom_function_pre=cfp, 
              custom_function_gradient=cfg)
    minim.finalize()

    if verbose:
        lbl = "timer_julia.json"
        if not minim.use_julia:
            lbl = "timer_python.json"
        minim.timer.save_json(lbl)

    # Check the differences in the atomic positions
    delta_s = np.max(np.abs(dyn_target.structure.coords - minim.dyn.structure.coords))
    starting_delta_s = np.max(np.abs(minim.dyn.structure.coords - dyn_start.structure.coords))
    starting_delta = np.max(np.abs(minim.dyn.dynmats[0] - dyn_start.dynmats[0]))
    
    
    # Compare the final dynamical matrices
    delta = np.max(np.abs(dyn_target.dynmats[0] - minim.dyn.dynmats[0]))

    if verbose:
        print("The difference with the original struct is {}".format(starting_delta_s))
        print("The difference on the structure is {}".format(delta_s))
        print("The difference on the original FC is {}".format(starting_delta))
        print("The difference on FC is {}".format(delta))
        
    EPS = 1e-5
    assert delta_s < EPS, "Error, the structure is displaced by {} > {}".format(delta_s, EPS)
    assert delta < EPS, "Error, the difference is: {} > {}".format(delta, EPS)



if __name__ == "__main__":
    test_update_weights(True)
    
    use_julia=True
    if len(sys.argv) > 1:
        use_julia = sys.argv[1] == "julia"
    
    test_simple_relax(True, use_julia)
