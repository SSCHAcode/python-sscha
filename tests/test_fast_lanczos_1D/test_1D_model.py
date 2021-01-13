from __future__ import print_function

import FF1D

import ase, ase.atoms
import warnings

import numpy as np

import cellconstructor as CC
import cellconstructor.Phonons
import sscha, sscha.Ensemble, sscha.SchaMinimizer
import sscha.Relax, sscha.DynamicalLanczos

import sys, os

def test_sscha_1D():
    total_path = os.path.dirname(os.path.abspath(__file__))
    os.chdir(total_path)

    # Get a standard dynamical matrix
    atoms = ase.atoms.Atoms("H2", [(0.1,0.1,0.1), (0.2,0.12,0.13)], )

    struct = CC.Structure.Structure()
    struct.generate_from_ase_atoms(atoms)
    struct.unit_cell = np.eye(3) * 10
    struct.unit_cell[1,1] = 9
    struct.unit_cell[2,2] = 8
    struct.has_unit_cell = True
    

    NEW_MINIM = True
    if os.path.exists("dyn1"):
        NEW_MINIM = False
        dyn = CC.Phonons.Phonons("dyn")
    else:
        dyn = CC.Phonons.Phonons(struct)

        # Randomize the starting dynamical matrix
        dyn.dynmats[0][:,:] = np.random.uniform(0, 0.2, size = dyn.dynmats[0].shape) 
        dyn.dynmats[0] += dyn.dynmats[0].T

        dyn.Symmetrize()
        dyn.ForcePositiveDefinite()

    # Setup the calculator
    calc = FF1D.ToyModel1D()
    

    # Start the SSCHA
    ensemble = sscha.Ensemble.Ensemble(dyn, 0)
    minim = sscha.SchaMinimizer.SSCHA_Minimizer(ensemble)
    minim.min_step_dyn = 1e-4
    minim.minim_struct = False
    minim.meaningful_factor = 0.01
    #minim.max_ka = 50

    minim.print_info()

    relax = sscha.Relax.SSCHA(minim, ase_calculator = calc,
                              N_configs = 100, max_pop = 1)

    if NEW_MINIM:
    
        relax.relax(start_pop = 1)

        relax.minim.min_step_dyn = 5e-3
        relax.max_pop = 2
        relax.relax(start_pop = 1)

        relax.minim.min_step_dyn = 1e-2
        relax.max_pop = 2
        relax.relax(start_pop = 1)

        relax.N_configs = 1000
        relax.minim.minim_struct = True
        relax.minim.min_step_dyn = 0.05
        relax.minim.min_step_struc = 0.05
        relax.max_pop = 5
        relax.relax(start_pop = 1)

        relax.minim.dyn.save_qe("dyn")


    relax.N_configs = 10000
    relax.minim.minim_struct = True
    relax.minim.min_step_dyn = 0.2
    relax.minim.min_step_struc = 0.2
    relax.max_pop = 20
    relax.relax(start_pop = 1)

    relax.minim.dyn.save_qe("dyn_final")


def solve_finite_temperature(T = 1000):

    total_path = os.path.dirname(os.path.abspath(__file__))
    os.chdir(total_path)

    final_dynname = "dyn_final_T_{:d}_".format(T)

    # Load the dynamical matrix
    dyn = CC.Phonons.Phonons("dyn_final")
    if os.path.exists(final_dynname + "1"):
        dyn = CC.Phonons.Phonons(final_dynname)

    
    # Setup the calculator
    calc = FF1D.ToyModel1D()

    # Start the SSCHA
    ensemble = sscha.Ensemble.Ensemble(dyn, T)
    minim = sscha.SchaMinimizer.SSCHA_Minimizer(ensemble)
    minim.min_step_dyn = 5e-1
    minim.minim_struct = True
    minim.min_step_struc = 5e-1
    minim.meaningful_factor = 1e-6
    #minim.max_ka = 50

    minim.print_info()


    # # First relaxation
    relax = sscha.Relax.SSCHA(minim, ase_calculator = calc,
                              N_configs = 1000, max_pop = 4)

    # relax.relax(start_pop = 1)

    # Improve the number of configurations
    relax.N_configs = 50000
    relax.relax(start_pop = 1)

    # Save the dynamical matrix
    relax.minim.dyn.save_qe(final_dynname)


def test_static_odd_ft(T = 20):
    dynname = "dyn_final_T_{:d}_".format(T)

    if not os.path.exists(dynname + "1"):
        solve_finite_temperature(T)

    dyn = CC.Phonons.Phonons(dynname)
    ens = sscha.Ensemble.Ensemble(dyn, T)

    # Generate the random configurations
    N_RANDOM = 10000
    ens.generate(N_RANDOM)

    # Setup the calculator
    calc = FF1D.ToyModel1D()

    # Compute forces and energies for the ensemble
    ens.compute_ensemble(calc, compute_stress = False)

    # Get odd correction with v4
    hessian = ens.get_free_energy_hessian(include_v4 = True, use_symmetries = True)
    hessian.save_qe("hessian_v4_T_{}_".format(T))

    # Setup the Lanczos (selecting only the interesting mode)
    lanc = sscha.DynamicalLanczos.Lanczos(ens, unwrap_symmetries = False, select_modes = np.array([True, True, True, False, False, True]))
    lanc.init()
    print("N_mod: {}".format(lanc.n_modes))

    lanc.ignore_v4 = False
    lanc.ignore_v3 = False
    #lanc.get_full_L_operator_FT(symmetrize = False)

    hessian_2 = lanc.run_biconjugate_gradient(save_g = "G_lanc.npy", tol = 1e-6, maxiter = 500, use_preconditioning = True, algorithm = "bicgstab")
    np.savetxt("hessian_lanc_T_{}.dat".format(T), hessian_2)

    # Get the frequencies
    ws2, pols = np.linalg.eigh(hessian_2)

    ws = np.sign(ws2) * np.sqrt(np.abs(ws2)) * CC.Units.RY_TO_CM

    print("Frequencies from BICG:")
    print("\n".join(["{:.4f} cm-1".format(x) for x in ws]))


def test_transposition(T = 500):
    total_path = os.path.dirname(os.path.abspath(__file__))
    os.chdir(total_path)
    
    # Prepare the SSCHA
    dynname = "dyn_final_T_{:d}_".format(T)

    if not os.path.exists(dynname + "1"):
        solve_finite_temperature(T)

    dyn = CC.Phonons.Phonons(dynname)
    ens = sscha.Ensemble.Ensemble(dyn, T)

    # Generate the random configurations
    N_RANDOM = 10000
    ens.generate(N_RANDOM)

    # Setup the calculator
    calc = FF1D.ToyModel1D()

    # Compute forces and energies for the ensemble
    ens.compute_ensemble(calc, compute_stress = False)

    # Setup the Lanczos
    lanc = sscha.DynamicalLanczos.Lanczos(ens, unwrap_symmetries = False, select_modes = np.array([True, True, True, False, True, True]))
    lanc.init()
    print("N_modes: {}".format(lanc.n_modes))

    lanc.ignore_v4 = True
    lanc.ignore_v3 = False
    lanc.ignore_harmonic = False
    #lanc.get_full_L_operator_FT(symmetrize = False)

    # Get the full L matrix and its transposed
    L = sscha.DynamicalLanczos.get_full_L_matrix(lanc, False)
    Lt = sscha.DynamicalLanczos.get_full_L_matrix(lanc, True)

    disp = np.max(np.abs(L - Lt.T)) / np.max(np.abs(L))

    print("Distance: {}".format(disp))

    np.savetxt("L.dat", L)
    np.savetxt("Lt.dat", Lt)

    assert disp < 1e-6

    
def test_harmonic_inversion(T = 500):
    """
    This subroutine test the preconditioner M matrix at finite temperature
    """
    total_path = os.path.dirname(os.path.abspath(__file__))
    os.chdir(total_path)
    
    # Prepare the SSCHA
    dynname = "dyn_final_T_{:d}_".format(T)

    if not os.path.exists(dynname + "1"):
        solve_finite_temperature(T)

    dyn = CC.Phonons.Phonons(dynname)
    ens = sscha.Ensemble.Ensemble(dyn, T)

    # Generate the random configurations
    N_RANDOM = 10000
    ens.generate(N_RANDOM)

    # Setup the calculator
    calc = FF1D.ToyModel1D()

    # Compute forces and energies for the ensemble
    ens.compute_ensemble(calc, compute_stress = False)

    # Setup the Lanczos
    lanc = sscha.DynamicalLanczos.Lanczos(ens, unwrap_symmetries = False, select_modes = np.array([True, True, True, False, False, True]))
    lanc.init()
    print("N_modes: {}".format(lanc.n_modes))

    lanc.ignore_v4 = True
    lanc.ignore_v3 = False
    lanc.ignore_harmonic = False
    #lanc.get_full_L_operator_FT(symmetrize = False)

    # Test if the preconditioning  diverges with one mode (it should not)
    n_tot = len(lanc.psi)
    M_direct = np.zeros((n_tot, n_tot), dtype= np.double)

    v = np.zeros(lanc.psi.shape, dtype = np.double)
    for i in range(n_tot):
        v[:] = 0
        v[i] = 1
        M_direct[:, i] = lanc.M_linop.matvec(v)

    np.savetxt("M.dat", M_direct, header = "The inverse of the L matrix")
        
    
    
if __name__ == "__main__":
    #test_harmonic_inversion()
    #test_sscha_1D()
    #test_preconditioning()
    test_transposition()

    exit()
    T_START = 0
    T_END = 10000
    DT = 500

    t = T_START
    while t <= T_END:
        test_static_odd_ft(T = t)
        t += DT
