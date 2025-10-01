import numpy as np
import cellconstructor as CC, cellconstructor.Phonons
import sscha, sscha.Ensemble, sscha.SchaMinimizer
import sscha.Relax
import ase, ase.calculators, ase.calculators.emt


def test_1d_asr_relax(verbose=False):
    # Set the current working directory
    total_path = os.path.dirname(os.path.abspath(__file__))
    os.chdir(total_path)


    # Load the dynamical matrix
    dyn = CC.Phonons.Phonons("1ddyn_asr")

    # Impose it is a 1D system
    dyn.structure.one_dim_axis = 2

    # Apply a small random perturbation to the dynamical matrix
    dyn.structure.coords += np.random.normal(0, 0.01, dyn.structure.coords.shape)

    # Impose the symmetries and the ASR of a 1D system
    dyn.Symmetrize()

    # Check that we have 4 acoustic modes
    w, p = dyn.DiagonalizeSupercell()
    acoustic = dyn.structure.get_asr_modes(p)
    assert np.sum(acoustic.astype(int)) == 4, "Starting dynamical matrix does not have 4 acoustic modes"

    if verbose:
        print("Starting acoustic modes")
        print(acoustic)

    ens = sscha.Ensemble.Ensemble(dyn, 100)
    minim = sscha.SchaMinimizer.SSCHA_Minimizer(ens)
    
    calc = ase.calculators.emt.EMT()

    relax = sscha.Relax.SSCHA(minim, 
                              ase_calculator=calc,
                              N_configs=10,
                              max_pop=2)
    relax.relax()


    # Check if we have
    w, p = relax.minim.dyn.DiagonalizeSupercell()
    acoustic = relax.minim.dyn.structure.get_asr_modes(p)
    assert np.sum(acoustic.astype(int)) == 4, "Final dynamical matrix does not have 4 acoustic modes"

    if verbose:
        print("Final acoustic modes")
        print(acoustic)


if __name__ == "__main__":
    test_1d_asr_relax(verbose=True)
