import cellconstructor as CC, cellconstructor.Phonons
import sscha, sscha.Ensemble
import sscha.SchaMinimizer
import sys, os


def test_load_noncomputed_ensemble():
    # Move the current directory to the test directory
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    # Load the dynamical matrix
    dyn = CC.Phonons.Phonons("dyn", 4)

    # Load the ensemble
    ens = sscha.Ensemble.Ensemble(dyn, 300)

    # Load the ensemble from the output of the calculator
    # In this case, the pwo files are output of the quantum espresso program.
    # Any output file that ASE is able to read can be used to load the ensemble.
    ens.load_from_calculator_output(directory="data", out_ext=".pwo")


    # Run the minimization
    minim = sscha.SchaMinimizer.SSCHA_Minimizer(ens)
    minim.init()
    minim.set_minimization_step(0.01)
    minim.run()
    minim.finalize()


if __name__ == "__main__":
    test_load_noncomputed_ensemble()
