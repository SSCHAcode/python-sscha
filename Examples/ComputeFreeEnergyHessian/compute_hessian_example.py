from __future__ import print_function
from __future__ import division

# Import the modules to read the dynamical matrix
import cellconstructor as CC
import cellconstructor.Phonons

# Import the SCHA modules
import sscha, sscha.Ensemble


# Here the input information
DATA_DIR = "../ensemble_data_test"
N_RANDOM = 10000
DYN_PREFIX = "../ensemble_data_test/dyn"
FINAL_DYN = "../ensemble_data_test/dyn"
SAVE_PREFIX = "dyn_plus_odd"
NQIRR = 1
Tg = 0
T = 0
POPULATION = 2
INCLUDE_V4 = False

INFO = """
In this example we compute the free energy hessian for the ice XI.

The ensemble has been generated with the dynamical matrix at:
{}

And to compute the hessian we will use reweighting at:
{}

The original temperature was {} K, we use reweighting to {} K.
The ensemble is located at: {} with id = {}
The number of configuration is {}.
Do we include the v4 in the calculation? {}

The result will be saved in: {}

""".format(DYN_PREFIX, FINAL_DYN, Tg, T, DATA_DIR, POPULATION,
           N_RANDOM, INCLUDE_V4, SAVE_PREFIX)


print(INFO)
print()
print(" ======== RUNNING ======== ")
print()

print("Loading the original dynamical matrix...")
dyn = CC.Phonons.Phonons(DYN_PREFIX, NQIRR)
print("Symmetrizing...")
dyn.Symmetrize()

print("Loading the current dynamical matrix...")
final_dyn = CC.Phonons.Phonons(FINAL_DYN, NQIRR)

print("Loading the ensemble...")
ens = sscha.Ensemble.Ensemble(dyn, Tg, dyn.GetSupercell())
ens.load(DATA_DIR, POPULATION, N_RANDOM)
# If the ensemble was saved in binary format, load it with
#ens.load_bin(DATA_DIR, POPULATION)

print("Updating the importance sampling...")
ens.update_weights(final_dyn, T)

print("")

print("Computing the free energy hessian...")
print("(This may take a while)")
# Set get_full_hessian to false to have only the odd correction
# Usefull if you want to study the convergence with the number of configuration
dyn_hessian = ens.get_free_energy_hessian(include_v4 = INCLUDE_V4,
                                          get_full_hessian = True,
                                          verbose = True)

print("Saving the hessian to {}...".format(SAVE_PREFIX))
dyn_hessian.save_qe(SAVE_PREFIX)
print("Done.")

