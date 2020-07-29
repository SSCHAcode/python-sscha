# -*- coding: utf-8 -*-

"""
This example the static phonons are computed using the odd3 correction
introduced by Bianco et al. in "Phys. Rev. B 96, 014111".
We use the ice population as before to test if the odd3 correction is working properly
"""
import numpy as np
import cellconstructor as CC
import cellconstructor.Phonons
import sscha, sscha.Ensemble
from sscha.SchaMinimizer import __RyToCm__ as RyCm
import time

# Load the dynamical matrix
dyn = CC.Phonons.Phonons("../ensemble_data_test/dyn")
dyn.Symmetrize()
w,p = dyn.DyagDinQ(0)

# Load the ensemble
ens = sscha.Ensemble.Ensemble(dyn, 0)
ens.load("../ensemble_data_test", 2, 1000)

# Unwrap the ensemble
ens.unwrap_symmetries()

# Compute the odd correction
t1 = time.time()
new_dyn_supercell = ens.get_odd_correction()
t2 = time.time()

#other_new_dyn = ens.get_odd_correction(store_v3 = False, progress = True)
other_new_dyn = ens.get_odd_correction(include_v4 = True, progress = True, v4_conv_thr = 1e-8)
t3 = time.time()

# Copy the matrix into the Phonons class
dyn.dynmats[0] = new_dyn_supercell
dyn2 = dyn.Copy()
dyn2.dynmats[0] = other_new_dyn

# Perform the symmetrization
dyn.Symmetrize()
dyn2.Symmetrize()

# Save the dynamical matrix with the odd3 corrections
dyn.save_qe("dyn_plus_odd_new")

odd = dyn.Copy()
odd.dynmats[0] = np.load("odd_corr.npy")
odd.save_qe("odd_new_nosym")
odd.Symmetrize()
odd.save_qe("odd_new_sym")


print "Elapsed time to compute odd3:", t2 - t1, "s"
print "Elapsed time to compute odd3 without storing it:", t3- t2, "s"
print "Difference between the result of storing and not storing it:"
print np.sqrt(np.sum( (new_dyn_supercell - other_new_dyn)**2))
print ""
print "Results saved into 'dyn_plus_odd_new'"

print
# Get the phonon frequencies and print them
w1,p = dyn.DyagDinQ(0)
w2,p = dyn2.DyagDinQ(0)
print "    %10s  |%10s  |%10s" % ("SCHA", "ODD3 v1", "ODD3 v2")
print "\n".join(["%2d) %10.4f  |%10.4f  |%10.4f  cm-1" % (i+1, w[i] * RyCm, w1[i] * RyCm, w2[i] * RyCm) for i in range(len(w))])
