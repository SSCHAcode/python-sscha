# -*- coding: utf-8 -*-
from __future__ import print_function
from __future__ import division

import os
import tempfile
import threading

import cellconstructor as CC
import cellconstructor.Phonons

import sscha
import sscha.Cluster
import sscha.Ensemble
import sscha.Relax
import sscha.SchaMinimizer
import sscha.Utilities

"""
Regression test for issue #114: save_binary failed with
TypeError: cannot pickle '_thread.lock' object
when the relax object holds a cluster that already ran a calculation
(Cluster.compute_ensemble_batch stores a threading.Lock on the cluster).
"""


def test_save_binary_relax_with_cluster(verbose=False):
    total_path = os.path.dirname(os.path.abspath(__file__))
    os.chdir(total_path)

    DATA_PATH = "../../Examples/ensemble_data_test/"

    dyn = CC.Phonons.Phonons(os.path.join(DATA_PATH, "dyn"))

    ens = sscha.Ensemble.Ensemble(dyn, 0, dyn.GetSupercell())
    ens.load(DATA_PATH, 2, 10)

    minim = sscha.SchaMinimizer.SSCHA_Minimizer(ens)

    cluster = sscha.Cluster.Cluster(hostname="localhost")
    relax = sscha.Relax.SSCHA(minim, N_configs=10, max_pop=2,
                              cluster=cluster)

    # Cluster.compute_ensemble_batch leaves a threading.Lock on the
    # cluster after the ensemble calculation; reproduce that state.
    relax.cluster.lock = threading.Lock()

    with tempfile.TemporaryDirectory() as tmpdir:
        filename = os.path.join(tmpdir, "relax.bin")
        sscha.Utilities.save_binary(relax, filename)

        loaded = sscha.Utilities.load_binary(filename)

    # The lock is transient runtime state and must come back unset.
    assert loaded.cluster.lock is None
    assert loaded.cluster.hostname == "localhost"
    assert loaded.N_configs == relax.N_configs
    assert loaded.minim.ensemble.N == ens.N

    if verbose:
        print("save_binary/load_binary round trip succeeded")


if __name__ == "__main__":
    test_save_binary_relax_with_cluster(True)
