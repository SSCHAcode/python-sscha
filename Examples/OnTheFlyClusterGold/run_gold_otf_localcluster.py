"""On-the-fly FLARE learning on gold (EMT) through LocalCluster.

Recipe B from the on-the-fly FLARE cluster plan: exercises the *exact*
remote-cluster code path (input files, tar, submission script, result
retrieval) mocked locally with bash + shutil (no ssh/scp, no SLURM).

    python run_gold_otf_localcluster.py

Requires: sscha (editable), flare, ase. No pw.x, no AiiDA, no SLURM.
"""
import os
import numpy as np

import cellconstructor as CC, cellconstructor.Phonons
import sscha, sscha.Ensemble, sscha.SchaMinimizer, sscha.Relax
import sscha.LocalCluster
import sscha.ClusterCalculators as cluster_calcs

from ase.calculators.emt import EMT


# Keep kernel alive (flare bug, see tests/test_otf/test_otf_base.py).
_KERNEL_REFS = []


def get_sgp_calc_au():
    """Return an empty FLARE SGP calculator for gold (single species)."""
    from flare.bffs.sgp._C_flare import NormalizedDotProduct, B2
    from flare.bffs.sgp import SGP_Wrapper
    from flare.bffs.sgp.calculator import SGP_Calculator

    cutoff = 4.0
    kernel = NormalizedDotProduct(sigma=2.0, power=2)
    b2 = B2("chebyshev", "quadratic", [0.0, cutoff], [],
            [1, 4, 3], cutoff * np.ones((1, 1)))
    sgp = SGP_Wrapper(
        [kernel], [b2], cutoff,
        sigma_e=0.1, sigma_f=0.1, sigma_s=0.1,
        species_map={79: 0},
        single_atom_energies={0: 0.0},
        variance_type="local",
        energy_training=True, force_training=True, stress_training=True,
        opt_method="L-BFGS-B", max_iterations=5,
    )
    _KERNEL_REFS.append(kernel)
    return SGP_Calculator(sgp)


def get_gold_dyn(supercell=(2, 2, 2)):
    """Harmonic Au (EMT), 8 atoms with the default supercell."""
    from ase.build import bulk
    struct = CC.Structure.Structure()
    struct.generate_from_ase_atoms(bulk("Au", "fcc", a=4.0782, cubic=False))
    dyn = CC.Phonons.compute_phonons_finite_displacements(struct, EMT(), supercell=supercell)
    dyn.Symmetrize()
    dyn.ForcePositiveDefinite()
    return dyn


def get_otf_ensemble(dyn, temperature=300, output_name="otf_gold"):
    np.random.seed(0)
    ensemble = sscha.Ensemble.Ensemble(dyn, temperature)
    ensemble.set_otf(
        get_sgp_calc_au(),
        std_tolerance_factor=100,
        max_atoms_added=-1,
        update_style="add_n",
        update_threshold=None,
        train_hyps=(1, np.inf),
        output_name=output_name,
    )
    return ensemble


def main_localcluster():
    dyn = get_gold_dyn()
    ensemble = get_otf_ensemble(dyn)

    calc = cluster_calcs.ASEFileCalculator(EMT())

    cluster = sscha.LocalCluster.LocalCluster("localhost")
    cluster.workdir = os.path.abspath("otf_cluster/remote")
    cluster.local_workdir = os.path.abspath("otf_cluster/local") + "/"
    cluster.submit_command = "bash"       # run the script directly (blocking)
    cluster.nonblocking_command = False   # do NOT poll squeue (G8)
    cluster.use_nodes = False
    cluster.use_cpu = False
    cluster.use_time = False
    cluster.use_account = False
    cluster.job_number = 1                # 1 structure per job
    cluster.batch_size = 2                # 2 ab-initio structures per learning cycle
    cluster.binary = calc.command         # python -m sscha.ASEClusterRunner ...
    cluster.mpi_cmd = ""                  # the runner is a serial process
    cluster.setup_workdir()               # mkdir -p workdir (locally)

    minim = sscha.SchaMinimizer.SSCHA_Minimizer(ensemble)
    minim.set_minimization_step(0.01)
    minim.meaningful_factor = 1e-3
    minim.kong_liu_ratio = 0.5

    relax = sscha.Relax.SSCHA(minimizer=minim, ase_calculator=calc,
                              cluster=cluster, N_configs=8, max_pop=3,
                              save_ensemble=False)
    relax.relax(get_stress=True)
    relax.minim.finalize()
    relax.minim.dyn.save_qe("sscha_gold_otf_dyn")


if __name__ == "__main__":
    main_localcluster()
