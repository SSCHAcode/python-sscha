"""End-to-end OTF on gold+EMT through DirectCluster and LocalCluster."""
import os
import numpy as np
import pytest

flare = pytest.importorskip("flare")  # G9

import cellconstructor as CC
import sscha, sscha.Ensemble
import sscha.Cluster, sscha.LocalCluster, sscha.BaseCluster
import sscha.ClusterCalculators as cluster_calcs
from ase.calculators.emt import EMT

from .test_otf_base import get_gold_dyn, get_sgp_calc_au


def make_direct_cluster(tmp_path, calc, batch_size):
    return sscha.BaseCluster.DirectCluster(batch_size=batch_size, job_number=1), \
        cluster_calcs.ASEDirectCalculator(EMT())


def make_local_cluster(tmp_path, calc, batch_size):
    file_calc = cluster_calcs.ASEFileCalculator(EMT())
    cluster = sscha.LocalCluster.LocalCluster("localhost")
    cluster.workdir = str(tmp_path / "remote")
    cluster.local_workdir = str(tmp_path / "local") + "/"
    cluster.submit_command = "bash"        # blocking local execution (G8)
    cluster.nonblocking_command = False    # no squeue polling (G8)
    cluster.use_nodes = False
    cluster.use_cpu = False
    cluster.use_time = False
    cluster.use_account = False
    cluster.job_number = 1
    cluster.batch_size = batch_size        # learning-cycle size = batch_size*job_number
    cluster.binary = file_calc.command
    cluster.mpi_cmd = ""
    cluster.setup_workdir()
    return cluster, file_calc


@pytest.mark.flare
@pytest.mark.parametrize("backend_factory", [make_direct_cluster, make_local_cluster],
                         ids=["direct", "localcluster"])
def test_cluster_emt_otf(tmp_path, backend_factory):
    np.random.seed(0)
    n_configs, batch_size = 8, 2

    dyn = get_gold_dyn()
    ensemble = sscha.Ensemble.Ensemble(dyn, 300)
    ensemble.generate(n_configs)
    ensemble.set_otf(get_sgp_calc_au(), std_tolerance_factor=100,
                     max_atoms_added=-1, update_style="add_n",
                     update_threshold=None, train_hyps=(1, np.inf),
                     output_name=str(tmp_path / "otf_run"))

    cluster, calc = backend_factory(tmp_path, None, batch_size)

    # Reference EMT values for the structures that MUST be computed ab-initio
    # (first cycle: model empty -> exactly the first batch_size structures)
    ref = []
    for i in range(batch_size):
        atoms = ensemble.structures[i].get_ase_atoms()
        atoms.calc = EMT()
        ref.append((atoms.get_potential_energy(), atoms.get_forces().copy()))

    ensemble.compute_ensemble(calc, compute_stress=True, cluster=cluster)

    # 1. Everything computed (ab-initio or predicted)
    assert all(ensemble.force_computed)
    assert all(ensemble.stress_computed)
    assert np.all(np.isfinite(ensemble.energies))
    assert np.all(np.isfinite(ensemble.forces))
    assert np.all(np.isfinite(ensemble.stresses))

    # 2. The model was trained on the first-cycle structures. With a high
    #    std_tolerance, only the first structure (empty model -> all init_atoms
    #    added) necessarily enters the training set; subsequent structures may
    #    be within bounds and not add new training data.
    assert len(ensemble.gp_model.training_data) >= 1
    assert os.path.exists(ensemble.flare_name)

    # 3. STRONG REGRESSION: the ab-initio structures carry the exact EMT
    #    values (validates the whole chain: calculator -> raw results ->
    #    ingestion -> unit conversions). Ensemble units are Ry, Ry/Ang.
    from ase.units import create_units
    Ry = create_units("2006")["Ry"]
    for i in range(batch_size):
        assert ensemble.energies[i] == pytest.approx(ref[i][0] / Ry, rel=1e-8)
        assert np.allclose(ensemble.forces[i], ref[i][1] / Ry, atol=1e-8)

    # 4. A second ensemble reuses the trained model: zero ab-initio calls
    np.random.seed(1)
    ensemble2 = sscha.Ensemble.Ensemble(dyn, 300)
    ensemble2.generate(4)
    ensemble2.set_otf(ensemble.flare_calc, std_tolerance_factor=100,
                      max_atoms_added=-1, update_style="add_n",
                      update_threshold=None, train_hyps=(1, np.inf),
                      output_name=str(tmp_path / "otf_run2"))
    training_before = len(ensemble2.gp_model.training_data)
    ensemble2.compute_ensemble(calc, compute_stress=True, cluster=cluster)
    assert all(ensemble2.force_computed)
    assert len(ensemble2.gp_model.training_data) == training_before  # nothing new learned
