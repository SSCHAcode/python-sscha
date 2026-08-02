"""Tests for the layered Cluster/calculator interface (no flare needed)."""
import os
import threading
import subprocess

import numpy as np
import pytest

import cellconstructor as CC
from cellconstructor.calculators import Espresso
from ase.build import bulk
from ase.calculators.emt import EMT

import sscha, sscha.Ensemble
import sscha.Cluster, sscha.LocalCluster, sscha.BaseCluster
import sscha.ClusterCalculators as cluster_calcs

from .test_otf_base import get_gold_dyn  # dyn fixture builder (no flare use)


def get_gold_structure():
    struct = CC.Structure.Structure()
    struct.generate_from_ase_atoms(bulk("Au", "fcc", a=4.0782, cubic=False))
    return struct


class MockDirectCalculator(cluster_calcs.DirectCalculator):
    """Deterministic fake: energy = index, zero forces/stress, echo structure."""
    def copy(self):
        return MockDirectCalculator()
    def compute(self, structure):
        nat = structure.N_atoms
        return {"energy": 1.0, "forces": np.zeros((nat, 3)),
                "stress": np.zeros(6), "structure": structure.copy(),
                "mock_extra": 42}


def test_directcluster_driver_no_otf():
    """The generic driver fills the ensemble through a DirectCalculator."""
    np.random.seed(0)
    ensemble = sscha.Ensemble.Ensemble(get_gold_dyn(), 0)
    ensemble.generate(4)
    cluster = sscha.BaseCluster.DirectCluster(batch_size=2, job_number=1)
    ensemble.compute_ensemble(MockDirectCalculator(), compute_stress=True,
                              cluster=cluster)
    assert all(ensemble.force_computed)
    assert all(ensemble.stress_computed)
    assert np.allclose(ensemble.energies, 1.0 / 13.605698066)  # eV -> Ry
    assert all(p["mock_extra"] == 42 for p in ensemble.all_properties)


def test_pairing_rules():
    """DirectCluster rejects file calculators; Cluster rejects direct ones."""
    cluster = sscha.BaseCluster.DirectCluster()
    with pytest.raises(TypeError, match="cannot be used"):
        cluster._check_calculator_interface(Espresso(
            input_data={"control": {}, "system": {}}, pseudopotentials={"Au": "Au.upf"}))

    remote = sscha.Cluster.Cluster(hostname="localhost")
    with pytest.raises(TypeError, match="cannot be used"):
        remote._check_calculator_interface(MockDirectCalculator())

    # correct pairings pass
    cluster._check_calculator_interface(MockDirectCalculator())
    remote._check_calculator_interface(Espresso(
        input_data={"control": {}, "system": {}}, pseudopotentials={"Au": "Au.upf"}))


def test_espresso_prepare_input_file_backward_compatible(tmp_path):
    """Phase 3.2 regression: QE calculators must behave exactly as before."""
    cluster = sscha.LocalCluster.LocalCluster("localhost")
    cluster.local_workdir = str(tmp_path) + "/"
    cluster.lock = threading.Lock()
    calc = Espresso(
        input_data={
            "control": {"tprnfor": True, "tstress": True},
            "system": {"ecutwfc": 30, "ecutrho": 240, "occupations": "fixed"},
            "electrons": {"conv_thr": 1e-8},
        },
        pseudopotentials={"Au": "Au.upf"},
        kpts=(1, 1, 1),
    )
    inputs, outputs = cluster.prepare_input_file([get_gold_structure()], calc, ["ESP_0"])
    assert inputs == ["ESP_0.pwi"]
    assert outputs == ["ESP_0.pwo"]
    assert os.path.exists(tmp_path / "ESP_0.pwi")


def test_espresso_calculator_adapter():
    """The espresso-specific calculator keeps class through copy()."""
    calc = cluster_calcs.EspressoCalculator(
        input_data={"control": {}, "system": {}}, pseudopotentials={"Au": "Au.upf"})
    assert calc.input_extension == ".pwi"
    assert calc.output_extension == ".pwo"
    assert "PREFIX.pwi" in calc.get_execution_command()
    assert isinstance(calc.copy(), cluster_calcs.EspressoCalculator)


def test_ase_file_calculator_prepare_input(tmp_path):
    """The ASE file bridge uses its own extensions and ships calculator.pkl."""
    cluster = sscha.LocalCluster.LocalCluster("localhost")
    cluster.local_workdir = str(tmp_path) + "/"
    cluster.lock = threading.Lock()
    calc = cluster_calcs.ASEFileCalculator(EMT())
    inputs, outputs = cluster.prepare_input_file([get_gold_structure()], calc, ["ESP_0"])
    assert inputs == ["ESP_0_input.json", "calculator.pkl"]
    assert outputs == ["ESP_0.json"]
    assert os.path.exists(tmp_path / "ESP_0_input.json")
    assert os.path.exists(tmp_path / "calculator.pkl")


def test_ase_file_calculator_local_roundtrip(tmp_path):
    """Run locally exactly what the cluster would run; compare with bare EMT."""
    calc = cluster_calcs.ASEFileCalculator(EMT())
    calc.set_directory(str(tmp_path))
    calc.set_label("ESP_0")
    struct = get_gold_structure()
    calc.write_input(struct)

    # Execute the runner exactly as Cluster.get_execution_command would
    cmd = calc.command.replace("PREFIX", os.path.join(str(tmp_path), "ESP_0"))
    subprocess.run(cmd, shell=True, check=True, cwd=str(tmp_path))

    calc.read_results()

    atoms = struct.get_ase_atoms()
    atoms.calc = EMT()
    assert calc.results["energy"] == pytest.approx(atoms.get_potential_energy(), rel=1e-10)
    assert np.allclose(calc.results["forces"], atoms.get_forces(), atol=1e-10)
    assert np.allclose(calc.results["stress"], atoms.get_stress(), atol=1e-10)
    assert np.allclose(calc.structure.coords, struct.coords, atol=1e-10)


def test_ase_direct_calculator_matches_bare_ase():
    calc = cluster_calcs.ASEDirectCalculator(EMT())
    struct = get_gold_structure()
    res = calc.compute(struct)
    atoms = struct.get_ase_atoms()
    atoms.calc = EMT()
    assert res["energy"] == pytest.approx(atoms.get_potential_energy(), rel=1e-12)
    assert np.allclose(res["forces"], atoms.get_forces(), atol=1e-12)
    assert np.allclose(res["stress"], atoms.get_stress(), atol=1e-12)
    # per-thread copies are independent
    assert calc.copy().ase_calc is not calc.ase_calc


def test_nonpicklable_calculator_fails_early():
    class Unpicklable:
        def __getstate__(self):
            raise TypeError("no pickle")
    with pytest.raises(TypeError):
        cluster_calcs.ASEFileCalculator(Unpicklable())


def test_directcluster_stress_less_calculator():
    """A calculator without stress works via DirectCluster with get_stress=False."""
    class NoStressEMT(EMT):
        def get_stress(self, atoms=None, voigt=True):
            raise NotImplementedError("NoStressEMT does not support stress")

    np.random.seed(0)
    ensemble = sscha.Ensemble.Ensemble(get_gold_dyn(), 0)
    ensemble.generate(4)
    cluster = sscha.BaseCluster.DirectCluster(batch_size=2, job_number=1)
    ensemble.compute_ensemble(cluster_calcs.ASEDirectCalculator(NoStressEMT()),
                              compute_stress=False, cluster=cluster)

    assert all(ensemble.force_computed)
    assert not any(ensemble.stress_computed)
    from ase.units import create_units
    Ry = create_units("2006")["Ry"]
    atoms = ensemble.structures[0].get_ase_atoms()
    atoms.calc = NoStressEMT()
    assert ensemble.energies[0] == pytest.approx(atoms.get_potential_energy() / Ry, rel=1e-10)
    assert np.allclose(ensemble.forces[0], atoms.get_forces() / Ry, atol=1e-10)


_NO_STRESS_MODULE = """
from ase.calculators.emt import EMT


class NoStressEMT(EMT):
    \"\"\"EMT that refuses to compute the stress (like stress-less codes).\"\"\"
    def get_stress(self, atoms=None, voigt=True):
        raise NotImplementedError("NoStressEMT does not support stress")
"""


def test_ase_file_calculator_stress_less_roundtrip(tmp_path):
    """The file bridge tolerates calculators without stress (get_stress=False).

    The runner must not crash on a missing stress implementation: the output
    JSON simply omits the ``stress`` key and ``read_results`` keeps working.
    """
    import sys as _sys

    (tmp_path / "_nostress_calc.py").write_text(_NO_STRESS_MODULE)
    _sys.path.insert(0, str(tmp_path))
    try:
        from _nostress_calc import NoStressEMT

        calc = cluster_calcs.ASEFileCalculator(NoStressEMT())
        calc.set_directory(str(tmp_path))
        calc.set_label("ESP_0")
        struct = get_gold_structure()
        calc.write_input(struct)

        # Execute the runner exactly as Cluster.get_execution_command would
        cmd = calc.command.replace("PREFIX", os.path.join(str(tmp_path), "ESP_0"))
        subprocess.run(cmd, shell=True, check=True, cwd=str(tmp_path))

        calc.read_results()
        assert "energy" in calc.results
        assert "forces" in calc.results
        assert "stress" not in calc.results

        atoms = struct.get_ase_atoms()
        atoms.calc = NoStressEMT()
        assert calc.results["energy"] == pytest.approx(atoms.get_potential_energy(), rel=1e-10)
        assert np.allclose(calc.results["forces"], atoms.get_forces(), atol=1e-10)
    finally:
        _sys.path.remove(str(tmp_path))
