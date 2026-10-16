"""Tests for the backend-agnostic OTF helpers on the base Ensemble class.

.. note:: The ensemble is built inside each test (not in a pytest fixture)
   because the flare C extension segfaults when SGP objects are created in
   a fixture and later used in the test body (a known flare/pytest
   interaction bug). Creating the objects inside the test function works
   reliably.
"""
import os
import numpy as np
import pytest

flare = pytest.importorskip("flare")  # G9: skip the whole module without flare

from flare.bffs.sgp.calculator import SGP_Calculator

import cellconstructor as CC, cellconstructor.Phonons
import sscha, sscha.Ensemble
from ase.calculators.emt import EMT
from ase.build import bulk


def get_gold_dyn(supercell=(2, 2, 2)):
    """Harmonic Au (EMT), 8 atoms with the default supercell."""
    struct = CC.Structure.Structure()
    struct.generate_from_ase_atoms(bulk("Au", "fcc", a=4.0782, cubic=False))
    dyn = CC.Phonons.compute_phonons_finite_displacements(struct, EMT(), supercell=supercell)
    dyn.Symmetrize()
    dyn.ForcePositiveDefinite()
    return dyn


# Keep kernel objects alive for the lifetime of the process: SGP_Wrapper does
# not store a Python reference to the kernels (only to the descriptor
# calculators), so without this the pybind11 wrapper is garbage-collected when
# the helper returns, leaving the C++ SparseGP with a dangling pointer →
# segfault in add_training_structure.
_KERNEL_REFS = []


def get_sgp_calc_au():
    """Empty SGP calculator for gold."""
    from flare.bffs.sgp._C_flare import NormalizedDotProduct, B2
    from flare.bffs.sgp import SGP_Wrapper
    cutoff = 4.0
    kernel = NormalizedDotProduct(2.0, 2)
    b2 = B2("chebyshev", "quadratic", [0.0, cutoff], [], [1, 4, 3], cutoff * np.ones((1, 1)))
    sgp = SGP_Wrapper([kernel], [b2], cutoff, 0.1, 0.1, 0.1, {79: 0},
                      single_atom_energies={0: 0.0}, variance_type="local",
                      opt_method="L-BFGS-B", max_iterations=5)
    _KERNEL_REFS.append(kernel)  # flare bug workaround
    return SGP_Calculator(sgp)


def _make_ensemble(tmp_path, n_configs=8):
    """Build an OTF ensemble for gold+EMT (call inside the test, not a fixture)."""
    np.random.seed(0)
    ens = sscha.Ensemble.Ensemble(get_gold_dyn(), 300)
    ens.generate(n_configs)
    ens.set_otf(get_sgp_calc_au(), std_tolerance_factor=100,
                max_atoms_added=-1, update_style="add_n", update_threshold=None,
                train_hyps=(1, np.inf), output_name=str(tmp_path / "otf_run"))
    return ens


def test_split_propagates_otf_state(tmp_path):
    """Phase 2 regression: get_noncomputed() must carry the OTF state."""
    ensemble = _make_ensemble(tmp_path)
    sub = ensemble.get_noncomputed()
    assert sub.gp_model is ensemble.gp_model
    assert sub.flare_calc is ensemble.flare_calc
    assert sub.std_tolerance == ensemble.std_tolerance
    assert sub.max_atoms_added == ensemble.max_atoms_added
    assert sub.update_style == ensemble.update_style
    assert sub.train_hyps == ensemble.train_hyps
    assert sub.output is ensemble.output
    assert sub.flare_name == ensemble.flare_name


def test_otf_setup_defaults(tmp_path):
    ensemble = _make_ensemble(tmp_path)
    ensemble._otf_setup_defaults(8)
    assert ensemble.max_atoms_added == 8
    assert ensemble.init_atoms == list(range(8))


def test_update_and_predict_cycle(tmp_path):
    """update GP with one EMT structure -> train -> predict the rest."""
    ensemble = _make_ensemble(tmp_path)
    ensemble._otf_setup_defaults(8)

    # empty model: nothing is predicted
    remaining = list(range(ensemble.N))
    ensemble._otf_predict(ensemble.structures, remaining)
    assert remaining == list(range(ensemble.N))

    # update with one EMT reference (empty model -> all init_atoms added)
    struct = ensemble.structures[0]
    atoms = struct.get_ase_atoms()
    atoms.calc = EMT()
    ensemble._otf_update_from_structure(
        struct, atoms.get_potential_energy(), atoms.get_forces(),
        atoms.get_stress())
    assert len(ensemble.gp_model.training_data) >= 1

    ensemble._otf_maybe_train_and_write(dft_counts=1)
    assert os.path.exists(ensemble.flare_name)

    # trained model with huge tolerance: everything predicted
    remaining = list(range(ensemble.N))
    ensemble._otf_predict(ensemble.structures, remaining)
    assert remaining == []
    assert all(ensemble.force_computed)
    assert np.all(np.isfinite(ensemble.energies))
    assert np.all(np.isfinite(ensemble.forces))


def test_clean_runs_generic_label(tmp_path, capsys):
    ensemble = _make_ensemble(tmp_path)
    ensemble.force_computed[:] = True
    ensemble._clean_runs(dft_counts=2)
    out, _ = capsys.readouterr()
    assert "SUMMARY CALCULATIONS" in out
    assert "Steps using OTF-ML model :" in out
