"""Regression tests for the opt-in qspace_light Ensemble mode.

Covers:
- standard 3D light mode does not retain (3N,3N) supercell pols and forbids the
  real-space update_weights (clear error instead of silent misbehavior);
- default mode is unchanged (eager pols, real-space update works);
- light update_weights_fourier reproduces the default fourier rho BITWISE,
  also when the new dyn reorders the sorted supercell modes (the ASR mode
  maps must be refreshed together with the dyn they describe);
- the Gamma-block ASR mask equals Structure.get_asr_modes ground truth;
- 1D systems (one_dim_axis set) take the legacy rotational ASR branch in
  light mode too (fallback with transient full pols);
- the fourier sscha_energies (init / init_from_structures) match the
  real-space formula in value, sign and units;
- _refresh_qspace_caches_from_real_space keeps the current (possibly
  relaxed) centroids and has no side effects on the real-space state;
- the provenance marker goes dirty when an update is interrupted midway;
- get_fourier_forces refuses a q grid where Gamma is not the first point.
"""
import os
import pytest
import numpy as np
import cellconstructor as CC
import cellconstructor.Phonons
import sscha, sscha.Ensemble

_JULIA = bool(getattr(sscha.Ensemble, "__JULIA_EXT__", False))
needs_julia = pytest.mark.skipif(
    not _JULIA, reason="Julia Fourier backend not available")


def _get_dyn(supercell=(3, 3, 3)):
    """Small 2-atom EMT crystal with a real q grid (deterministic)."""
    import ase.calculators.emt
    np.random.seed(0)
    struct = CC.Structure.Structure(2)
    a_param = 4.0
    struct.unit_cell = np.eye(3) * a_param
    struct.atoms[0] = "Au"
    struct.atoms[1] = "Ag"
    struct.coords[1, :] = np.ones(3) * a_param / 2 + 0.2
    struct.build_masses()
    calc = ase.calculators.emt.EMT()
    dyn = CC.Phonons.compute_phonons_finite_displacements(
        struct, calc, supercell=supercell)
    dyn.AdjustQStar()
    dyn.Symmetrize()
    dyn.ForcePositiveDefinite()
    return dyn


def _reordered_dyn(dyn, factor=25.0):
    """Copy of dyn with one non-Gamma q block scaled: same ASR, but the
    global sorted order of the supercell modes interleaves differently."""
    new_dyn = dyn.Copy()
    iq = int(np.argmax([np.linalg.norm(q) for q in dyn.q_tot]))
    assert np.linalg.norm(dyn.q_tot[iq]) > 1e-8
    new_dyn.dynmats[iq][:, :] *= factor
    return new_dyn


def test_light_no_pols_and_rs_update_raises():
    dyn = _get_dyn()
    T = 300.0
    np.random.seed(1)
    ens = sscha.Ensemble.Ensemble(dyn, T, qspace_light=True)
    ens.generate(4)
    assert ens.pols_0 is None
    assert ens.current_pols is None
    with pytest.raises(RuntimeError):
        ens.update_weights(dyn.Copy(), T)
    # still no (3N,3N) after the failed call
    assert ens.pols_0 is None
    assert ens.current_pols is None


def test_default_mode_unchanged():
    dyn = _get_dyn()
    T = 300.0
    np.random.seed(1)
    ens = sscha.Ensemble.Ensemble(dyn, T)
    ens.generate(4)
    nat_sc = dyn.structure.N_atoms * int(np.prod(dyn.GetSupercell()))
    assert ens.pols_0 is not None and ens.pols_0.shape == (3 * nat_sc, 3 * nat_sc)
    assert ens.current_pols is not None
    ens.update_weights(_reordered_dyn(dyn, 1.1), T)
    assert ens._last_weight_update_fourier is False
    ens.__dict__.pop("_last_weight_update_fourier")  # emulate a pickle made before the Codex marker
    ens.update_weights(_reordered_dyn(dyn, 1.05), T)
    assert ens._last_weight_update_fourier is False
    assert np.all(np.isfinite(ens.rho))


@needs_julia
def test_light_ft_rho_bitwise_vs_default_ft():
    dyn = _get_dyn()
    T = 300.0
    new_dyn = _reordered_dyn(dyn, 25.0)  # reorders the sorted supercell modes

    np.random.seed(2)
    ens_d = sscha.Ensemble.Ensemble(dyn, T)
    ens_d.generate(8)
    np.random.seed(2)
    ens_l = sscha.Ensemble.Ensemble(dyn, T, qspace_light=True)
    ens_l.generate(8)

    ens_d.update_weights_fourier(new_dyn, T)
    ens_l.update_weights_fourier(new_dyn, T)
    assert ens_d._last_weight_update_fourier is True
    assert ens_l._last_weight_update_fourier is True

    assert np.array_equal(ens_l.rho, ens_d.rho), \
        "light fourier rho differs from the default fourier rho"
    assert ens_l.pols_0 is None and ens_l.current_pols is None


@needs_julia
def test_light_asr_mask_matches_get_asr_modes():
    dyn = _get_dyn()
    T = 300.0
    new_dyn = _reordered_dyn(dyn, 25.0)

    np.random.seed(3)
    ens = sscha.Ensemble.Ensemble(dyn, T, qspace_light=True)
    ens.generate(4)
    ens.update_weights_fourier(new_dyn, T)

    # The stored maps must be coherent with current_dyn (freshly recomputed)
    _, iq_ref, band_ref, _, _ = ens.current_dyn.DiagonalizeSupercell(q_only=True)
    assert np.array_equal(ens.mode_iq_current, iq_ref)
    assert np.array_equal(ens.mode_band_current, band_ref)

    # Gamma-block ASR mask == full get_asr_modes ground truth
    mask_light = ens._get_asr_modes_gamma(
        ens.pols_q_current, ens.mode_iq_current, ens.mode_band_current,
        ens.current_dyn)
    full_pols = ens.current_dyn.DiagonalizeSupercell(return_qmodes=True)[1]
    super_struct = ens.current_dyn.structure.generate_supercell(
        ens.current_dyn.GetSupercell())
    mask_truth = super_struct.get_asr_modes(full_pols)
    assert np.array_equal(mask_light, mask_truth)


@needs_julia
def test_light_1d_rotational_asr_fallback():
    total_path = os.path.dirname(os.path.abspath(__file__))
    dyn = CC.Phonons.Phonons(os.path.join(total_path, "..", "test_1d_mat",
                                          "1ddyn_asr"))
    dyn.structure.one_dim_axis = 2
    dyn.Symmetrize()
    dyn.ForcePositiveDefinite()

    # Sanity: the 1D system has 4 ASR modes (3 translations + 1 rotation)
    w, p = dyn.DiagonalizeSupercell()
    assert np.sum(dyn.structure.get_asr_modes(p).astype(int)) == 4

    T = 100.0
    new_dyn = dyn.Copy()
    for i in range(len(new_dyn.dynmats)):
        new_dyn.dynmats[i][:, :] *= 1.1

    np.random.seed(4)
    ens_d = sscha.Ensemble.Ensemble(dyn, T)
    ens_d.generate(8)
    np.random.seed(4)
    ens_l = sscha.Ensemble.Ensemble(dyn, T, qspace_light=True)
    ens_l.generate(8)

    ens_d.update_weights(new_dyn, T)          # legacy real-space reference
    ens_l.update_weights_fourier(new_dyn, T)  # light fourier with 1D fallback
    assert ens_d._last_weight_update_fourier is False
    assert ens_l._last_weight_update_fourier is True

    # Without the one_dim_axis fallback the light path would keep the
    # rotational mode in the reweighting (translations-only mask) and rho
    # would be wildly different; with it, only cross-path fp noise remains.
    assert np.allclose(ens_l.rho, ens_d.rho, rtol=1e-6), \
        "1D light fourier rho does not match the real-space reference"


@needs_julia
def test_generate_fourier_sscha_energies_and_forces():
    """After generate(), the fourier-side sscha_energies must equal the
    real-space SSCHA formula (they used to have the wrong sign) and
    sscha_forces_qspace must be stored in Ry/A (a conversion was missing)."""
    dyn = _get_dyn()
    T = 300.0
    np.random.seed(5)
    ens = sscha.Ensemble.Ensemble(dyn, T)
    ens.generate(4)

    e_ref, f_ref = ens.dyn_0.get_energy_forces(None, displacement=ens.u_disps)
    assert np.allclose(ens.sscha_energies, e_ref, rtol=1e-8, atol=1e-12)

    f_back = sscha.Ensemble.julia.Main.vector_q2r(
        ens.sscha_forces_qspace, ens.q_grid, ens.itau, ens.r_lat)
    assert np.allclose(f_back.reshape(f_ref.shape), f_ref,
                       rtol=1e-6, atol=1e-10)


@needs_julia
def test_init_fourier_sscha_energies_match_realspace():
    """Ensemble.init() on the fourier path must compute sscha_energies from
    the SSCHA forces (it used to use the ab-initio forces_qspace)."""
    dyn = _get_dyn()
    T = 300.0
    np.random.seed(6)
    ens = sscha.Ensemble.Ensemble(dyn, T)
    ens.generate(4)
    ens.fourier_gradient = True
    ens.init()
    assert ens._last_weight_update_fourier is False

    e_ref, f_ref = ens.dyn_0.get_energy_forces(None, displacement=ens.u_disps)
    assert np.allclose(ens.sscha_energies, e_ref, rtol=1e-8, atol=1e-12)
    assert np.allclose(ens.sscha_forces, f_ref, rtol=1e-6, atol=1e-10)


@needs_julia
def test_refresh_qspace_after_realspace_update_relaxed_structure():
    """After a real-space update_weights with moved centroids, the targeted
    q-space refresh must describe the displacements w.r.t. the NEW centroids
    and leave rho, sscha_energies and the dyn_0 baseline untouched
    (Ensemble.init() would reset the displacements to dyn_0's centroids and
    overwrite the real-space state)."""
    dyn = _get_dyn()
    T = 300.0
    np.random.seed(7)
    ens = sscha.Ensemble.Ensemble(dyn, T)
    ens.generate(6)

    new_dyn = _reordered_dyn(dyn, 1.2)
    new_dyn.structure.coords[1, :] += 0.03  # relaxed centroid
    ens.update_weights(new_dyn, T)
    assert ens._last_weight_update_fourier is False

    before = {
        "rho": ens.rho.copy(),
        "sscha_energies": ens.sscha_energies.copy(),
        "sscha_forces": ens.sscha_forces.copy(),
        "u_disps": ens.u_disps.copy(),
        "u_disps_original": ens.u_disps_original.copy(),
    }
    ens._refresh_qspace_caches_from_real_space()

    # No side effects on the real-space state
    for key, val in before.items():
        assert np.array_equal(getattr(ens, key), val), \
            "refresh modified ensemble.%s" % key

    # Round trip: the q-space displacements are the CURRENT ones
    u_back = sscha.Ensemble.julia.Main.vector_q2r(
        ens.u_disps_qspace, ens.q_grid, ens.itau, ens.r_lat)
    assert np.allclose(u_back, ens.u_disps, atol=1e-10)
    # ...and NOT the dyn_0-referenced ones
    assert not np.allclose(u_back, ens.u_disps_original, atol=1e-6)

    # SSCHA forces coherent with what update_weights computed (current_dyn)
    f_back = sscha.Ensemble.julia.Main.vector_q2r(
        ens.sscha_forces_qspace, ens.q_grid, ens.itau, ens.r_lat)
    assert np.allclose(f_back.reshape(ens.sscha_forces.shape),
                       ens.sscha_forces, atol=1e-10)

    # The dyn_0 baseline is transformed from u_disps_original, not clobbered
    u0_back = sscha.Ensemble.julia.Main.vector_q2r(
        ens.u_disps_original_qspace, ens.q_grid, ens.itau, ens.r_lat)
    assert np.allclose(u0_back, ens.u_disps_original, atol=1e-10)


@needs_julia
def test_marker_dirty_on_interrupted_update(monkeypatch):
    """If an update dies halfway, the provenance marker must be left dirty so
    that QSpaceLanczos re-derives the caches instead of trusting them."""
    dyn = _get_dyn()
    T = 300.0
    np.random.seed(8)
    ens = sscha.Ensemble.Ensemble(dyn, T)
    ens.generate(4)
    ens.update_weights_fourier(dyn.Copy(), T)
    assert ens._last_weight_update_fourier is True

    def boom(*args, **kwargs):
        raise RuntimeError("interrupted")

    monkeypatch.setattr(sscha.Ensemble.Ensemble, "update_displacements", boom)
    with pytest.raises(RuntimeError):
        ens.update_weights_fourier(_reordered_dyn(dyn, 1.3), T)
    assert ens._last_weight_update_fourier is False


@needs_julia
def test_get_fourier_forces_requires_gamma_first():
    dyn = _get_dyn()
    T = 300.0
    np.random.seed(9)
    ens = sscha.Ensemble.Ensemble(dyn, T)
    ens.generate(4)
    ens.get_fourier_forces(get_error=False)  # Gamma first: fine
    ens.q_grid = np.roll(ens.q_grid, 1, axis=0)
    with pytest.raises(ValueError):
        ens.get_fourier_forces(get_error=False)


@needs_julia
def test_fourier_update_rejects_reordered_q_and_changed_cell():
    """update_weights_fourier contracts dynq with u_disps_qspace index by
    index: a dyn with permuted q_tot or a different unit cell must be
    rejected loudly instead of silently mixing q points."""
    dyn = _get_dyn()
    T = 300.0
    np.random.seed(10)
    ens = sscha.Ensemble.Ensemble(dyn, T)
    ens.generate(4)

    # Permuted q ordering (swap two non-Gamma points)
    perm_dyn = dyn.Copy()
    perm_dyn.q_tot[1], perm_dyn.q_tot[2] = \
        perm_dyn.q_tot[2].copy(), perm_dyn.q_tot[1].copy()
    perm_dyn.dynmats[1], perm_dyn.dynmats[2] = \
        perm_dyn.dynmats[2].copy(), perm_dyn.dynmats[1].copy()
    with pytest.raises(ValueError):
        ens.update_weights_fourier(perm_dyn, T)

    # Changed unit cell (uniform 1% strain)
    cell_dyn = dyn.Copy()
    cell_dyn.structure.unit_cell = cell_dyn.structure.unit_cell * 1.01
    with pytest.raises(ValueError):
        ens.update_weights_fourier(cell_dyn, T)

    # A tiny strain must be rejected too: the guard uses rtol=0, so the
    # stated atol=1e-8 is the actual tolerance (np.allclose's default
    # rtol=1e-5 would silently let ~1e-5 relative strains through)
    tiny_dyn = dyn.Copy()
    tiny_dyn.structure.unit_cell = tiny_dyn.structure.unit_cell * (1 + 3e-6)
    with pytest.raises(ValueError):
        ens.update_weights_fourier(tiny_dyn, T)

    # The unmodified dyn still works after the failed attempts
    ens.update_weights_fourier(_reordered_dyn(dyn, 1.1), T)
    assert np.all(np.isfinite(ens.rho))
