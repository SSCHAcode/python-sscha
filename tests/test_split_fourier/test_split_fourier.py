# -*- coding: utf-8 -*-
from __future__ import print_function

import os
import pytest
import numpy as np

import cellconstructor as CC, cellconstructor.Phonons
import sscha, sscha.Ensemble
import ase.calculators.emt


def get_ensemble(n_configs=40):
    """Build a small ensemble of a low symmetry Au-Ag crystal with EMT."""
    np.random.seed(0)

    struct = CC.Structure.Structure(2)
    a_param = 4
    struct.unit_cell = np.eye(3) * a_param
    struct.atoms[0] = "Au"
    struct.atoms[1] = "Ag"
    struct.coords[1, :] = np.ones(3) * a_param / 2 + 0.2
    struct.build_masses()

    calculator = ase.calculators.emt.EMT()
    dyn = CC.Phonons.compute_phonons_finite_displacements(struct, calculator, supercell=(2, 2, 2))
    dyn.AdjustQStar()
    dyn.Symmetrize()
    dyn.ForcePositiveDefinite()

    ensemble = sscha.Ensemble.Ensemble(dyn, 0)
    ensemble.generate(n_configs)
    ensemble.compute_ensemble(calculator)
    return ensemble


@pytest.mark.julia
def test_split_fourier():
    """
    The split ensemble must have the same q space forces and gradient
    of the original one (issue #436)
    """
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    ensemble = get_ensemble()
    assert ensemble.fourier_gradient

    # Move the dynamical matrix away from dyn_0, as during a minimization
    new_dyn = ensemble.current_dyn.Copy()
    for iq in range(len(new_dyn.q_tot)):
        new_dyn.dynmats[iq] *= 1.05
    ensemble.update_weights_fourier(new_dyn, 0)

    splitted = ensemble.split(np.ones(ensemble.N, dtype=bool))

    assert np.max(np.abs(splitted.forces_qspace)) > 0, "Error, the split ensemble has no forces in q space"
    assert np.allclose(splitted.forces_qspace, ensemble.forces_qspace)
    assert np.allclose(splitted.sscha_forces_qspace, ensemble.sscha_forces_qspace)
    assert np.allclose(splitted.rho, ensemble.rho)

    grad, err = ensemble.get_fourier_gradient()
    grad_split, err_split = splitted.get_fourier_gradient()
    assert np.allclose(grad, grad_split), "Error, the gradient of the split ensemble is wrong"
    assert np.isclose(err, err_split)


@pytest.mark.julia
def test_fourier_sscha_energies():
    """
    The sscha energies computed in fourier space (init and init_from_structures)
    must match the real space ones.
    """
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    ensemble = get_ensemble()
    assert ensemble.fourier_gradient

    energies_real, _ = ensemble.dyn_0.get_energy_forces(None, displacement=ensemble.u_disps)

    ensemble.init()
    assert np.allclose(ensemble.sscha_energies, energies_real), "Error, wrong sscha energies after init"

    new_ensemble = sscha.Ensemble.Ensemble(ensemble.dyn_0, 0)
    new_ensemble.init_from_structures(ensemble.structures)
    assert np.allclose(new_ensemble.sscha_energies, energies_real), \
        "Error, wrong sscha energies after init_from_structures"


if __name__ == "__main__":
    test_split_fourier()
    test_fourier_sscha_energies()
