# -*- coding: utf-8 -*-
from __future__ import print_function

"""
Test the gradient computed splitting the configurations among processors
(get_preconditioned_gradient_parallel).

Run it also in parallel with:
    mpirun -np 4 python test_parallel_gradient.py
"""

import os
import numpy as np

import cellconstructor as CC, cellconstructor.Phonons, cellconstructor.Settings
import sscha, sscha.Ensemble, sscha.Parallel
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

    # Move the dynamical matrix away from dyn_0, so that the weights are not trivial
    new_dyn = dyn.Copy()
    for iq in range(len(new_dyn.q_tot)):
        new_dyn.dynmats[iq] *= 1.05
    ensemble.update_weights(new_dyn, 0)
    return ensemble


def check_parallel_gradient(ensemble, preconditioned):
    grad, err = ensemble.get_preconditioned_gradient(True, True, preconditioned=preconditioned)
    grad_par, err_par = ensemble.get_preconditioned_gradient_parallel(True, True, preconditioned=preconditioned)

    delta = np.max(np.abs(grad - grad_par)) / np.max(np.abs(grad))
    assert delta < 1e-8, "Error, the parallel gradient differs from the serial one (relative {})".format(delta)

    # The error must be a real stochastic error (not a constant), of the same order of the serial one
    assert np.all(np.isfinite(err_par))
    norm_err = np.sqrt(np.sum(np.abs(err)**2))
    norm_err_par = np.sqrt(np.sum(np.abs(err_par)**2))
    assert not np.allclose(np.abs(err_par), 1), "Error, the error of the gradient is a constant"
    assert 0.2 < norm_err_par / norm_err < 5, "Error, wrong error of the parallel gradient ({} vs {})".format(norm_err_par, norm_err)
    return delta, norm_err, norm_err_par


def test_parallel_gradient(monkeypatch=None):
    """
    Split the configurations in many chunks (each chunk is computed as by a different processor)
    and check that the gradient matches the serial one.
    """
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    ensemble = get_ensemble()

    # Emulate 4 processors: in serial, GoParallel computes all the chunks and sums them
    def split_in_four(n_configs):
        edges = np.linspace(0, n_configs, 5).astype(int)
        return [(edges[i], edges[i+1]) for i in range(4)]

    if monkeypatch is not None:
        monkeypatch.setattr(CC.Settings, "split_configurations", split_in_four)

    for preconditioned in [1, 0]:
        check_parallel_gradient(ensemble, preconditioned)


if __name__ == "__main__":
    # Real parallel check (run with mpirun)
    ensemble = get_ensemble()
    for preconditioned in [1, 0]:
        delta, err, err_par = check_parallel_gradient(ensemble, preconditioned)
        sscha.Parallel.pprint("NPROC = {} preconditioned = {}: relative difference = {:.2e}, |err| = {:.4e} (serial {:.4e})".format(
            CC.Settings.GetNProc(), preconditioned, delta, err_par, err))
