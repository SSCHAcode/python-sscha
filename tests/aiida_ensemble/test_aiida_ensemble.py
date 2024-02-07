"""Tests for :mod:`sscha.aiida_ensemble`."""
import pytest
import numpy as np

from sscha.aiida_ensemble import AiiDAEnsemble


def get_ensemble() -> AiiDAEnsemble:
    """Return an AiiDAEnsemble instance."""
    import os
    from cellconstructor.Phonons import Phonons

    path = os.path.dirname(os.path.abspath(__file__))

    return AiiDAEnsemble(Phonons(os.path.join(path,"dyn"), 3), 0, (2,1,2))


def test_clean_runs():
    """Test the :func:`sscha.aiida_ensemble.AiiDAEnsemble._clean_runs` method."""
    ensemble = get_ensemble()
    num_configs, num_atoms = 4, 1
    ensemble.generate(num_configs)

    ensemble.energies = np.ones((num_configs,)) # (configs,)
    ensemble.forces = np.ones((num_configs, num_atoms, 3)) # (configs, atoms, force index)
    ensemble.stresses = np.ones((num_configs, 3, 3)) # (configs, 3, 3)
    ensemble.force_computed = np.array([True, False, True, True], dtype=bool)
    ensemble.stress_computed = np.copy(ensemble.force_computed)
    ensemble._clean_runs()

    assert all(ensemble.force_computed)
    assert len(ensemble.force_computed) == 3
    assert len(ensemble.stress_computed) == 3
    assert np.all(np.isclose(ensemble.forces, np.ones((num_configs-1, num_atoms, 3))))


@pytest.mark.usefixtures('aiida_profile')
def test_get_running_workchains(generate_workchain_pw_node):
    """Test the :func:`sscha.aiida_ensemble.get_running_workchains` method."""
    from plumpy import ProcessState
    from sscha.aiida_ensemble import get_running_workchains

    finished = ProcessState.FINISHED
    running  = ProcessState.RUNNING

    workchains = [
        generate_workchain_pw_node(process_state=finished, exit_status=0, label='T_300_id_0'),
        generate_workchain_pw_node(process_state=finished, exit_status=300, label='T_300_id_1'),
    ]
    success = [False, False]

    wcs_left = get_running_workchains(workchains=workchains, success=success)

    assert not wcs_left
    assert success == [True, False]

    workchains = [
        generate_workchain_pw_node(process_state=running, label='T_300_id_0'),
        generate_workchain_pw_node(process_state=finished, label='T_300_id_1'),
        generate_workchain_pw_node(process_state=running, label='T_300_id_2'),
    ]
    success = [False, False, False]

    wcs_left = get_running_workchains(workchains=workchains, success=success)

    assert wcs_left == [workchains[0], workchains[2]]
    assert success == [False, True, False]


@pytest.mark.usefixtures('aiida_profile')
def test_submit_and_get_workchains(fixture_code):
    """Test the :func:`sscha.aiida_ensemble.submit_and_get_workchains` method."""
    from cellconstructor.Structure import Structure
    from sscha.aiida_ensemble import submit_and_get_workchains

    pw_code = fixture_code('quantumespresso.pw')
    structures = [Structure(nat=1) for _ in range(5)]

    workchains = submit_and_get_workchains(
        structures=structures,
        pw_code=pw_code,
        temperature=300,
    )

    assert len(workchains) == 5
