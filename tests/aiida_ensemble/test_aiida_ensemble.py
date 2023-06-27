"""Tests for :mod:`sscha.aiida_ensemble`."""
import pytest


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
    
    