"""Tests for :mod:`sscha.aiida_ensemble`."""
import pytest
import numpy as np

from ase import Atoms
from ase.calculators.lj import LennardJones

from flare.atoms import FLARE_Atoms
from flare.bffs.sgp.calculator import SGP_Calculator
from .get_sgp import get_sgp_calc, get_random_atoms, get_empty_sgp


@pytest.fixture
def generate_structure():
    """Return an :class:`cellconstructor.Structure.Structure` instance."""
    
    def _generate_structure(a=2.0, sc_size=2, numbers=[6, 8]):
        """Return an :class:`cellconstructor.Structure.Structure` instance."""
        import numpy as np
        from ase.build import make_supercell
        from cellconstructor.Structure import Structure

        cell = np.eye(3) * a
        positions = np.array([[0, 0, 0], [0.5 , 0.5, 0.5]])
        unit_cell = Atoms(cell=cell, scaled_positions=positions, numbers=numbers, pbc=True)
        multiplier = np.identity(3) * sc_size
        atoms = make_supercell(unit_cell, multiplier)

        structure = Structure()
        structure.generate_from_ase_atoms(atoms)

        return structure
        
    return _generate_structure


@pytest.fixture
def generate_ensemble(generate_structure):
    """Return an AiiDAEnsemble instance."""
    
    def _generate_ensemble(temperature: float = 0.0):
        """Return an AiiDAEnsemble instance."""
        from cellconstructor.Phonons import compute_phonons_finite_displacements
        from sscha.aiida_ensemble import AiiDAEnsemble

        dyn = compute_phonons_finite_displacements(generate_structure(), LennardJones(), supercell=[1,1,1])
        dyn.Symmetrize()
        dyn.ForcePositiveDefinite()

        return AiiDAEnsemble(dyn, temperature)
    
    return _generate_ensemble


def test_no_otf(generate_ensemble):
    """Test the :class:`sscha.aiida_ensemble.AiiDAEnsemble` initialization parameters w/o OTF.
    
    .. note:: this test shouldn't be here, but in the tests for generic Ensemble.
    """
    ensemble = generate_ensemble()
    
    assert ensemble.gp_model is None
    assert ensemble.flare_calc is None
    assert ensemble.std_tolerance is None
    assert ensemble.max_atoms_added is None
    assert ensemble.update_style is None
    assert ensemble.update_threshold is None
    assert ensemble.build_mode is None
    assert ensemble.output is None
    assert ensemble.output_name is None
    assert ensemble.checkpt_name is None
    assert ensemble.flare_name is None
    assert ensemble.atoms_name is None
    assert ensemble.checkpt_files is None
    assert ensemble.write_model is None 
    assert ensemble.init_atoms is None 
    assert ensemble.train_hyps is None 


def test_set_otf(generate_ensemble):
    """Test the :func:`sscha.aiida_ensemble.AiiDAEnsemble.set_otf` method."""
    ensemble = generate_ensemble()
    flare_calc = get_sgp_calc()
    ensemble.set_otf(flare_calc, max_atoms_added=-1)
    
    assert ensemble.gp_model is not None
    assert ensemble.train_hyps == (100,120)


def test_compute_properties(generate_ensemble):
    """Test the :func:`sscha.aiida_ensemble.AiiDAEnsemble._compute_properties` method."""
    ensemble = generate_ensemble()
    ensemble.generate(20)
    atoms = ensemble.structures[0].get_ase_atoms()
    ensemble.set_otf(get_sgp_calc(), max_atoms_added=-1)
    
    ensemble._compute_properties(atoms)


def test_predict_with_model(generate_ensemble):
    """Test the :func:`sscha.aiida_ensemble.AiiDAEnsemble._predict_with_model` method."""
    num_configs = 20
    ensemble = generate_ensemble()
    ensemble.generate(num_configs)
    ensemble.set_otf(get_sgp_calc(), max_atoms_added=-1)
    
    dft_indices = np.arange(0,num_configs,1).tolist()
    
    ensemble._predict_with_model(
        ensemble.structures,
        dft_indices
    )
    
    assert dft_indices == []
    assert len(ensemble.force_computed) == num_configs
    assert all(ensemble.force_computed)


def test_write_model(generate_ensemble):
    """Test the :func:`sscha.aiida_ensemble.AiiDAEnsemble._write_model` method.
    
    .. note:: in principle here we should double check with a regression test that
        the json file is formatted as expected.
    """
    ensemble = generate_ensemble()
    ensemble.set_otf(get_sgp_calc(), max_atoms_added=-1)
    
    ensemble._write_model()


@pytest.mark.parametrize('flare_calc', (get_sgp_calc(), SGP_Calculator(get_empty_sgp())))
def test_update_gp(generate_ensemble, flare_calc):
    """Test the :func:`sscha.aiida_ensemble.AiiDAEnsemble._update_gp` method."""
    ensemble = generate_ensemble()
    ensemble.set_otf(flare_calc, max_atoms_added=-1)
    ensemble.init_atoms = [1]
    
    atoms = get_random_atoms()
    atoms.calc = LennardJones()
    forces = atoms.get_forces()
    energy = atoms.get_potential_energy()
    stress = atoms.get_stress()

    ensemble._update_gp(
        FLARE_Atoms.from_ase_atoms(atoms),
        forces,
        energy,
        stress
    )

    assert len(ensemble.gp_model.training_data) in [1, 2]

def test_train_gp(generate_ensemble):
    """Test the :func:`sscha.aiida_ensemble.AiiDAEnsemble._train_gp` method."""
    ensemble = generate_ensemble()
    ensemble.set_otf(get_sgp_calc(), max_atoms_added=-1)

    ensemble._train_gp()