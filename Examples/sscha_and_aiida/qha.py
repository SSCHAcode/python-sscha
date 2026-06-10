"""Script to run QHA with FLARE calculator."""
import numpy as np

from ase import units, Atoms
from ase.io import read
from ase.optimize import BFGS
from ase.constraints import ExpCellFilter

from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms

from flare.bffs.sgp.calculator import SGP_Calculator
from flare.atoms import FLARE_Atoms

structure_filepath     = 'Si.pwi'
model_filepath         = 'model.json'
initial_vcrelax        = True
fmax                   = 0.001
scale_i                = 0.8
scale_f                = 1.2
scale_step             = 12
supercell_matrix       = [4, 4, 4]
distance               = 0.01
t_min                  = 0
t_max                  = 300
symmetrize             = True
primitive_matrix       = None
mesh                   = 100
unitcell0              = read(structure_filepath)

flare_calc, _ = SGP_Calculator.from_file(model_filepath)
unitcell0.calc = flare_calc

def run() -> None:
    """Run QHA calculation."""
    # Initial optimization
    print("Running initial geometry optimization...", flush=True)
    if initial_vcrelax:
        print("Initial volume: ", unitcell0.get_volume(), flush=True)
    unitcell = geometry_optimization(
        unitcell0, 
        vcrelax=initial_vcrelax,
        fmax=fmax,
    )
    if initial_vcrelax:
        print("Final volume: ", unitcell.get_volume(), flush=True)

    # Run equation of state
    print("Running equation of state...", flush=True)
    all_scaled_atoms, energies = run_eos(unitcell)
    
    print("Running phonons for each volume...")
    start_index = -energies.index(min(energies))
    for scaled_atoms in all_scaled_atoms:
        phonons(scaled_atoms, filename=f'./thermal_properties.yaml-{start_index}')
        start_index += 1
    print("Run completed")
    
def run_eos(unitcell: Atoms) -> tuple[tuple[Atoms], list]:
    """Run equation of states."""
    energies, volumes = [], []
    all_scaled_atoms = []
    
    scale_factors = np.linspace(scale_i, scale_f, scale_step)

    for scale_factor in scale_factors:
        scaled_atoms = scale_and_relax(unitcell, scale_factor, fmax=fmax)
        all_scaled_atoms.append(scaled_atoms)
        energies.append(scaled_atoms.get_potential_energy())
        volumes.append(scaled_atoms.get_volume())
    
    np.savetxt('./e-v.dat', np.array([volumes, energies]).T)

    return all_scaled_atoms, energies
    
def phonons(atoms: Atoms, filename: str = None) -> None:
    """Run phonons using Phonopy."""
    unitcell = PhonopyAtoms(
        symbols=atoms.get_chemical_symbols(),
        numbers=atoms.get_atomic_numbers(),
        scaled_positions=atoms.get_scaled_positions(),
        cell=atoms.get_cell(),
    )

    ph = Phonopy(
        unitcell=unitcell,
        primitive_matrix=primitive_matrix,
        supercell_matrix=supercell_matrix,
    )

    ph.generate_displacements(distance=distance)
    supercells = ph.supercells_with_displacements

    sets_of_forces = []
    for supercell in supercells:
        cell, scaled_positions, numbers = supercell.totuple()
        supercell_atoms = Atoms(
            cell=cell,
            scaled_positions=scaled_positions,
            numbers=numbers,
            calculator=flare_calc,
        )
        sets_of_forces.append(supercell_atoms.get_forces())

    ph.forces = sets_of_forces
    ph.produce_force_constants()

    if symmetrize:
        ph.symmetrize_force_constants()
        ph.symmetrize_force_constants_by_space_group()
        
    ph.run_mesh(mesh=mesh)
    ph.run_thermal_properties(t_min=t_min, t_max=t_max)
    ph.thermal_properties.write_yaml(filename)


def geometry_optimization(atoms: Atoms, vcrelax: bool = False, fmax: float = 0.001) -> Atoms:
    """Optimize geometry of the given atoms."""
    if vcrelax:
        ecf = ExpCellFilter(atoms, scalar_pressure=0)
        optimizer = BFGS(ecf)
    else:
        optimizer = BFGS(atoms=atoms)

    optimizer.run(fmax=fmax)

    if vcrelax:
        return optimizer.atoms.atoms
    return optimizer.atoms
   
 
def scale_and_relax(atoms: Atoms, scale_factor: float, fmax: float = 0.001) -> Atoms:
    """Scale and relax the atoms of an ASE Atoms object."""
    ase = atoms.copy()
    ase.calc = atoms.calc
    ase.set_cell(ase.get_cell() * float(scale_factor) ** (1 / 3), scale_atoms=True)
    
    return geometry_optimization(ase, vcrelax=False, fmax=fmax)


if __name__ == "__main__":
    run()
    