# -*- coding: utf-8 -*-
"""Module for handling automated calculation via aiida-quantumespresso."""
from __future__ import annotations

from typing import Literal
from copy import copy, deepcopy
import time
import sys

from ase import units
from cellconstructor.Structure import Structure
import numpy as np
from numpy import ndarray

from .Ensemble import Ensemble

try:
    from aiida.orm import WorkChainNode
    from qe_tools import CONSTANTS

    gpa_to_rybohr3 = 1.0 / (CONSTANTS.ry_si / CONSTANTS.bohr_si**3 / 1.0e9)  # GPa -> Ry/Bohr^3
    ase_stress_units = -1.0 * gpa_to_rybohr3 * units.Ry / units.Bohr**3  # convention as in ASE (sign and eV/Ang^3)
except ImportError:
    import warnings
    warnings.warn('aiida or aiida-quantumespresso are not installed')

try:
    from flare.atoms import FLARE_Atoms
    from flare.io.output import compute_mae
    from flare.learners.utils import get_env_indices, is_std_in_bound
except ImportError:
    pass


class AiiDAEnsemble(Ensemble):
    """Ensemble subclass to interface SSCHA with aiida-quantumespresso."""

    def compute_ensemble( # pylint: disable=arguments-renamed
        self,
        pw_code: str,
        protocol: Literal['fast', 'balanced', 'stringent'] = 'balanced',
        options: dict | None  = None,
        overrides: dict | None = None,
        group_label: str | None = None,
        waiting_time: int | float = 2.5,
        batch_number: int = 1,
        check_time: int = 60,
        **kwargs
    ) -> None:
        """Compute ensemble properties.

        Args:
        ----
            pw_code: The string associated with the AiiDA code for `pw.x`
            protocol: The protocol to be used; available protocols are 'fast', 'balanced' and 'stringent'
            options: The options for the calculations, such as the resources, wall-time, etc.
            overrides: The overrides for the :func:`aiida_quantumespresso.workflows.pw.base.PwBaseWorkChain.get_builder_from_protocol`
            group_label: The group label where to add the submitted nodes for eventual future inspection
            waiting_time: Time delay in seconds for WorkChain submission; usefull for many configurations
            batch_number: Number of batches used to split the submission of all the structures, one after the other.
                For example: 2 would submit two batches, computing the first one, then the second.
                This is particularly useful when performing on-the-fly simulations, so that the ML potential
                can be trained on previous batches and (hopefully) predict on the following batches.
            check_time: Seconds to wait before checking the status of the submitted workchains
            kwargs: The kwargs for the get_builder_from_protocol

        """
        from aiida.orm import load_group

        try:
            group = None if group_label is None else load_group(group_label)
        except: # NotExsistent
            from aiida.orm import Group
            group = Group(group_label)
            group.store()

        # Check if not all the calculation needs to be done
        if self.force_computed is None:
            self.force_computed = np.array([False] * self.N, dtype=bool)

        self.has_stress = True  # by default we calculate stresses with the `get_builder_from_protocol`
        if overrides:
            try:
                tstress = overrides['pw']['parameters']['CONTROL']['tstress']
                self.has_stress = tstress
            except KeyError:
                pass

        structures = copy(self.structures)
        dft_counts = 0
        dft_indices_batches = split_array(list(range(len(structures))), batch_number)  # store here the indices to run with DFT/AiiDA
        
        if batch_number > 1:
            print(f"Submission in batches is active. Number of batches that will be submitted: {batch_number}")
        
        for batch_n, dft_indices in enumerate(dft_indices_batches):
            dft_indices = dft_indices.tolist()
            if batch_number > 1:
                print(f"Batch submitted: {batch_n+1}/{batch_number}")
            
            # ================ FLARE SECTION ================= #
            # If a model is specified and it's not empty, try to predict.
            # Predict only the ones that are within uncertainty, the rest do via DFT/AiiDA.
            if self.gp_model is not None:
                number_of_atoms = structures[0].get_ase_atoms().get_global_number_of_atoms()
                
                if self.max_atoms_added < 0:
                    self.max_atoms_added = number_of_atoms

                if self.init_atoms is None:
                    self.init_atoms = list(range(number_of_atoms))
                
                if len(self.gp_model.training_data) > 0:
                    self._predict_with_model(structures, dft_indices)
                
            dft_counts += len(dft_indices)

            # ================= AIIDA SECTION  ================ #
            if len(dft_indices) > 0:
                workchains = submit_and_get_workchains(
                    structures=[structures[i] for i in dft_indices],
                    pw_code=pw_code,
                    temperature=self.current_T,
                    dft_indices=dft_indices,
                    protocol=protocol,
                    options=options,
                    overrides=overrides,
                    waiting_time=waiting_time,
                    **kwargs
                )

                sys.stdout.flush()

                if group:
                    group.add_nodes(workchains)

                workchains_copy = copy(workchains)
                while workchains_copy:
                    workchains_copy = get_running_workchains(workchains_copy, self.force_computed)
                    if workchains_copy:
                        time.sleep(check_time)  # wait before checking again
                
                # ================ UPDATE SECTION ================ #
                for i, is_computed in enumerate(self.force_computed):
                    if is_computed and i in dft_indices:
                        dft_stress = None
                        wc = workchains[dft_indices.index(i)]

                        dft_energy = wc.outputs.output_parameters.dict.energy
                        dft_forces = wc.outputs.output_trajectory.get_array('forces')[-1]

                        self.energies[i] = dft_energy / CONSTANTS.ry_to_ev # eV     -> Ry
                        self.forces[i] = dft_forces / CONSTANTS.ry_to_ev   # eV/Ang -> Ry/Ang

                        if self.has_stress:
                            stress = wc.outputs.output_trajectory.get_array('stress')[-1]

                            self.stresses[i, :, :] = stress * gpa_to_rybohr3 # GPa -> Ry/(Bohr^3)

                            dft_stress = ase_stress_units * np.array([
                                stress[0, 0], stress[1, 1], stress[2, 2], 
                                stress[1, 2], stress[0, 2], stress[0, 1],
                            ]) # GPa -> -eV/(Ang^3)

                        if self.gp_model is not None:
                            self._update_gp(
                                FLARE_Atoms.from_ase_atoms(wc.inputs.pw.structure.get_ase()),
                                dft_frcs=dft_forces,
                                dft_energy=dft_energy,
                                dft_stress=dft_stress,
                            )

                # ================ TRAIN SECTION ================ #
                if self.gp_model is not None:
                    if dft_counts > 0:
                        if self.train_hyps[0] <= len(self.gp_model.training_data) <= self.train_hyps[1]:
                            self._train_gp()
                        self._write_model()
                
                sys.stdout.flush()

        # ================ FINALIZE ================ #
        # if self.has_stress:
        #     self.stress_computed = copy(self.force_computed)

        self._clean_runs(dft_counts)
        self.init()

    def _predict_with_model(
        self,
        structures: list[Structure],
        dft_indices: list[int],
    ) -> None:
        """Predict on all the structures and estimate errors.

        This is used to remove the structures indecis to not compute via AiiDA/DFT.

        Args:
        ----
            structures: list of :class:`~cellconstructor.Structure.Structure` to simulate
            dft_indices: list of integers related to the structures

        """
        sub_indices = deepcopy(dft_indices)

        for index in sub_indices:
            structure = structures[index]
            atoms = FLARE_Atoms.from_ase_atoms(structure.get_ase_atoms())
            self._compute_properties(atoms)

            # get max uncertainty atoms
            if self.build_mode == 'bayesian':
                env_selection = is_std_in_bound
            elif self.build_mode == 'direct':
                env_selection = get_env_indices

            tic = time.time()

            std_in_bound, _ = env_selection(
                self.std_tolerance,
                self.gp_model.force_noise,
                atoms,
                max_atoms_added=self.max_atoms_added,
                update_style=self.update_style,
                update_threshold=self.update_threshold,
            )

            self.output.write_wall_time(tic, task='Env Selection')

            if not std_in_bound:
                print(f"[DFT CALLED] For structure with id={index}")
            else:
                print(f"[BFFS  USED] For structure with id={index}")
                dft_indices.remove(index)  # remove index computed via ML-FF
    
                self.energies[index] = deepcopy(atoms.get_potential_energy()) / units.Ry # eV -> Ry
                self.forces[index] = deepcopy(atoms.get_forces()) / units.Ry # eV/Ang -> Ry/Ang
                if self.has_stress:
                    self.stresses[index, :, :] = -1 * deepcopy(atoms.get_stress(voigt=False)) * (units.Bohr**3 / units.Ry) # -eV/(Ang^3) -> Ry/(Bohr^3)
                    self.stress_computed[index] = True

                self.force_computed[index] = True
            
            sys.stdout.flush()
    

    def _compute_properties(self, atoms: FLARE_Atoms) -> None:
        """Compute energies, forces, stresses, and their uncertainties.

        The FLARE ASE calculator is used, and write the results.

        Args:
        ----
            atoms: a :class:`flare.atoms.FLARE_Atoms` instance for which to compute properties

        """
        tic = time.time()

        atoms.calc = self.flare_calc
        atoms.calc.calculate(atoms)

        self.output.write_wall_time(tic, task='Compute Properties')

    def _write_model(self) -> None:
        """Write the current model in a JSON file."""
        self.flare_calc.write_model(self.flare_name)

    def _update_gp(
        self,
        atoms: FLARE_Atoms,
        dft_frcs: ndarray,
        dft_energy: float | None = None,
        dft_stress: ndarray | None = None,
    ) -> None:
        """Update the current GP model.

        Args:
        ----
            atoms (FLARE_Atoms): :class:`flare.atoms.FLARE_Atoms`` instance whose
                local environments will be added to the training set.
            dft_frcs (np.ndarray): DFT forces on all atoms in the structure, in eV/Angstrom.
            dft_energy (float): total energy of the entire structure, in eV.
            dft_stress (np.ndarray): DFT stress on structure.
                Sign as in ASE (-1 in respect with QE), units in eV/Angstrom^3,
                and in Voigt notation, i.e. (xx, yy, zz, yz, xz, xy).

        """
        from ase.calculators.singlepoint import SinglePointCalculator

        tic = time.time()
        is_empty_model = len(self.gp_model.training_data) == 0

        # Here we make the decision to skip adding environments, if the stds 
        # are within the user-defined boundaries, even if the ab-initio calculation
        # was performed. This avoids slowing down the model, while the SSCHA
        # is feeded with the DFT results.
        if is_empty_model:
            std_in_bound = False
            train_atoms = self.init_atoms
        else:
            self._compute_properties(atoms)

            # get max uncertainty atoms
            if self.build_mode == 'bayesian':
                env_selection = is_std_in_bound
            elif self.build_mode == 'direct':
                env_selection = get_env_indices

            tic = time.time()

            std_in_bound, train_atoms = env_selection(
                self.std_tolerance,
                self.gp_model.force_noise,
                atoms,
                max_atoms_added=self.max_atoms_added,
                update_style=self.update_style,
                update_threshold=self.update_threshold,
            )

            self.output.write_wall_time(tic, task='Env Selection')

            # compute mae and write to output
            e_mae, e_mav, f_mae, f_mav, s_mae, s_mav = compute_mae(
                atoms,
                self.output.basename,
                atoms.potential_energy,
                atoms.forces,
                atoms.stress,
                dft_energy,
                dft_frcs,
                dft_stress,
                False,
            )

        if not std_in_bound:
            if not is_empty_model:
                stds = self.flare_calc.results.get('stds', np.zeros_like(dft_frcs))
                self.output.add_atom_info(train_atoms, stds)

            # Convert ASE stress (xx, yy, zz, yz, xz, xy) to FLARE stress
            # (xx, xy, xz, yy, yz, zz).
            flare_stress = None
            if dft_stress is not None:
                flare_stress = -np.array([
                    dft_stress[0],
                    dft_stress[5],
                    dft_stress[4],
                    dft_stress[1],
                    dft_stress[3],
                    dft_stress[2],
                ])

            results = {
                'forces': dft_frcs,
                'energy': dft_energy,
                'free_energy': dft_energy,
                'stress': dft_stress,
            }

            atoms.calc = SinglePointCalculator(atoms, **results)

            # update gp model
            self.gp_model.update_db(
                atoms,
                dft_frcs,
                custom_range=train_atoms,
                energy=dft_energy,
                stress=flare_stress,
            )

            self.gp_model.set_L_alpha()
            self.output.write_wall_time(tic, task='Update GP')
        

    def _train_gp(self) -> None:
        """Optimize the hyperparameters of the current GP model."""
        tic = time.time()

        self.gp_model.train(logger_name=self.output_name + 'hyps.dat')

        self.output.write_wall_time(tic, task='Train Hyps')

        hyps, labels = self.gp_model.hyps_and_labels
        if labels is None:
            labels = self.gp_model.hyp_labels

        self.output.write_hyps(
            labels,
            hyps,
            tic, # actually here there should be the actual start time of the entire simulation
            self.gp_model.likelihood,
            self.gp_model.likelihood_gradient,
            hyps_mask=self.gp_model.hyps_mask,
        )
        
    def _clean_runs(self, dft_counts: int) -> None:
        """Clean the failed runs and print summary.
        
        Args:
        ----
            dft_counts (int): number of performed DFT calculations.
 
        """
        n_calcs = np.sum(self.force_computed.astype(int))
        print('=============== SUMMARY AIIDA CALCULATIONS =============== \n')
        print('Total structures included: ', n_calcs)
        print('Structures not included  : ', self.N-n_calcs)
        if self.gp_model is not None:
            print('Steps using OTF-ML model : ', self.N-dft_counts)
        print()
        print('===================== END OF SUMMARY ===================== \n')
        if n_calcs != self.N:
            self.remove_noncomputed()


def get_running_workchains(workchains: list[WorkChainNode], success: list[bool]) -> list:
    """Get the running workchains popping the finished ones.

    Two extra array should be given to populate the successfully finished runs.

    Args:
    ----
        workchains: list of :class:`~aiida.orm.WorkChainNode`
        success: list where to store whether the workchains finished successfully or not.

    """
    wcs_left = copy(workchains)

    for workchain in workchains:
        if workchain.is_terminated:
            if workchain.is_failed:
                print(f'[FAILURE] for <PwBaseWorkChain> with PK={workchain.pk}')
            else:
                index = int(workchain.label.split('_')[-1])
                success[index] = True
                print(f'[SUCCESS] for <PwBaseWorkChain> with PK={workchain.pk}')

            wcs_left.remove(workchain)  # here it may be critical

    return wcs_left


def submit_and_get_workchains(
    structures: list[Structure],
    pw_code: str,
    temperature: float | int,
    dft_indices: list[int],
    protocol: str = 'moderate',
    options: dict = None,
    overrides: dict = None,
    waiting_time: int | float = 2.5,
    **kwargs
) -> list[WorkChainNode]:
    """Submit and return the workchains for a list of :class:`~cellconstructor.Structure.Structure`.

    Args:
    ----
        structures: a list of :class:`~cellconstructor.Structure.Structure` to run via PwBaseWorkChain.
        pw_code: The string associated with the AiiDA code for `pw.x`
        temperature: The temperature corresponding to the structures ensemble
        dft_indices: The indices of the compute ensemble related to the structures.
        protocol: The protocol to be used; available protocols are 'fast', 'moderate' and 'precise'
        options: The options for the calculations, such as the resources, wall-time, etc.
        overrides: The overrides for the get_builder_from_protocol
        waiting_time: Time delay in seconds for WorkChain submission; usefull for many configurations
        kwargs: The kwargs for the get_builder_from_protocol

    """
    from aiida.engine import submit
    from aiida.orm import StructureData
    from aiida.plugins import WorkflowFactory

    PwBaseWorkChain = WorkflowFactory('quantumespresso.pw.base')

    structures_data = [StructureData(ase=cc.get_ase_atoms()) for cc in structures]
    workchains = []

    for i, structure in zip(dft_indices, structures_data):
        builder = PwBaseWorkChain.get_builder_from_protocol(
            code=pw_code, structure=structure, protocol=protocol, options=options, overrides=overrides, **kwargs
        )
        builder.metadata.label = f'T_{temperature}_id_{i}'
        workchains.append(submit(builder))
        print(f'Launched <PwBaseWorkChain> with id={i} PK={workchains[-1].pk}')
        time.sleep(waiting_time)

    return workchains

def split_array(array: list, n: int) -> list[list]:
    """Split a generic array into N subarrays.
    
    .. note:: if `n` is larger then len(array)

    Args:
    ----
        array: a flat array to split into (semi)equal pieces.
        n: number of pieces

    """
    array = np.array(array)
    # Ensure N does not exceed the number of elements in the array
    n = min(n, len(array))
    
    # Calculate the size of each chunk
    chunk_sizes = np.full(n, len(array) // n)
    chunk_sizes[:len(array) % n] += 1

    # Generate the indices at which to split the array
    indices = np.cumsum(chunk_sizes)

    # Split the array at the calculated indices
    return np.split(array, indices[:-1])