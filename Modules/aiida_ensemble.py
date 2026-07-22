# -*- coding: utf-8 -*-
"""Module for handling automated calculation via aiida-quantumespresso."""
from __future__ import annotations

from typing import Literal
from copy import copy
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
                self._otf_setup_defaults(number_of_atoms)
                self._otf_predict(structures, dft_indices)
                
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
                    self._otf_maybe_train_and_write(dft_counts)
                
                sys.stdout.flush()

        # ================ FINALIZE ================ #
        # if self.has_stress:
        #     self.stress_computed = copy(self.force_computed)

        self._clean_runs(dft_counts, label="AIIDA CALCULATIONS")
        self.init()


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