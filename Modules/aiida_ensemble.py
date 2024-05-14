"""Module for handling automated calculation via aiida-quantumespress."""
from __future__ import annotations

import time
import numpy as np
from copy import copy
import cellconstructor

from .Ensemble import Ensemble
try:
    import aiida
    import aiida_quantumespresso
    from qe_tools import CONSTANTS
    
    gpa_to_rybohr3 = 1. / (CONSTANTS.ry_si/ CONSTANTS.bohr_si**3 / 1.0e9) # GPa -> Ry/Bohr^3
except ImportError:
    import warnings
    warnings.warn('aiida or aiida-quantumespresso are not installed')


class AiiDAEnsemble(Ensemble):
    """Ensemble subclass to interface SSCHA with aiida-quantumespresso."""

    def compute_ensemble(
        self, 
        pw_code: str,
        protocol: str = 'moderate',
        options: dict = None,
        overrides: dict = None,
        group_label: str = None,
        **kwargs
    ):
        """Get ensemble properties.

        All the parameters refer to the 
        :func:`aiida_quantumespresso.workflows.pw.base.PwBaseWorkChain.get_builder_from_protocol`
        method.

        Parameters
        ---------
            pw_code:
                The string associated with the AiiDA code for `pw.x`
            protocol:
                The protocol to be used; available protocols are 'fast', 'moderate' and 'precise'
            options:
                The options for the calculations, such as the resources, wall-time, etc.
            overrides:
                The overrides for the get_builder_from_protocol
            group_label:
                The group label where to add the submitted nodes for eventual future inspection
            kwargs:
                The kwargs for the get_builder_from_protocol
        """
        from aiida.orm import load_group
        
        group = None if group_label is None else load_group(group_label)        
        
        # Check if not all the calculation needs to be done
        if self.force_computed is None:
            self.force_computed = np.array([False] * self.N, dtype=bool)

        self.has_stress = True # by default we calculate stresses with the `get_builder_from_protocol`
        if overrides:
            try:
                tstress = overrides['pw']['parameters']['CONTROL']['tstress']
                self.has_stress = tstress
            except KeyError:
                pass
            
        # ============= AIIDA SECTION ============= #
        workchains = submit_and_get_workchains(
            structures=self.structures,
            pw_code=pw_code,
            temperature=self.current_T,
            protocol=protocol,
            options=options,
            overrides=overrides,
            **kwargs
        )
        
        if group:
            group.add_nodes(workchains)

        workchains_copy = copy(workchains)
        while(workchains_copy):
            workchains_copy = get_running_workchains(workchains_copy, self.force_computed)
            if workchains_copy:
                time.sleep(60) # wait before checking again
        
        for i, is_computed in enumerate(self.force_computed):
            if is_computed:
                out = workchains[i].outputs
                self.energies[i] = out.output_parameters.dict.energy / CONSTANTS.ry_to_ev
                self.forces[i] = out.output_trajectory.get_array('forces')[-1] / CONSTANTS.ry_to_ev
                if self.has_stress:
                    self.stresses[i] = out.output_trajectory.get_array('stress')[-1] * gpa_to_rybohr3
        # ============= AIIDA SECTION ============= #

        if self.has_stress:
            self.stress_computed = copy(self.force_computed)

        self._clean_runs()
        self.init()

    def _clean_runs(self) -> None:
        """Clean the failed runs and print summary."""
        n_calcs = np.sum(self.force_computed.astype(int))
        print('=============== SUMMARY AIIDA CALCULATIONS =============== \n')
        print('Total structures included: ', n_calcs)
        print('Structures not included  : ', self.N-n_calcs, '\n')
        print('===================== END OF SUMMARY ===================== \n')
        if n_calcs != self.N:
            self.remove_noncomputed()


def get_running_workchains(workchains: list, success: list[bool]) -> list:
    """Get the running workchains popping the finished ones.
    
    Two extra array should be given to populate the successfully finished runs.
    """
    wcs_left = copy(workchains)

    for workchain in workchains:
        if workchain.is_finished:
            if workchain.is_failed:
                print(f'[FAILURE] for <PwBaseWorkChain> with PK={workchain.pk}')
            else:
                index = int(workchain.label.split('_')[-1])
                success[index] = True
                print(f'[SUCCESS] for <PwBaseWorkChain> with PK={workchain.pk}')
                
            wcs_left.remove(workchain) # here it may be critical
    
    return wcs_left


def submit_and_get_workchains(
    structures: list[cellconstructor.Structure.Structure],
    pw_code: str,
    temperature: float | int,
    protocol: str = 'moderate',
    options: dict = None,
    overrides: dict = None,
    **kwargs
):
    """Submit and return the workchains for a list of :class:`~cellconstructor.Structure.Structure`.

    Parameters
    ---------
        structures:
            a list of :class:`~cellconstructor.Structure.Structure` to run via PwBaseWorkChain.
        pw_code:
            The string associated with the AiiDA code for `pw.x`
        temperature:
            The temperature corresponding to the structures ensemble
        protocol:
            The protocol to be used; available protocols are 'fast', 'moderate' and 'precise'
        options:
            The options for the calculations, such as the resources, wall-time, etc.
        overrides:
            The overrides for the get_builder_from_protocol
        kwargs:
            The kwargs for the get_builder_from_protocol
    """
    from aiida.engine import submit
    from aiida.plugins import WorkflowFactory
    from aiida.orm import StructureData

    PwBaseWorkChain = WorkflowFactory('quantumespresso.pw.base')

    structures_data = [StructureData(ase=cc.get_ase_atoms()) for cc in structures]
    workchains = []

    for i, structure in enumerate(structures_data):
        builder = PwBaseWorkChain.get_builder_from_protocol(
            code=pw_code,
            structure=structure,
            protocol=protocol,
            options=options,
            overrides=overrides,
            **kwargs
        )
        builder.metadata.label = 'T_{}_id_{}'.format(temperature, i)
        workchains.append(submit(builder))
        print(f'Launched <PwBaseWorkChain> with PK={workchains[-1].pk}')
    
    return workchains