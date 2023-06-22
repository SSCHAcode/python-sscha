"""Module for handling automated calculation via aiida-quantumespress."""
from __future__ import annotations

import time
import numpy as np

from .Ensemble import Ensemble
try:
    import aiida
    import aiida_quantumespresso
except ImportError:
    import warnings
    warnings.warn('aiida or aiida-quantumespresso are not installed')


class AiiDAEnsemble(Ensemble):
    """Ensemble subclass to interface SSCHA with aiida-quantumespresso."""

    def compute_ensemble(
        self, 
        pw_code: str,
        profile: str | None = None, 
        protocol: str = 'moderate',
        options: dict = None,
        overrides: dict = None,
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
            profile: 
                AiiDA profile to use during the calculations. If None, the default is used
            protocol:
                The protocol to be used; available protocols are 'fast', 'moderate' and 'precise'
            options:
                The options for the calculations, such as the resources, wall-time, etc.
            overrides:
                The overrides for the get_builder_from_protocol
            kwargs:
                The kwargs for the get_builder_from_protocol
        """
        from aiida import load_profile
        from aiida.engine import submit
        from aiida.plugins import WorkflowFactory
        from aiida.orm import StructureData
        
        PwBaseWorkChain = WorkflowFactory('quantumespresso.pw')
    
        if profile:
            load_profile(profile)
        else:
            load_profile()
        
        # Check if not all the calculation needs to be done
        n_calcs = np.sum( self.force_computed.astype(int)) 
        computing_ensemble = self

        if overrides:
            try:
                tstress = overrides['pw']['parameters']['CONTROL']['tstress']
                self.has_stress = tstress
            except KeyError:
                self.has_stress = True # by default we calculate stresses with the `get_builder_from_protocol`

        # Check wheter compute the whole ensemble, or just a small part
        should_i_merge = False
        if n_calcs != self.N:
            should_i_merge = True
            computing_ensemble = self.get_noncomputed()
            self.remove_noncomputed()
            
        # ============= HERE AIIDA PART ============= #
        structures_data = [StructureData(ase=cc.get_ase_atoms()) for cc in self.structures]
        workchains = []
        
        for i, structure in enumerate(structures_data):
            builder = PwBaseWorkChain.get_builder_from_protocol(
                code=pw_code,
                structure=structure,
                options=options,
                overrides=overrides,
                **kwargs
            )
            builder.label = str(i)
            workchains.append(submit(builder))
        
        workchains_success = []
        workchains_failed = []

        while(workchains):
            for workchain in workchains:
                if workchain.is_finished:
                    if workchain.is_failed:
                        workchains_failed.append(workchain)
                    else:
                        workchains_success.append(workchain)
                    workchains.remove(workchain) # here it may be critical

            time.sleep(60) # wait before checking again
        
        # ---- TO DO
        #   * Take energies, forces, stresses
        #   * Print which one was successfull
        #   * Add correct labels to workchains to gather afterwards
        #   * For sure something else
        #   * Add the possibility of adding wc to group
        #   * Keep track is some better way (printing?) of the workchains submitted
        # ============= HERE AIIDA PART ============= #

        print("CE BEFORE MERGE:", len(self.force_computed))

        if should_i_merge:
            # Remove the noncomputed ensemble from here, and merge
            self.merge(computing_ensemble)
        print("CE AFTER MERGE:", len(self.force_computed))

        print('ENSEMBLE ALL PROPERTIES:', self.all_properties)

