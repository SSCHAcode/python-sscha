# -*- coding: utf-8 -*-

"""
This file contains the SSCHA minimizer tool
It is possible to use it to perform the anharmonic minimization 
"""

#import Ensemble
import numpy as np

class SSCHA_Minimizer:
    
    def __init__(self, ensemble, root_representation = "normal"):
        """
        This class create a minimizer to perform the sscha minimization.
        It performs the sscha minimization.
        
        Parameters
        ----------
            ensemble : Ensemble.Ensemble()
                This is the Ensemble. This class contains the pool of configurations to be
                used in the minimizations.
            root_representation : string
                Chose between "normal", "sqrt" and "root4". These are the nonlinear change
                of variable to speedup the code.
        """
        
        self.ensemble = ensemble
        self.root_representation = root_representation
        
        
        # The symmetries
        self.symmetries = None
        
        # The minimization step
        self.min_step_dyn = 1
        self.min_step_struc = 1
        
        self.dyn = self.ensemble.current_dyn.Copy()
        
        # Projection. This is chosen to fix some constraint on the minimization
        self.projector_dyn = None
        self.projector_struct = None
        
        # The gradient before the last step was performed (Used for the CG)
        self.prev_grad = None
        
        
    def minimization_step(self, algorithm = "sdes"):
        """
        Perform the single minimization step.
        This modify the self.dyn matrix and updates the ensemble
    
        
        Parameters
        ----------
            algorithm : str
                The minimization algorithm. By default it is steepest descent.
                Supported altorithms:
                    - "sdes"
        """
        
        if algorithm != "sdes":
            raise ValueError("Error, %s algorithm is not supported." % algorithm)
        
        # Get the gradient of the free-energy respect to the dynamical matrix
        dyn_grad = self.ensemble.get_free_energy_gradient_respect_to_dyn()
        
        # TODO: use the nonlinear change of variable to minimize correctly the dynamical matrix
        
        # Get the gradient of the free-energy respect to the structure
        struct_grad = - self.ensemble.get_average_forces()
        
        # Apply the translational symmetry on the structure
        struct_grad -= np.tile(np.einsum("ij->j", struct_grad), (self.ensemble.current_dyn.structure.N_atoms,1))
        
        
        current_dyn = self.ensemble.current_dyn
        current_struct = self.ensemble.current_dyn.structure
        
        # Perform the step for the dynamical matrix
        new_dyn = current_dyn
        new_dyn.dynmats[0] = current_dyn.dynmats[0] - self.min_step_dyn * dyn_grad
        
        # Perform the step for the structure
        current_struct.coords -= struct_grad * self.min_step_dyn
        
        # Symmetrize the structure after the gradient is applied
        if self.symmetries is not None:
            current_struct.impose_symmetries(self.symmetries, verbose = False)
        
        new_dyn.structure = current_struct
        
        # Update the ensemble
        self.ensemble.update_weights(new_dyn, self.ensemble.current_T)
        self.dyn = new_dyn
        
        # Update the previous gradient
        self.prev_grad = dyn_grad
        
    def update(self):
        """
        UPDATE IMPORTANCE SAMPLING
        ==========================
        
        
        This methods makes the self.dyn coincide with self.ensemble.current_dyn, and overwrites the stochastic
        weights of the current_dyn.
        
        Call this method each time you modify the dynamical matrix of the minimization to avoid errors.
        
        NOTE: it is equivalent to call self.ensemble.update_weights(self.dyn, self.ensemble.current_T)
        """
        
        self.ensemble.update_weights(self.dyn, self.ensemble.current_T)
        
        
    def get_free_energy(self, return_error = False):
        """
        SSCHA FREE ENERGY
        =================
        
        Obtain the SSCHA free energy for the system.
        This is done by integrating the free energy along the hamiltonians, starting
        from current_dyn to the real system.
        
        The result is in Rydberg.
        
        NOTE: this method just recall the self.ensemble.get_free_energy function.
        
        .. math::
            
            \\mathcal F = \\mathcal F_0 + \\int_0^1 \\frac{d\\mathcal F_\\lambda}{d\\lambda} d\\lambda
        
        Where :math:`\\lambda` is the parameter for the adiabatic integration of the hamiltonian.
        
        .. math::
            
            H(\\lambda) = H_0 + (H - H_0) \\lambda
        
        here :math:`H_0` is the sscha harmonic hamiltonian, while :math:`H_1` is the real hamiltonian 
        of the system.
        
        Returns
        -------
            float
                The free energy in the current dynamical matrix and at the ensemble temperature
        
        """
        
        #TODO: CHECK THE CONSISTENCY BETWEEN THE DYNAMICAL MATRICES
        # Check if the dynamical matrix has correctly been updated
        #if np.sum( self.dyn != self.ensemble.current_dyn):
        #    raise ValueError("Error, the ensemble dynamical matrix has not been updated. You forgot to call self.update() before")
        
        return self.ensemble.get_free_energy(return_error = return_error)

def get_root_dyn(dyn_fc, root_representation):
    """
    Get the root dyn matrix
    
    
    This method computes the root equivalent of the dynamical matrix
    """
    # TODO: To be ultimated
    pass
