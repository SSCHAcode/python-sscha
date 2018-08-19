# -*- coding: utf-8 -*-

"""
This file contains the SSCHA minimizer tool
It is possible to use it to perform the anharmonic minimization 
"""

#import Ensemble
import numpy as np
import matplotlib.pyplot as plt


# Rydberg to cm-1 and meV conversion factor
__RyToCm__  = 109691.40235
__RyTomev__ = 13605.698066


class SSCHA_Minimizer:
    
    def __init__(self, ensemble, root_representation = "normal",
                 kong_liu_ratio = 0.5, meaningful_factor = 0.1,
                 minimization_algorithm = "sdes"):
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
            kong_liu_ratio : float
                The ration of the Kong-Liu effective sample size below which
                the minimization is stopped and a new ensemble needs to be
                generated to proceed.
            meaningful_factor : float
                The ration between the gradient and its error below which
                the minimization is considered to be converged.
            minimization_algorithm : string
                The minimization algoirthm used. One between 'sdes', 'cgrf' or
                'auto'. They behave as follow:
                    - 'sdes' => Steepest Descent
                    - 'cgrf' => Congjugate gradient
                    - 'auto' => Uses cgrf if the error lower than the gradient, 'sdes' otherwise.
                
                NOTE: Only sdes is currently implemented.
        """
        
        self.ensemble = ensemble
        self.root_representation = root_representation
        
        
        # The symmetries
        self.symmetries = None
        
        # The minimization step
        self.min_step_dyn = 1e-4
        self.min_step_struc = 1
        
        self.dyn = self.ensemble.current_dyn.Copy()
        
        # Projection. This is chosen to fix some constraint on the minimization
        self.projector_dyn = None
        self.projector_struct = None
        
        # The gradient before the last step was performed (Used for the CG)
        self.prev_grad = None
        
        # Setup the statistical threshold
        self.kong_liu_ratio = kong_liu_ratio
        
        # Setup the meaningful_factor
        self.meaningful_factor = meaningful_factor
        
        # Setup the minimization algorithm
        self.minimization_algorithm = minimization_algorithm
        
        # Initialize the variable for convergence
        self.__converged__ = False
        
        # Initialize all the variables to store the minimization
        self.__fe__ = []
        self.__fe_err__ = []
        self.__gc__ = []
        self.__gc_err__ = []
        self.__gw__ = []
        self.__gw_err__ = []
        self.__KL__ = []
        
        
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
        dyn_grad, err = self.ensemble.get_free_energy_gradient_respect_to_dyn()
        
        # Store the gradient in the minimization
        self.__gc__.append(np.trace(dyn_grad.dot(dyn_grad)))
        self.__gc_err__.append(np.trace(err.dot(err)))
        
        # TODO: use the nonlinear change of variable to minimize correctly the dynamical matrix
        
        # Get the gradient of the free-energy respect to the structure
        struct_grad = - self.ensemble.get_average_forces()
        
        # Apply the translational symmetry on the structure
        # TODO: implement symmetries directly on the gradient of wyckoff
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
        self.dyn = new_dyn
        
        # Update the ensemble
        self.update()
        
        # Update the previous gradient
        self.prev_grad = dyn_grad
        
        
    def is_converged(self):
        """
        Simple method to check if the simulation is converged or
        requires a new population to be runned.
        
        Result
        ------
            bool : 
                True if the simulation ended for converging.
        """
        return self.__converged__
        
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
    
    def run(self, verbose = 1):
        """
        RUN THE SSCHA MINIMIZATION
        ==========================
        
        This function uses all the setted up parameters to run the minimization
        
        The minimization is stopped only when one of the stopping criteria are met.
        
        The verbose level can be chosen.
        
        Parameters
        ----------
            verbose : int
                The verbosity level.
                    - 0 : Noting is printed
                    - 1 : For each step only the free energy, the modulus of the gradient and 
                        the Kong-Liu effective sample size is printed.
                    - 2 : The dynamical matrix at each step is saved on output with a progressive integer
        """
        
        # Eliminate the convergence flag
        self.__converged__ = False
        
        # TODO: Activate a new pipe to avoid to stop the execution of the python 
        #       code when running the minimization. This allows for interactive plots
        running = True
        while running:
            # Compute the free energy and its error
            fe, err = self.get_free_energy(True)
            self.__fe__.append(fe)
            self.__fe_err__.append(err)
            
            
            # Compute the KL ratio
            self.__KL__.append(self.ensemble.get_effective_sample_size())
            
            # Perform the minimization step
            self.minimization_step(self.minimization_algorithm)
            
            
            # Print the step
            if verbose >= 1:
                print "Step ka = ", len(self.__fe__)
                print "Free energy = %16.8f +- %16.8f meV" % (self.__fe__[-1] * __RyTomev__, 
                                                              self.__fe_err__[-1] * __RyTomev__)
                print "FC gradient modulus = %16.8f +- %16.8f meV/A" % (self.__gc__[-1] * __RyTomev__, 
                                                                       self.__gc_err__[-1] * __RyTomev__)
                print "Kong-Liu effective sample size = ", self.__KL__[-1]
            
            if verbose >= 2:
                # Print the dynamical matrix at each step
                ka = len(self.__fe__)
                self.dyn.save_qe("minim_dyn_step%d_" % ka)
                
            # Get the stopping criteria
            running = not self.check_stop()
            
        
    def check_stop(self):
        """
        CHECK THE STOPPING CONDITION
        ============================
        
        Check the stopping criteria and returns True if the stopping
        condition is satisfied
        
        Result
        ------
            bool : 
                True if the minimization must be stopped, False otherwise
        
        """
        
        # Check the gradient
        last_gc = self.__gc__[-1]
        last_gc_err = self.__gc_err__[-1]
        
        if last_gc < last_gc_err * self.meaningful_factor:
            self.__converged__ = True
            return True
        
        # Check the KL
        kl = self.ensemble.get_effective_sample_size()
        
        if kl / float(self.ensemble.N) < self.kong_liu_ratio:
            self.__converged__ = False
            return True
        
        return False
            
    def plot_results(self):
        """
        PLOT RESULTS
        ============
        
        This usefull methods uses matplotlib to generate a plot of the
        minimization.
        
        """
        
        # Convert the data in numpy arrays
        fe = np.array(self.__fe__) * __RyTomev__
        fe_err = np.array(self.__fe_err__) * __RyTomev__
        
        gc = np.array(self.__gc__) * __RyTomev__
        gc_err = np.array(self.__gc_err__) * __RyTomev__
        
        kl = np.array(self.__KL__)
        
        steps = np.arange(len(fe))
        
        # Plot
        plt.figure()
        plt.title("Free energy")
        plt.errorbar(steps, fe, yerr = fe_err, label = "Free energy")
        plt.ylabel(r"$F$ [meV]")
        plt.xlabel("steps")
        
        plt.figure()
        plt.title("Gradient")
        plt.errorbar(steps, gc, yerr = gc_err, label = "gradient")
        plt.ylabel(r"$|\vec g|$ [meV / A]")
        plt.xlabel("steps")
        
        plt.figure()
        plt.title("Kong-Liu effective sample size")
        plt.plot(steps, kl)
        plt.ylabel(r"$\frac{N_{eff}}{N_0}$")
        plt.xlabel("steps")
        
        plt.show()
            
    

def get_root_dyn(dyn_fc, root_representation):
    """
    Get the root dyn matrix
    
    
    This method computes the root equivalent of the dynamical matrix
    """
    # TODO: To be ultimated
    pass
