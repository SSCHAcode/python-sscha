import numpy as np 
import sys, os

class Minimizer:
    def __init__(self, minim_struct = True, algorithm = "sdes", root_representation = False, step = 1):
        """
        This class is a minimizer that performs the optimization given the gradient of the dynamical matrix.

        Parameters
        ----------
            minim_struc : bool
                If true minimizes also the structure
            algorithm : string
                The algorithm used for the minimization
            root_representation : bool
                If true the minimization is performed in the root space
            step : double
                The initial step of the minimization
        """

        self.minim_struct = minim_struct
        self.algorithm = algorithm
        self.step = step
        self.setp_index = 0 # The index of the minimization

        self.old_direction = None 
        self.current_direction = None

        self.current_x = None 
        self.old_x = None
        self.old_kl = None 
        self.dyn = None

        self.new_direction = True


        # Some parameter to optimize the evolution

        # Maximum step of the kl ratio for the minimization
        self.kl_ratio_thr = 0.95
        # How much try to increase the step when a good step is found (must be > 1)
        self.increment_step = 2 
        # How much decrease the step when it is too big (must be < 1)
        self.decrement_step = 0.9

    def init(self, dyn, kl_ratio):
        """
        Initialize the minimizer with the current dynamical matrix and structure.

        Parameters
        ----------
            dyn : CC.Phonons.Phonons
                The dynamical matrix that identifies the starting point for the minimization.
            kl_ratio : float
                The initial Kong-Liu effective sample size
        """
        self.step_index = 0

        # create a vector of the dyn shape
        nq = len(dyn.q_tot)
        nmodes = dyn.structure.N_atoms * 3

        x_dyn = np.zeros( nq * nmodes**2, dtype = np.complex128)

        pos = 0
        for i in range(nq):
            x_dyn[pos : pos + nmodes**2] = dyn.dynmats[i].ravel()
        
        if self.minim_struct:
            x_struct = np.zeros(nmodes, dtype = np.complex128)
            x_struct[:] = dyn.structure.coords.ravel()
            self.current_x  = np.concatenate((x_dyn, x_struct))
        else:
            self.current_x = x_dyn
        
        self.old_x = self.current_x.copy()
        self.old_kl = kl_ratio
        self.new_direction = True

    def transform_gradients(self, dyn_gradient, structure_gradient = None):
        """
        Transform the gradients from dynamical matrix and structure
        to a single vector.
        The result of this function is what is needed by the run_step method.

        Parameters
        ----------
            dyn_gradient : ndarray( nq, nmodes, nmodes)
                The gradient of the dynamical matrix
            structure_gradient : ndarray(nmodes), optional
                The gradient of the structure (only needed if self.minim_struct = True)

        Results
        -------
            gradient : ndarray
                1D array where the gradients are collapsed all togheter.
        """

        if not self.minim_struct:
            return dyn_gradient.ravel()
        
        return np.concatenate( (dyn_gradient.ravel(), structure_gradient.ravel()) )

    def run_step(self, gradinet, kl_new):
        """
        Perform the minimization step with the line minimization
        """
        # Check consistency
        assert self.increment_step > 1 and self.decrement_step < 1

        if self.new_direction:
            # A new direction, update the position with the last one
            self.old_kl = kl_new
            if self.algorithm.lower() in ["sdes", "sd", "steepest descend"]:
                self.direction = gradient.copy()
            else:
                raise NotImplementedError("Error, algorithm '{}' not implemented".format(self.algorithm))
            
            self.old_x = self.current_x.copy()
            self.new_direction = False 

            # Enlarge the step
            self.step *= self.increment_step
        else:
            # Proceed with the line minimization

            # Compute the scalar product between the gradient and the direction
            # (Considering real and imaginary part as independent)
            scalar = np.dot(np.real(self.direction), np.real(gradient))
            scalar += np.dot(np.imag(self.direction), np.imag(gradient))

            kl_ratio = kl_new / self.old_kl

            # Check if the step was too big (scalar is negative or kl_ratio is below the threshold)
            # and decrement the step if needed
            if (scalar < 0) or (kl_ratio < self.kl_ratio_thr):
                self.step *= self.decrement_step
            else:
                # The step is good, therefore next step perform a new direction
                self.new_direction = True
        
        # Perform the minimiziation step
        self.current_x = self.old_x - self.step * self.direction


    def update_dyn(self, new_kl_ratio, dyn_gradient, structure_gradient = None):
        """
        Update the dynamical matrix.

        Parameters
        ----------
            new_kl_ratio : float
                Kong Liu effective sample size ratio of the current dynamical matrix
            dyn_gradient : ndarray( nq, nmodes, nmodes)
                The gradient of the dynamical matrix
            structure_gradient : ndarray(nmodes), optional
                The gradient of the structure (only needed if self.minim_struct = True)

        Results
        -------
            dyn : CC.Phonons.Phonons
                The updated dynamical matrix with the correct minimization step.
        """
        assert new_kl_ratio <= 1, "Error, the kl_ratio is defined between (0, 1], {} given.".format(new_kl_ratio)

        # TODO: Transform the dyn gradient in the root gradient

        grad_vector = self.transform_gradients(dyn_gradient, structure_gradient)
        self.run_step(grad_vector, new_kl_ratio)



        
        

