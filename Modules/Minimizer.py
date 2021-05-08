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
        self.nq = 0
        self.n_modes = 0

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
        self.nq = len(dyn.q_tot)
        self.n_modes = dyn.structure.N_atoms * 3

        x_dyn = np.zeros( self.nq * self.n_modes**2, dtype = np.complex128)

        pos = 0
        for i in range(nq):
            x_dyn[pos : pos + self.n_modes**2] = dyn.dynmats[i].ravel()

        # Transform in root space if needed.
        x_dyn, = get_root_dyn_grad(x_dyn, np.zeros(x_dyn.shape, dtype = np.complex128), self.root_representation)
        
        if self.minim_struct:
            x_struct = np.zeros(self.n_modes, dtype = np.complex128)
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

    def get_dyn_struct(self):
        """
        From the current position of the minimization,
        get back the dynamical matrix and the structure
        """

        dynq = np.zeros( (self.nq, self.n_modes, self.n_modes), dtype = np.complex128)

        for i in range(self.nq):
            dynq[i,  :, :] = self.current_x[ i * self.n_modes**2 : (i+1) * self.n_modes**2].reshape( (self.n_modes, self.n_modes))
        
        # Revert the dynamical matrix if we are in the root representation
        dynq = get_standard_dyn(dynq, self.root_representation)

        struct = self.dyn.structure.coords.copy()
        if self.minim_struct:
            struct = self.current_x[self.nq * self.n_modes * self.n_modes :].reshape((self.n_modes, 3))
        return dynq, struct

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

        current_dyn, = self.get_dyn_struct()

        # Now we can obtain the gradient in the root representation
        root_dyn, root_grad = get_root_dyn_grad(current_dyn, dyn_gradient, self.root_representation)

        grad_vector = self.transform_gradients(root_grad, structure_gradient)
        self.run_step(grad_vector, new_kl_ratio)

        



def get_root_dyn_grad(dyn_q, grad_q, root_representation = "sqrt"):
    """
    ROOT MINIMIZATION STEP
    ======================
    
    As for the [Monacelli, Errea, Calandra, Mauri, PRB 2017], the nonlinear
    change of variable is used to perform the step.
    
    It works as follows:
    
    .. math::
        
        \\Phi \\rightarrow \\sqrt{x}{\\Phi}
        
        \\frac{\\partial F}{\\partial \\Phi} \\rightarrow \\frac{\\partial F}{\\partial \\sqrt{x}{\\Phi}}
        
        \\sqrt{x}{\\Phi^{(n)}} \\stackrel{\\frac{\\partial F}{\\partial \\sqrt{x}{\\Phi}}}{\\longrightarrow} \\sqrt{x}{\\Phi^{(n+1)}}
        
        \\Phi^{(n+1)} = \\left(\\sqrt{x}{\\Phi^{(n+1)}})^{x}
        
    Where the specific update step is determined by the minimization_algorithm, while the :math:`x` order 
    of the root representation is determined by the root_representation argument.
    
    Parameters
    ----------
        dyn_q : ndarray( NQ x 3nat x 3nat )
            The dynamical matrix in q space. The Nq are the total number of q.
        grad_q : ndarray( NQ x 3nat x 3nat )
            The gradient of the dynamical matrix.
        step_size : float
            The step size for the minimization
        root_representation : string
            choice between "normal", "sqrt" and "root4". The value of :math:`x` will be, respectively,
            1, 2, 4.
        minimization_algorithm : string
            The minimization algorithm to be used for the update.
    
    Result
    ------
        new_dyn_q : ndarray( NQ x 3nat x 3nat )
            The updated dynamical matrix in q space
    """
    ALLOWED_REPRESENTATION = ["normal", "root4", "sqrt", "root2"]

    # Avoid case sensitiveness
    root_representation = root_representation.lower()
    
    if not root_representation in ALLOWED_REPRESENTATION:
        raise ValueError("Error, root_representation is %s must be one of '%s'." % (root_representation,
                                                                                  ", ".join(ALLOWED_REPRESENTATION)))
    # Allow also for the usage of root2 instead of sqrt
    if root_representation == "root2":
        root_representation = "sqrt"
    
    # Apply the diagonalization
    nq = np.shape(dyn_q)[0]
    
    # Check if gradient and dyn_q have the same number of q points
    if nq != np.shape(grad_q)[0]:
        raise ValueError("Error, the number of q point of the dynamical matrix %d and gradient %d does not match!" % (nq, np.shape(grad_q)[0]))
    
    # Create the root representation
    new_dyn = np.zeros(np.shape(dyn_q), dtype = np.complex128)
    new_grad = np.zeros(np.shape(dyn_q), dtype = np.complex128)

    # Copy
    new_dyn[:,:,:] = dyn_q
    new_grad[:,:,:] = grad_q

    # Cycle over all the q points
    if root_representation != "normal":
        for iq in range(nq):
            # Dyagonalize the matrix
            eigvals, eigvects = np.linalg.eigh(dyn_q[iq, :, :])
            
            # Regularize acustic modes
            if iq == 0:
                eigvals[eigvals < 0] = 0.
            
            # The sqrt conversion
            new_dyn[iq, :, :] = np.einsum("a, ba, ca", np.sqrt(eigvals), eigvects, np.conj(eigvects))
            new_grad[iq, :, :] = new_dyn[iq, :, :].dot(grad_q[iq, :, :]) + grad_q[iq, :, :].dot(new_dyn[iq, :, :])
            
            # If root4 another loop needed
            if root_representation == "root4":                
                new_dyn[iq, :, :] = np.einsum("a, ba, ca", np.sqrt(np.sqrt(eigvals)), eigvects, np.conj(eigvects))
                new_grad[iq, :, :] = new_dyn[iq, :, :].dot(new_grad[iq, :, :]) + new_grad[iq, :, :].dot(new_dyn[iq, :, :])

    return new_dyn, new_grad

def get_standard_dyn(root_dyn, root_representation):
    """
    Get the standard dynamical matrix from the root representation

    Parameters
    ----------
        root_dyn : ndarray (nq x 3nats x 3nats)
            The dynamical matrix in the root representation
        root_representation : string
            The kind of root representation
    """

    new_dyn = root_dyn.copy()
    if root_representation != "normal":
        for iq in range(nq):
            # The square root conversion
            new_dyn[iq, :, :] = new_dyn[iq, :,:].dot(new_dyn[iq, :,:])
            
            # The root4 conversion
            if root_representation == "root4":
                new_dyn[iq, :, :] = new_dyn[iq, :,:].dot(new_dyn[iq, :,:])
    
    # Return the matrix
    return new_dyn
            
            
