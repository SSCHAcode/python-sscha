"""
Here the optimizer for the cell and the gradient
"""
import numpy as np
import SCHAModules

class UC_OPTIMIZER:
    """
    This class is used as father of all the 
    other subclasses
    """
    def __init__(self):
        
        self.last_grad = np.zeros(9, dtype = np.float64)
        self.last_direction = np.zeros(9, dtype = np.float64)
        self.alpha = 1
        self.n_step = 0

        
    def mat_to_line(self, matrix):
        """
        From a 3x3 unit cell (or stress) matrix
        obtain the 6 degrees of freedom.
        """
        # line = np.zeros(6, dtype= np.float64)
        # line[:3] = matrix[0, :]
        # line[3:5] = matrix[1, 1:]
        # line[5] = matrix[2,2]
        return matrix.ravel(order = "C")
    
    def line_to_mat(self, line):
        """
        from a line 6 element array to a 3x3 matrix
        """
        # matrix = np.zeros( (3,3), dtype = np.float64)
        # matrix[0,:] = line[:3]
        # matrix[:,0] = line[:3]
        # matrix[1, 1:] = line[3:5]
        # matrix[1:, 1] = line[3:5]
        # matrix[2,2] = line[5]
        return line.reshape((3,3), order = "C")

    def get_line_step(self, grad):
        """
        LINE MINIMIZATION
        =================

        This function is in common for all the minimizer
        it implements the line minimization.
        """
        
        if self.n_step != 0:
            y0 = self.last_direction.dot(self.last_grad)
            y1 =  self.last_direction.dot(grad)
            factor = y0 / (y0 - y1)
            
            # Regularization (avoid run away)
            #self.alpha *= 1 + np.tanh( (factor - 1) )
            self.alpha *= factor

    def perform_step(self, x_old, grad):
        """
        The step, hierarchical structure.
        Here a standard steepest descent
        """
        #self.get_line_step(grad)
        self.last_direction = grad
        self.last_grad = grad

        x_new = x_old - grad * self.alpha
        self.n_step += 1
        return x_new

            
    def UpdateCell(self, unit_cell, stress_tensor):
        """
        PERFORM THE CELL UPDATE
        =======================
        
        This method perform the unit cell update using
        the BFGS algorithm.
        
        Parameters
        ----------
            unit_cell : 3x3 matrix
                The unit cell of the system. It will be
                updated after the minimization step.
            stress_tensor : 3x3 matrix
                The stress tensor according to which you want
                to optimize the gradient.
        """
        
        # Convert the stress tensor into the Free energy gradient with respect
        # to the unit cell
        volume = np.linalg.det(unit_cell)
        uc_inv = np.linalg.inv(unit_cell)
        grad_mat = - volume * np.transpose(uc_inv).dot(stress_tensor)

        # Other grad mat
        new_grad = SCHAModules.cell_force(np.transpose(uc_inv), stress_tensor, volume, 0, 1)
        new_grad = np.transpose(new_grad)

        print "NEW GRADIENT FORTRAN:"
        print new_grad
        
        
        x_old = self.mat_to_line(unit_cell)
        grad = self.mat_to_line(grad_mat)

        print "GRAD MAT:"
        print grad_mat
        
        x_new = self.perform_step(x_old, grad)
        
        unit_cell[:,:] = self.line_to_mat(x_new)


class BFGS_UC(UC_OPTIMIZER):
    def __init__(self):
        """
        This is the BFGS algorithm adapted to
        optimize the unit cell.
        """

        # Initialize the standard methods in the UC optimizer
        UC_OPTIMIZER.__init__(self)

        # BFGS estimates also the hessian
        self.hessian = np.eye(9, dtype = np.float64)
    
    def get_direction(self, grad):
        """
        Get the searching direction from the gradient.
        """
        
        p_vec =  - np.linalg.inv(self.hessian).dot(grad)
        return p_vec
    
    def setup_hessian_from_bulk_modulus(self, unit_cell, static_bm):
        """
        Setup the hessian matrix from a consant uniform bulk modulus
        """
        # TODO: To be implemented
        pass
            
    def get_hessian(self, grad):
        
        if self.n_step != 0:
            s_vec = self.alpha * self.last_direction
            y_vec = grad - self.last_grad
            
            U = np.outer(y_vec, y_vec) / y_vec.dot(s_vec)
            V1 = self.hessian.dot(s_vec)
            V = np.outer(V1, V1) / s_vec.dot(V1)
            
            self.hessian += U - V
            
    def perform_step(self, old_x, grad):
        """
        Perform the BFGS STEP
        """
        self.get_hessian(grad)
        self.get_line_step(grad)
        new_dir = self.get_direction(grad)
        
        x = old_x + self.alpha * new_dir
        self.n_step += 1
        self.last_direction = new_dir
        self.last_grad = grad
        
        return x
            
