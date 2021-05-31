"""
Here the optimizer for the cell and the gradient
"""
from __future__ import print_function
import numpy as np
import SCHAModules
import sys, os


from sscha.Parallel import pprint as print

__EPSILON__ = 1e-8

class UC_OPTIMIZER:
    """
    This class is used as father of all the 
    other subclasses
    """
    
    # Avoid to reduce the step less than this value
    # Usefull, because in stochastic determination of the gradient it can
    # happen a negative bigger gradient that wuould bring alpha to 0.
    # NOTE: the maximum is always 2
    min_alpha_factor = 0.5
    
    # This variable set up how the minimum total decreasing alpha step
    # with respect to the initial one
    min_alpha_step = 1e-1
    
    alpha0 = 1
    
    # If true the strain is resetted each step
    # It should regolarize the minimization. 
    # If false the minimization will not go too far from the initial cell
    
    def __init__(self, starting_unit_cell):
        """
        To be initialized it needs to get the starting unit cell
        """
        
        self.last_grad = np.zeros(9, dtype = np.float64)
        self.last_direction = np.zeros(9, dtype = np.float64)
        self.alpha = np.float64(1)
        self.n_step = 0
        self.uc_0 = np.float64(starting_unit_cell.copy())
        self.uc_0_inv = np.linalg.inv(self.uc_0)
        self.x_start = None
        self.reset_strain = False

        self.uc_old = None

        # Allowed keys are sd or cg
        self.algorithm = "sd" 

        # If true the line minimization is employed during the minimization
        self.use_line_step = True
        
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
        good_step = True
        
        if self.n_step == 0:
            self.alpha0 = self.alpha
        
        if self.n_step != 0 and np.max(np.abs(self.last_direction)) > 1e-10:
            y0 = self.last_direction.dot(self.last_grad)
            y1 =  self.last_direction.dot(grad)

            print("[CELL] y0 = {} | y1 = {}".format(y0, y1))
            print("[CELL] grad = {}".format(grad))
            print("[CELL] lastgrad = {}".format(self.last_grad))
            sys.stdout.flush()
            factor = y0 / (y0 - y1)


            sys.stdout.flush()

            print("[CELL]  GRADIENT DOT DIRECTION = {}".format(factor))
            
            sys.stdout.flush()
            # Regularization (avoid run away)
            if factor > 2:
                factor = 2
                #good_step = False
            elif factor < 0 and y1 > 0:
                factor = 1.15
                # The minimization seem to have a negative curvature, continue with constant step
                print("[CELL] Warning: the cell minimization have a negative curvature.")
                print("                I do not know how to improve the minimization step.")
                print("                incrementing by a 15 %")
                
            elif factor < self.min_alpha_factor:
                factor = self.min_alpha_factor
                good_step = False

            if np.isnan(factor):
                factor = self.min_alpha_factor
                good_step = False
                
            self.alpha *= factor
        
            
            #if self.alpha < self.min_alpha_step * self.alpha0:
            #    self.alpha = self.min_alpha_step * self.alpha0
            #self.alpha *= factor
            
        return good_step

    def reset(self, unit_cell):
        """
        Reset the optimizer with the given unit cell
        (Usefull to reduce the numerical noise)
        """

        print("[CELL] resetting...")

        self.last_grad = np.zeros(9, dtype = np.double)
        self.last_direction = np.zeros(9, dtype = np.double)
        
        self.uc_0 = np.float64(unit_cell.copy())
        self.uc_0_inv = np.linalg.inv(self.uc_0)
        self.uc_old = self.uc_0.copy()
        self.x_start = np.zeros(9, dtype = np.double)
        self.reset_strain = False

    def check_reset(self, good_cell, thr = 10):
        """
        If the given good cell is much different from the original unit cell,
        reset the calculation.
        """

        if np.sqrt(np.sum( (self.uc_0 - good_cell)**2)) > thr:
            self.reset(good_cell)

        

    def perform_step(self, x_old, grad):
        """
        The step, hierarchical structure.
        Here a standard steepest descent
        """
        print()
        if self.x_start is None:
            self.x_start = x_old.copy()

        new_step = True
        if self.use_line_step:
            new_step = self.get_line_step(grad)


        if new_step:
            direction = self.get_new_direction(grad)
            x_new = x_old - direction * self.alpha
            self.x_start = x_old.copy()
            self.reset_strain = True
            print("[CELL] New step:")
            print("[CELL]    X_OLD = {}   | ALPHA = {}".format(x_old, self.alpha))
            print("[CELL]    DIRECTION = {}".format(direction))
        else:
            # Go back using the last direction selecting a smaller step
            x_new = self.x_start - self.last_direction * self.alpha
            print("[CELL] Step not good:")
            print("[CELL]    X_START = {}  | ALPHA = {}".format(self.x_start, self.alpha))
            print("[CELL]    DIRECTION = {}".format(self.last_direction))
            self.reset_strain = False

        self.n_step += 1
        print("[CELL]    GRADIENT = {}".format(grad))
        print("[CELL]    X_NEW = {}".format(x_new))
        print("[CELL]  Step number = {}".format(self.n_step))
        print()

        sys.stdout.flush()
        return x_new

    def get_new_direction(self, grad):
        # Steepest descent
        if self.algorithm == "sd":
            direction = grad
            self.last_direction = direction.copy()
            self.last_grad = grad.copy()
        elif self.algorithm == "cg":
            if np.max(np.abs(self.last_direction)) < 1e-10:
                # First step
                direction = grad
                self.last_direction = direction.copy()
                self.last_grad = grad.copy()
            else:
                gamma = np.dot(grad - self.last_grad, grad)
                gamma /= self.last_grad.dot(self.last_grad)

                direction = grad + gamma * self.last_direction
                if direction.dot(grad) < 0:
                    direction = grad
                    self.last_direction = grad.copy()
                    self.last_grad = grad.copy()
                    print("[CELL] Reset CG algorithm to SD")
                else:
                    self.last_direction = direction.copy()
                    self.last_grad = grad.copy()
        else:
            raise NotImplementedError("Error, cell algorithm setted to {}, but only 'sd' and 'cg' supported.".format(self.algorithm))
                    
            

        return direction

            
    def UpdateCell(self, unit_cell, stress_tensor, fix_volume = False, verbose = False):
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
            fix_volume : bool, optional
                If true, only the cell shape is affected by the update, the volume is taken as
                a constant
        """
        
        # Convert the stress tensor into the Free energy gradient with respect
        # to the unit cell
        volume = np.abs(np.linalg.det(unit_cell))
        #uc_inv = np.linalg.inv(unit_cell)
        I = np.eye(3, dtype = np.float64)

        uc_old = unit_cell.copy()
        if self.uc_old is None:
            self.uc_old = uc_old
            

        
        
        # Get the strain tensor up to know
        strain = np.transpose(self.uc_0_inv.dot(unit_cell) - I)
        
        if verbose:
            print ("ALPHA:", self.alpha)
            print ("VOLUME:", volume)
            
            print ("CURRENT STRAIN:")
            print (strain)
        
        # Get the gradient with respect to the strain
        grad_mat =  - volume * stress_tensor.dot(np.linalg.inv(I + strain.transpose()))
        
        # Modify the gradient if you need to fix the volume, in order to cancel the mean strain
        if fix_volume:
            grad_mat -= I * np.trace(grad_mat) / np.float64(3)
        
        #grad_mat = - volume * np.transpose(uc_inv).dot(stress_tensor)

        
        x_old = self.mat_to_line(strain)
        grad = self.mat_to_line(grad_mat)

        if verbose:
            print ("GRAD MAT:")
            print (grad_mat)
        
        x_new = self.perform_step(x_old, grad)
        
        strain_new  = self.line_to_mat(x_new)

        if verbose:
            print ("NEW STRAIN:")
            print (strain_new)
             
        unit_cell[:,:] = self.uc_0.dot( I + strain_new.transpose())
        
        if fix_volume:
            # Fix the volume
            unit_cell[:,:] *= (np.linalg.det(self.uc_0) / np.linalg.det(unit_cell))**(1/np.float64(3))
        
        if verbose:
            print ("NEW VOLUME:", np.linalg.det(unit_cell))
        


        # If the last step was good, check if the cell needs to be resetted.
        if self.reset_strain:
            self.check_reset(uc_old)
            self.reset_strain = False
        
            

class BFGS_UC(UC_OPTIMIZER):
    def __init__(self, unit_cell, bulk_modulus = 1, update_h_step = 0.2):
        """
        This is the BFGS algorithm adapted to
        optimize the unit cell.
        
        Parameters
        ----------
            unit_cell : array(3x3)
                The unit cell to initialize the minimizer
            bulk_modulus : float
                The static bulk modulus, initialize the step alpha
            update_h_step : float
                How much do you want to update the hessian matrix at each step?
                1 for standard BFGS
        """

        # Initialize the standard methods in the UC optimizer
        UC_OPTIMIZER.__init__(self, unit_cell)

        # BFGS estimates also the hessian
        volume = np.linalg.det(unit_cell)
        self.hessian = np.eye(9, dtype = np.float64)# * (3 *volume * bulk_modulus)
        self.alpha = 1 / (3 *volume * bulk_modulus)
        self.hessian_update = update_h_step
    
    def get_direction(self, grad):
        """
        Get the searching direction from the gradient.
        """
        
        p_vec =  - np.linalg.inv(self.hessian).dot(grad)
        return p_vec
    
            
    def get_hessian(self, grad):
        
        if self.n_step != 0:
            s_vec = self.alpha * self.last_direction
            y_vec = grad - self.last_grad
            
            U = np.outer(y_vec, y_vec) / y_vec.dot(s_vec)
            V1 = self.hessian.dot(s_vec)
            V = np.outer(V1, V1) / s_vec.dot(V1)
            
            self.hessian += self.hessian_update * (U - V)
            
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
        
        print ("HESSIAN eigvals:", 1 / np.linalg.eigvalsh(self.hessian))
        
        return x
            
class SD_PREC_UC(UC_OPTIMIZER):
    def __init__(self, unit_cell, bulk_modulus):
        """
        This is a quasi-Newton algorithm, where as a
        preconditioner is used the static bulk modulus
        """

        # Initialize the standard methods in the UC optimizer
        UC_OPTIMIZER.__init__(self, unit_cell)

        # Get the hessian matrix
        volume = np.abs(np.linalg.det(unit_cell))
        hessian = bulk_modulus * volume
    
        if np.shape(bulk_modulus)[0] != 9 or np.shape(bulk_modulus)[1] != 9:
            raise ValueError("Error, the bulk modulus must be a 9x9 matrix")

        # Correct adding 1 if the diagonal element are 0 (avoids signular matrix errors)
        for i in range(9):
            if np.abs(hessian[i,i]) < __EPSILON__:
                hessian[i,i] = 1

        self.preconditioner = np.linalg.inv(hessian)
        print ("PRECOND:")
        print (self.preconditioner)

    
    def perform_step(self, old_x, grad):
        """
        Perform the STEP
        """
        
        x = old_x - self.alpha * self.preconditioner.dot(grad)
        
        self.n_step += 1
        self.last_direction = - self.preconditioner.dot(grad)
        self.last_grad = grad
        return x
    
    
class CG_UC(UC_OPTIMIZER):
    """
    The optimizer that uses the conjugate gradient algoirthm.
    This converge faster in non isotropic bulk modulus.
    """
    
    def perform_step(self, old_x, grad):
        
        self.get_line_step(grad)
        
        if self.n_step >= 1:
            beta = grad.dot(grad) / self.last_grad.dot(self.last_grad)
            new_dir = grad + beta * self.last_direction
        else:
            new_dir = grad
        
        x = old_x - self.alpha * new_dir
        
        self.last_direction = new_dir
        self.last_grad = grad
        self.n_step += 1
        
        return x
        