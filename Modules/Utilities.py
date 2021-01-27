"""
This module contains the prefined custom function to
perform some standard variants of the minimization.

The function are classified as:
CGF : custom_gradient_function

"""
from __future__ import print_function
import difflib
import cellconstructor as CC
import cellconstructor.Phonons
import numpy as np

import pickle

__UTILS_NAMESPACE__ = "utils"
__UTILS_SAVEFREQ_FILENAME__ = "save_freq_filename"
__UTILS_SAVERHO_FILENAME__ = "save_rho_filename"
__UTILS_LOCKMODE_START__ = "mu_lock_start"
__UTILS_LOCKMODE_END__ = "mu_lock_end"
__UTILS_FREEMODE_START__ = "mu_free_start"
__UTILS_FREEMODE_END__ = "mu_free_end"
__UTILS_PROJECT_DYN__ = "project_dyn"
__UTILS_PROJECT_STRUCTURE__ = "project_structure"

__ALLOWED_KEYS__ = [__UTILS_SAVEFREQ_FILENAME__, __UTILS_SAVERHO_FILENAME__,
                    __UTILS_LOCKMODE_START__, __UTILS_LOCKMODE_END__,
                    __UTILS_FREEMODE_START__, __UTILS_FREEMODE_END__,
                    __UTILS_PROJECT_DYN__, __UTILS_PROJECT_STRUCTURE__]


NOT_INITIALIZED_ERROR = """
Error, the ModeLocking has not been initialized, 
please, consider running a function like: SetupFreeModes
before passing the mode locking function to the minimizer.
"""


def get_custom_functions_from_namelist(namelist, dyn):
    """
    GET CFUNCTIONS FROM NAMELIST
    ============================
    
    This method reads a namelist and returns the three pointers
    to thee custom functions.
    
    NOTE: Up to know the mode locking is not abilitated
    
    Parameters
    ----------
        namelist : dict
            The dictionary containing the parsed namelist
        dyn : CC.Phonons.Phonons
            The dynamical matrix, used if you want to setup the
            mode locking
    
    Results
    -------
        cf_pre, cf_grad, cf_post:
            Three custom function to be executed, respectively,
            before, during and after the minimization step.
    """
    if not __UTILS_NAMESPACE__ in namelist.keys():
        raise IOError("Error, the namespace %s is missing" % __UTILS_NAMESPACE__)
    
    c_info = namelist[__UTILS_NAMESPACE__]
    keys = c_info.keys()
    
    use_io = False
    use_modelocking = False
    
    io_info = IOInfo()
    
    # Check for unknown keys
    for k in keys:
        if not k in __ALLOWED_KEYS__:
            print ("Error with the key:", k)
            s =  "Did you mean something like:" + str( difflib.get_close_matches(k, __ALLOWED_KEYS__))
            print (s)
            raise IOError("Error in "+__UTILS_NAMESPACE__+" namespace: key '" + k +"' not recognized.\n" + s)
    
    
    
    if __UTILS_SAVEFREQ_FILENAME__ in keys:
        use_io = True
        io_info.SetupSaving(c_info[__UTILS_SAVEFREQ_FILENAME__])
    
    if __UTILS_SAVERHO_FILENAME__ in keys:
        use_io = True
        io_info.SetupWeights(c_info[__UTILS_SAVERHO_FILENAME__])
    
    def cfp(minim):
        if use_io:
            print ("SAVING")
            return io_info.CFP_SaveAll(minim)
        else:
            print ("NOT SAVING")
        
    # Setup the mode projection
    locking = False
    mu_start = 0
    mu_end = dyn.structure.N_atoms * 3
    if __UTILS_LOCKMODE_START__ in keys:
        use_modelocking = True
        locking = True
        mu_start = int(c_info[__UTILS_LOCKMODE_START__]) - 1
    if __UTILS_LOCKMODE_END__ in keys:
        use_modelocking = True
        locking = True
        mu_end = int(c_info[__UTILS_LOCKMODE_END__]) - 1
    if __UTILS_FREEMODE_START__ in keys:
        use_modelocking = True
        if locking:
            raise ValueError("Error, you cannot set both to free and lock modes.")
        
        mu_start = int(c_info[__UTILS_FREEMODE_START__]) - 1
    if __UTILS_FREEMODE_END__ in keys:
        use_modelocking = True
        if locking:
            raise ValueError("Error, you cannot set both to free and lock modes.")
        
        mu_end = int(c_info[__UTILS_FREEMODE_END__]) - 1
        
    project_structure = True
    project_dyn = True
    if __UTILS_PROJECT_DYN__ in keys:
        project_dyn = bool(c_info[__UTILS_PROJECT_DYN__])
    if __UTILS_PROJECT_STRUCTURE__ in keys:
        project_structure = bool(c_info[__UTILS_PROJECT_STRUCTURE__])
    
    
    # Get the pols
    nat = dyn.structure.N_atoms
    nq = len(dyn.q_tot)
    n_selected = mu_end - mu_start
    if n_selected < 0:
        raise ValueError("Error, the start mode cannot be smaller than the final one.")
        
    pols = np.zeros( (3*nat, n_selected, nq), dtype = np.complex128, order = "F")
    
    if use_modelocking:
        
        if mu_start < 0 or mu_start >= 3*nat:
            raise ValueError("Error, the modes specified for modelocking must be between %d and %d" % (1, 3*nat))
        
        if mu_end < 0 or mu_end >= 3*nat:
            raise ValueError("Error, the modes specified for modelocking must be between %d and %d" % (1, 3*nat))

        # for iq in range(nq):
        #     w, p_v = dyn.DyagDinQ(iq)
        #     pols[:, :, iq] = p_v[:, mu_start : mu_end]
        
        # # Setup the mode projection
        # ModProj = ModeProjection(pols, dyn.structure.get_masses_array())
        ModProj = ModeProjection(dyn)
        ModProj.SetupFreeModes(mu_start, mu_end)
    
        def modlock(dyn_grad, struct_grad):
            if project_dyn and project_structure:
                if locking:
                     # Project out the modes
                     dyn1 = dyn_grad.copy()
                     struc1 = struct_grad.copy()
                     
                     ModProj.CFG_ProjectOnModes(dyn1, struc1)
                     dyn_grad -= dyn1
                     struct_grad -= struc1
                else:
                     ModProj.CFG_ProjectOnModes(dyn_grad, struct_grad)
            elif project_dyn:
                if locking:
                     # Project out the modes
                     dyn1 = dyn_grad.copy()
                     
                     ModProj.CFG_ProjectDyn(dyn1, None)
                     dyn_grad -= dyn1
                else:
                     ModProj.CFG_ProjectDyn(dyn_grad, None)
            elif project_structure:
                if locking:
                     # Project out the modes
                     struc1 = struct_grad.copy()
                     
                     ModProj.CFG_ProjectStructure(None, struc1)
                     struct_grad -= struc1
                else:
                     ModProj.CFG_ProjectStructure(None, struct_grad)
            else:
                raise ValueError("Internal error, asked for mode locking but neither the dyn or the structure is locked.")
        
        return None, modlock, cfp
    
    return None, None, cfp
        


class ModeProjection:
    def __init__(self, dyn):#pols, masses_array):
        """
        This class is used to project the gradient in or out of
        a space generated by the mode in the pols array.
        
        The projection on the subspace is done in the following way.
        The force constant matrix is transformed in the dynamical matrix:
            
        .. math::
            
            D = M^{-\\frac 1 2} \\Phi M^{-\\frac 1 2}
            
            D_{proj} = \\sum_{\\mu} \\left|e_\\mu\\right> \\left<e_\\mu\\right| D \\left|e_\\mu\\right> \\left<e_\\mu\\right|
            
            \\Phi_{proj} = M^{\\frac 1 2} D_{proj} M^{\\frac 1 2}
            
        This means that the projection is done with the two operators
        
        .. math::
            
            \\Phi_{proj} = P^\\dagger \\Phi P
            
            P = \\sum_{\\mu} M^{-\\frac 12}\\left|e_\\mu\\right> \\left<e_\\mu\\right| M^{\\frac 12}
            
        The same matrix is used to project the vectors.

        Parameters
        ----------
            dyn : Phonons()
                The dynamical matrix
        """
        
        # Prepare the modes to be locked
        self.nat = dyn.structure.N_atoms 
        self.nmodes = 3*self.nat 
        self.nq = len(dyn.q_tot)

        # Get the polarization vectors for each q point
        self.pols = np.zeros( (3*self.nat, self.nmodes, self.nq), dtype = np.complex128, order = "F")
        for iq in range(self.nq):
            w, p = dyn.DyagDinQ(iq)
            self.pols[:, :, iq] = p 

            # Check normalization of p
            assert np.max(np.abs(np.conj(p.T).dot(p) - np.eye(self.nmodes))) < 1e-6, "Error, for some reason vectors are not normalized"

        self.masses = dyn.structure.get_masses_array()
        
        
        # Prepare the projector
        self.projector = np.zeros( (self.nq, 3*self.nat, 3*self.nat), dtype = np.complex128)
        self.projectorH = np.zeros( (self.nq, 3*self.nat, 3*self.nat), dtype = np.complex128)
        self.proj_vec = np.zeros( (3*self.nat, 3*self.nat), dtype = np.float64)

        self.mu_start = 0#index_mode_start 
        self.mu_end = 0#index_mode_end


        self.testing = False

        self.initialized = False

    def SetupFreeModes(self, index_mode_start, index_mode_end):
        """
        SETUP FREE MODES
        ================

        Constrain the minimization only in the modes between index_modes_start and
        index_mdoe_end in all the q points.

        For each q points only modes between these indices will be minimized.
        """
        self.initialized = True

        self.mu_start = index_mode_start 
        self.mu_end = index_mode_end

        # Generate the array for the masses aligned as the polarization vector
        _m_ = np.tile(self.masses, (3, 1)).ravel(order = "F")
        _msq_ = np.sqrt(_m_)
        
        # Setup the projector on the dynamical matrix
        for iq in range(self.nq):
            for mu in range(index_mode_start, index_mode_end):
                pvec = self.pols[:, mu, iq].copy()
                
                #pvec /= pvec.dot(pvec) # Normalization
                self.projector[iq, :, :] += np.outer(pvec / _msq_, np.conj(pvec) * _msq_ )
                self.projectorH[iq, :, :] += np.outer(pvec*_msq_, np.conj(pvec) / _msq_)
                
        # Prepare the projector on the structure
        for mu in range(index_mode_start, index_mode_end):
            pvec = np.real(self.pols[:, mu, 0])
            #pvec /= pvec.dot(pvec)
            self.proj_vec[:,:] += np.outer(pvec, pvec)
                

        # Impose the sum rule 
        # Note that in principle it should be satisfied, however the python diagonalization
        # sucks, therefore the polarization vector are not exactly orthonormal.
        CC.symmetries.CustomASR(self.projector[0, :, :])

    def CFG_ProjectOnModes(self, dyn_grad, struct_grad):
        """
        PROJECT THE GRADIENT IN THE SELECTED MODES
        ==========================================


        Function to be passed to the minimizer as 'custom_function_gradient'. 
        It project the gradients in the polarization vector subspace.
        As any custom_function_gradient, it takes as input the two gradients.
        """

        assert self.initialized, NOT_INITIALIZED_ERROR

        # Project the structure in the polarization vectors
        struct_grad_new = self.proj_vec.dot(struct_grad.ravel())
        struct_grad = struct_grad_new.reshape((self.nat, 3))

        _m_ = np.tile(self.masses, (3, 1)).T.ravel()

        # Do the same for the matrix
        for iq in range(self.nq):

            # remove the masses from the gradient
            grad_nomass = dyn_grad[iq, :, :] / np.sqrt( np.outer(_m_, _m_))

            # Transfer in polarization space
            grad_polbasis = np.conj(self.pols[:,:, iq]).T.dot( grad_nomass.dot(self.pols[:,:,iq]))

            # Fix the modes
            projected_grad = np.zeros(grad_polbasis.shape, dtype = np.complex128)
            #print("MU START:", self.mu_start, "MU END:", self.mu_end)
            projected_grad[self.mu_start:self.mu_end, self.mu_start:self.mu_end] = grad_polbasis[self.mu_start:self.mu_end, self.mu_start:self.mu_end] 

            # Go back in cartesian coordinates
            projected_grad_cart = self.pols[:,:, iq].dot(projected_grad.dot(np.conj(self.pols[:,:,iq]).T))

            # Put the masses again and overwrite the gradient
            new_dyngrad = projected_grad_cart * np.sqrt( np.outer(_m_, _m_))

            # Check if the grad increased 
            if self.testing:
                norm_old = np.sum(np.abs(dyn_grad[iq, :, :])**2)
                norm_new = np.sum(np.abs(new_dyngrad)**2) 

                if norm_new > norm_old:
                    print("Error on q = {}".format(iq))
                    print("Old norm: {} | New norm: {}".format(norm_old, norm_new)) 
                
                assert norm_new < norm_old

            dyn_grad[iq, :, :] = new_dyngrad
                    
            #dyn_grad[iq, :, :] = self.projectorH[iq, :, :].dot(dyn_grad[iq, :, :].dot(self.projector[iq, :, :]))
            # Lets check if the matrix satisfy the sum rule
            #print "DIAG:", np.linalg.eigvalsh(dyn_grad[iq, :, :])
            
    def CFG_ProjectStructure(self, dyn_grad, struct_grad):
        """
        PROJECT ONLY THE STRUCTURE GRADIENT IN THE SELECTED MODES
        =========================================================
        
        This subroutine constraints only the structure gradient, leaving the
        dynamical matrix to minimize on all the possible degrees of freedom.
        """
        
        assert self.initialized, NOT_INITIALIZED_ERROR

        # Project the structure in the polarization vectors
        struct_grad_new = self.proj_vec.dot(struct_grad.ravel())
        struct_grad = struct_grad_new.reshape((self.nat, 3))
            
    def CFG_ProjectDyn(self, dyn_grad, struct_grad):
        """
        PROJECT ONLY THE FC GRADIENT IN THE SELECTED MODES
        ==================================================
        
        This subroutine constrains only the dynamical matrix, leaving the structure
        to minimize on all the possible degrees of freedom.
        """
        assert self.initialized, NOT_INITIALIZED_ERROR
        
        for iq in range(self.nq):
            dyn_grad[iq, :, :] = self.projectorH[iq, :, :].dot(dyn_grad[iq, :, :].dot(self.projector[iq, :, :]))







class IOInfo:
    
    save_weights = False
    weights_file = "weights.dat"
    weights = []
    save_dynmats = False
    ka = 0
    save_dyn_prefix = "minim_dyn"
    
    def __init__(self):
        """
        This class is meant to deal with standard verbose I/O operation,
        like printing the frequencies as a function of the time step of a dynamical matrix (and so on)

        """

        self.total_freqs = []
        self.__save_fname = None
        self.__save_each_step = False

    def Reset(self):
        """
        Reset the data to empty.
        """

        self.__init__()
        
    def SetupWeights(self, fname, save_each_step = True):
        """
        Setup the weights saving
        
        Parameters
        ----------
            fname : string
                Path to the file to which save the frequencies.
        """
        
        self.weights_file = fname
        self.save_weights = True
        self.__save_each_step = save_each_step

    def SetupSaving(self, fname, save_each_step = True):
        """
        Setup the system to save the data each time the function is called.

        Parameters
        ----------
            fname : string 
                path to the file to save the frequencies vs time
            save_each_step : bool
                If true the file is saved (and updated) each time step.

        """

        self.__save_fname = fname
        self.__save_each_step = save_each_step

    def Save(self, fname= None):
        """
        Save the data on a file
        
        Parameters
        ----------
            fname : string, optional
                If given, the file will be saved in the specified location.
                Otherwise the default one is used (must be initialized by SetupSaving)
        """
        if fname is None:
            if self.__save_fname is None:
                raise IOError("Error, a filename must be specified to save the frequencies.")
            np.savetxt(self.__save_fname, self.total_freqs, header = "Time vs Frequencies")
        else:
            np.savetxt(fname, self.total_freqs, header = "Time vs Frequencies")
            
            
        if self.save_weights:
            np.savetxt(self.weights_file, np.transpose(self.weights), header = "Each row is a step containing all the weights of the configurations")
            
        
    def CFP_SaveAll(self, minim):
        """
        This method saves everithing stored in this class.
        
        It can be passed as custom_function_post to the run method of the SchaMinimizer.
        """
        
        # Get the weights if required
        if self.save_weights:
            self.weights.append(minim.ensemble.rho[:])
            
        if self.save_dynmats:
            minim.dyn.save_qe(self.save_dyn_prefix + "_ka%05d_" % self.ka)
            
        # This perform also the saving
        self.CFP_SaveFrequencies(minim)
        
        # Update the step
        self.ka += 1


    def CFP_SaveFrequencies(self, minim):
        """
        This custom method stores the total frequencies updating an exeternal file
        """
        
        # Generate the supercell in real space
        dyn_sc = minim.dyn.GenerateSupercellDyn( minim.ensemble.supercell )

        # Dyagonalize
        w, pols = dyn_sc.DyagDinQ(0)
        self.total_freqs.append(w)
        

        if self.__save_each_step:
            self.Save()



def get_fix_rotations_CFG(dyn):
    """
    FIX ROTATIONS
    =============

    This function returns a pointer to a function that allows to fix the rotational degrees of freedom.

    It can be used when minimizing nanostructures or particles.
    The dynamical matrix must be a gamma point matrix
    """

    assert len(dyn.q_tot) == 1, "Error, only Gamma calculations can fix the rotations."
    assert np.max(np.abs(dyn.q_tot[0])) < 1e-5, "Error, the dynamical matrix must be at Gamma"

    nat = dyn.structure.N_atoms 
    # TODO: To be ended
    
    

def save_binary(object, filename):
    """
    SAVE EVERYTHING
    ===============

    This method saves the whole status of a class (it may be the relax or even the minimizer)
    So that it can be used for analyzing the results. 
    It will recursivly contain all the data stored by the class, included the ensemble and the dynamical matrices.

    The file will be in binary.

    NOTE: There is no warranty that the fill will be readable when loaded with a different python version.

    For this reason, if you want to store the ensemble, the save_bin from ensemble is strongly suggested.

    This methods is just a wrapper for the python pickle utility.

    Parameters
    ----------
        object : anything
            The python object you want to save, it may even be a class. 
        filename : string (path to file)
            The filename on which you want to save the binary data.
    """

    pickle.dump(object, open(filename, "wb"))


def load_binary(filename):
    """
    LOAD EVERYTHING
    ===============

    This method loads the whole status of a class (it may be the relax or even the minimizer)
    So that it can be used for analyzing the results. 
    It will recursivly contain all the data stored by the class, included the ensemble and the dynamical matrices.

    It can read data generated with save_binary

    NOTE: There is no warranty that the fill will be readable when loaded with a different python version.

    This methods is just a wrapper for the python pickle utility.

    Parameters
    ----------
        filename : string (path to file)
            The filename you want to read
    """

    return pickle.load(open(filename, "rb"))