# -*- coding: utf-8 -*-
from __future__ import print_function

"""
This module performs the relax over more
SCHA minimization. It can be both a constant temperature and a 
constant pressure relax. 
"""
import numpy as np
import difflib
import sscha, sscha.Ensemble, sscha.SchaMinimizer
import sscha.Optimizer
import sscha.Calculator
import sscha.Cluster
import sscha.Utilities as Utilities
import cellconstructor as CC
import cellconstructor.symmetries


__EPSILON__ = 1e-5

__RELAX_NAMESPACE__ = "relax"
__RELAX_TYPE__ = "type"
__RELAX_NCONFIGS__ = "n_configs"
__RELAX_MAX_POP__ = "max_pop_id"
__RELAX_START_POP__ = "start_pop"
__RELAX_SAVE_ENSEMBLE__ = "ensemble_datadir"
__RELAX_GENERATE_FIRST_ENSEMBLE__ = "generate_ensemble"
__RELAX_TARGET_PRESSURE__ = "target_pressure"
__RELAX_FIXVOLUME__ = "fix_volume"
__RELAX_BULK_MODULUS__ = "bulk_modulus"

__TYPE_SINGLE__ = "sscha"
__TYPE_RELAX__ = "relax"
__TYPE_VCRELAX__ = "vc-relax"
__ALLOWED_RELAX_TYPES__ = [__TYPE_RELAX__, __TYPE_VCRELAX__]

__REQ_KEYS__ = [__RELAX_TYPE__, __RELAX_NCONFIGS__]
__ALLOWED_KEYS__ = [__RELAX_TYPE__, __RELAX_NCONFIGS__, __RELAX_MAX_POP__,
                    __RELAX_START_POP__, __RELAX_SAVE_ENSEMBLE__,
                    __RELAX_FIXVOLUME__, __RELAX_TARGET_PRESSURE__,
                    __RELAX_BULK_MODULUS__, __RELAX_GENERATE_FIRST_ENSEMBLE__]

class SSCHA(object):
    
    def __init__(self, minimizer = None, ase_calculator=None, N_configs=1, max_pop = 20, 
                 save_ensemble = False, cluster = None):
        """
        This module initialize the relaxer. It may perform
        constant volume or pressure relaxation using fully anharmonic potentials.
        
        Parameters
        ----------
            minimizer : SSCHA_Minimizer
                An initialized SCHA minimizer. Note that its ensemble must be initialized
                with the correct dynamical matrix.
            ase_calculator : ase.calculators...
                An initialized ASE calculator. This will be used to compute energies and forces
                for the relaxation of the SSCHA.
            N_configs : int
                The number of configuration to be used for each population
            max_pop: int, optional
                The maximum number of iteration (The code will stop)
            save_ensemble : bool, optional
                If True (default false) the ensemble is saved after each energy and forces calculation.
            cluster : Cluster.Cluster, optional
                If different from None, the ensemble force and energy calculations
                will be runned in the provided cluster.
        """
        
        if minimizer == None:
            self.minim = sscha.SchaMinimizer.SSCHA_Minimizer()
        else:
            self.minim = minimizer

        self.calc = ase_calculator
        self.N_configs = N_configs
        self.max_pop = max_pop
        self.cluster = cluster
        self.start_pop = 1
        
        # If the ensemble must be saved at each iteration.
        # 
        self.save_ensemble = save_ensemble
        self.data_dir = ""
        
        

        self.__cfpre__ = None
        self.__cfpost__ = None
        self.__cfg__ = None

        # The variable cell attributes
        self.bulk_modulus = 15
        self.target_pressure = 0
        self.fix_volume = False


        # Setup the attribute control
        self.__total_attributes__ = [item for item in self.__dict__.keys()]
        self.fixed_attributes = True # This must be the last attribute to be setted

    def __setattr__(self, name, value):
        """
        This method is used to set an attribute.
        It will raise an exception if the attribute does not exists (with a suggestion of similar entries)
        """

        
        if "fixed_attributes" in self.__dict__:
            if name in self.__total_attributes__:
                super(SSCHA, self).__setattr__(name, value)
            elif self.fixed_attributes:
                similar_objects = str( difflib.get_close_matches(name, self.__total_attributes__))
                ERROR_MSG = """
        Error, the attribute '{}' is not a member of '{}'.
        Suggested similar attributes: {} ?
        """.format(name, type(self).__name__,  similar_objects)

                raise AttributeError(ERROR_MSG)
        else:
            super(SSCHA, self).__setattr__(name, value)
        
        
        
    def setup_from_namelist(self, namelist):
        """
        SETUP THE RELAXER FROM THE NAMELIST
        ===================================
        
        Setup the SSCHA relaxer from the given namelist.
        
        Note the calculation will be also started by this method.
        
        Parameters
        ----------
            namelist : dict
                A dictionary that contains the namespaces
        """
        if not __RELAX_NAMESPACE__ in namelist.keys():
            raise IOError("Error, %s namespace not found" % __RELAX_NAMESPACE__)
        
        c_info = namelist[__RELAX_NAMESPACE__]
        keys = c_info.keys()
        
        # Load the minim and the ensemble if a inputscha namespace is found
        if "inputscha" in namelist.keys():
            self.minim.setup_from_namelist(namelist)
            if self.minim.ensemble is None:
                raise IOError("Error, please specify the input dynamical matrix.")
            
        
        # Load the cluster if any
        if sscha.Cluster.__CLUSTER_NAMELIST__ in namelist.keys():
            self.cluster = sscha.Cluster.Cluster()
            self.cluster.setup_from_namelist(namelist)
        
        # Load the calculator if any
        if sscha.Calculator.__CALCULATOR_HEAD__ in namelist.keys():
            self.calc = sscha.Calculator.prepare_calculator_from_namelist(namelist)
        
        # Get if you must use a already loaded ensemble or not
        restart_from_ens = False
        if __RELAX_GENERATE_FIRST_ENSEMBLE__ in keys:
            restart_from_ens = not bool(c_info[__RELAX_GENERATE_FIRST_ENSEMBLE__])
        
        # Get the number of configurations
        if __RELAX_NCONFIGS__ in keys:
            self.N_configs = int(c_info[__RELAX_NCONFIGS__])
        
        # Get the maximum population
        if __RELAX_MAX_POP__ in keys:
            self.max_pop = int(c_info[__RELAX_MAX_POP__])
            
        if __RELAX_START_POP__ in keys:
            self.start_pop = int(c_info[__RELAX_START_POP__])
        
        if __RELAX_SAVE_ENSEMBLE__ in keys:
            self.save_ensemble = True
            self.data_dir = c_info[__RELAX_SAVE_ENSEMBLE__]
        
        if __RELAX_BULK_MODULUS__ in keys:
            self.bulk_modulus = np.float64(c_info[__RELAX_BULK_MODULUS__])
            
        if __RELAX_TARGET_PRESSURE__ in keys:
            self.target_pressure = np.float64(c_info[__RELAX_TARGET_PRESSURE__])
            
        if __RELAX_FIXVOLUME__ in keys:
            self.fix_volume = bool(c_info[__RELAX_FIXVOLUME__])
        
        
        # Check the allowed keys
        for k in keys: 
            if not k in __ALLOWED_KEYS__:
                print ("Error with the key:", k)
                s = "Did you mean something like:" + str( difflib.get_close_matches(k, __ALLOWED_KEYS__))
                print (s)
                raise IOError("Error in calculator namespace: key '" + k +"' not recognized.\n" + s)
        
        # Check for mandatory keys
        for req_key in __REQ_KEYS__:
            if not req_key in keys:
                raise IOError("Error, the calculator configuration namelist requires the keyword: '" + req_key + "'")
        
        
        
        # Get the calculation type
        if __RELAX_TYPE__ in keys:
            rtype = c_info[__RELAX_TYPE__]
            if not rtype in __ALLOWED_RELAX_TYPES__:
                print ("Unknown relaxation option:", rtype)
                print ("Did you mean:", difflib.get_close_matches(rtype, __ALLOWED_RELAX_TYPES__))
                raise ValueError("Error with key %s" % __RELAX_TYPE__)
        
            # Setup custom functions
            if Utilities.__UTILS_NAMESPACE__ in namelist.keys():
                cfgs = Utilities.get_custom_functions_from_namelist(namelist, self.minim.dyn)
                self.setup_custom_functions(cfgs[0], cfgs[1], cfgs[2])
        
            # Initialize the minimizer
            #self.minim.init()
            self.minim.print_info()
        
            if rtype == __TYPE_RELAX__:
                self.relax(restart_from_ens, self.minim.ensemble.has_stress, self.data_dir,
                           self.start_pop)
            elif rtype == __TYPE_VCRELAX__:
                self.vc_relax(self.target_pressure, self.bulk_modulus, 
                              restart_from_ens, self.data_dir, self.start_pop, fix_volume=self.fix_volume)
        
    def setup_custom_functions(self, custom_function_pre = None,
                               custom_function_gradient = None,
                               custom_function_post = None):
        """
        This subroutine setup which custom functions should be called during the minimization.
        Look for the SCHA_Minimizer.run() method for other details.
        """
        
        self.__cfpre__ = custom_function_pre
        self.__cfpost__ = custom_function_post
        self.__cfg__ = custom_function_gradient
        
        
    def relax(self, restart_from_ens = False, get_stress = False,
              ensemble_loc = None, start_pop = None):
        """
        COSTANT VOLUME RELAX
        ====================
        
        This function performs the costant volume SCHA relaxation, by submitting several populations
        until the minimization converges (or the maximum number of population is reached)
        
        Parameters
        ----------
            restart_from_ens : bool, optional
                If True the ensemble is used to start the first population, without recomputing
                energies and forces. If False (default) the first ensemble is overwritten with
                a new one, and the minimization starts.
            get_stress : bool, optional
                If true the stress tensor is calculated. This may increase the computational
                cost, as it will be computed for each ab-initio configuration (it may be not available
                with some ase calculator)
            ensemble_loc : string
                Where the ensemble of each population is saved on the disk. If none, it will
                use the content of self.data_dir. If also self.data_dir is None, 
                the ensemble will not not be saved (useful to avoid disk I/O for force fields)
            start_pop : int, optional
                The starting index for the population, used only for saving the ensemble and the dynamical 
                matrix. If None, the content of self.start_pop will be used.
            
        Returns
        -------
            status : bool
                True if the minimization converged, False if the maximum number of 
                populations has been reached.
        """

        if ensemble_loc is None:
            ensemble_loc = self.data_dir
        
        if start_pop is None:
            start_pop = self.start_pop

        pop = start_pop
                
        running = True
        while running:
            # Generate the ensemble
            self.minim.ensemble.dyn_0 = self.minim.dyn.Copy()
            
            if pop != start_pop or not restart_from_ens:
                self.minim.ensemble.generate(self.N_configs)
            
                # Compute energies and forces
                self.minim.ensemble.compute_ensemble(self.calc, get_stress, 
                                                 cluster = self.cluster)
                #self.minim.ensemble.get_energy_forces(self.calc, get_stress)
            
                if ensemble_loc is not None and self.save_ensemble:
                    self.minim.ensemble.save_bin(ensemble_loc, pop)
            
            self.minim.population = pop
            self.minim.init(delete_previous_data = False)

            self.minim.run(custom_function_pre = self.__cfpre__,
                           custom_function_post = self.__cfpost__,
                           custom_function_gradient = self.__cfg__)
        
            
            self.minim.finalize()
            
            # Perform the symmetrization
            print ("Checking the symmetries of the dynamical matrix:")
            qe_sym = CC.symmetries.QE_Symmetry(self.minim.dyn.structure)
            qe_sym.SetupQPoint(verbose = True)
            
            print ("Forcing the symmetries in the dynamical matrix.")
            fcq = np.array(self.minim.dyn.dynmats, dtype = np.complex128)
            qe_sym.SymmetrizeFCQ(fcq, self.minim.dyn.q_stars, asr = "custom")
            for iq,q in enumerate(self.minim.dyn.q_tot):
                self.minim.dyn.dynmats[iq] = fcq[iq, :, :]

            # Save the dynamical matrix
            if self.save_ensemble:
                self.minim.dyn.save_qe("dyn_pop%d_" % pop)
        
            # Check if it is converged
            running = not self.minim.is_converged()
            pop += 1
            
            
            if pop > self.max_pop:
                running = False
                
        return self.minim.is_converged()
    
    
    def vc_relax(self, target_press = 0, static_bulk_modulus = 100,
                 restart_from_ens = False,
                 ensemble_loc = ".", start_pop = 1, stress_numerical = False,
                 cell_relax_algorithm = "sd", fix_volume = False):
        """
        VARIABLE CELL RELAX
        ====================
        
        This function performs a variable cell SCHA relaxation at constant pressure,
        It is similar to the relax calculation, but the unit cell is updated according
        to the anharmonic stress tensor at each new population. 
        The cell optimization is performed with the BFGS algorithm. 
        
        NOTE: 
            remember to setup the stress_offset variable of the SCHA_Minimizer,
            because in ab-initio calculation the stress tensor converges porly with the cutoff, 
            but stress tensor differences converges much quicker. Therefore, setup the
            stress tensor difference between a single very high-cutoff calculation and a
            single low-cutoff (the one you use), this difference will be added at the final
            stress tensor to get a better estimation of the true stress.
        
        Parameters
        ----------
            target_press : float, optional
                The target pressure of the minimization (in GPa). The minimization is stopped if the 
                target pressure is the stress tensor is the identity matrix multiplied by the
                target pressure, with a tollerance equal to the stochastic noise. By default 
                it is 0 (ambient pressure)
            static_bulk_modulus : float (default 100), or (9x9) matrix or string, optional
                The static bulk modulus, expressed in GPa. It is used to initialize the
                hessian matrix on the BFGS cell relaxation, to guess the volume deformation caused
                by the anharmonic stress tensor in the first steps. By default is 100 GPa (higher value
                are safer, since they means a lower change in the cell shape).
                It can be also the whole non isotropic matrix. If you specify a string, it 
                can be both:
                    - "recalc" : the static bulk modulus is recomputed with finite differences after
                        each step
                    - "bfgs" : the bfgs algorithm is used to infer the Hessian from previous calculations.
            restart_from_ens : bool, optional
                If True the ensemble is used to start the first population, without recomputing
                energies and forces. If False (default) the first ensemble is overwritten with
                a new one, and the minimization starts.
            ensemble_loc : string
                Where the ensemble of each population is saved on the disk. You can specify None
                if you do not want to save the ensemble (useful to avoid disk I/O for force fields)
            start_pop : int, optional
                The starting index for the population, used only for saving the ensemble and the dynamical 
                matrix.
            stress_numerical : bool
                If True the stress is computed by finite difference (usefull for calculators that 
                does not support it by default)
            cell_relax_algorithm : string
                This identifies the stress algorithm. It can be both sd (steepest-descent),
                cg (conjugate-gradient) or bfgs (Quasi-newton)
            fix_volume : bool, optional
                If true (default False) the volume is fixed, therefore only the cell shape is relaxed.
                
            
        Returns
        -------
            status : bool
                True if the minimization converged, False if the maximum number of 
                populations has been reached.
        """
        # Rescale the target pressure in eV / A^3
        target_press_evA3 = target_press / sscha.SchaMinimizer.__evA3_to_GPa__
        I = np.eye(3, dtype = np.float64)
        
        SUPPORTED_ALGORITHMS = ["sd", "cg", "bfgs"]
        if not cell_relax_algorithm in SUPPORTED_ALGORITHMS:
            raise ValueError("Error, cell_relax_algorithm %s not supported." %  cell_relax_algorithm)

        # Read the bulk modulus
        kind_minimizer = "SD"
        if type(static_bulk_modulus) == type(""):
            if static_bulk_modulus == "recalc":
                kind_minimizer = "RPSD"
            elif static_bulk_modulus == "none":
                kind_minimizer = "SD"
                static_bulk_modulus = 100
            elif static_bulk_modulus == "bfgs":
                static_bulk_modulus = 100
                kind_minimizer = "BFGS"
            else:
                raise ValueError("Error, value '%s' not supported for bulk modulus." % static_bulk_modulus)
        elif len(np.shape(static_bulk_modulus)) == 0:
            kind_minimizer = cell_relax_algorithm.upper()
        elif len(np.shape(static_bulk_modulus)) == 2:
            kind_minimizer = "PSD"
        else:
            raise ValueError("Error, the given value not supported as a bulk modulus.")
            
        
        

        if static_bulk_modulus == "recalc":
            # Rescale the static bulk modulus in eV / A^3
            static_bulk_modulus /= sscha.SchaMinimizer.__evA3_to_GPa__ 

        # initilaize the cell minimizer
        #BFGS = sscha.Optimizer.BFGS_UC(self.minim.dyn.structure.unit_cell, static_bulk_modulus)
        if kind_minimizer == "SD":
            BFGS = sscha.Optimizer.UC_OPTIMIZER(self.minim.dyn.structure.unit_cell)
            BFGS.alpha = 1 / (3 * static_bulk_modulus * np.linalg.det(self.minim.dyn.structure.unit_cell))
        if kind_minimizer == "CG":
            BFGS = sscha.Optimizer.CG_UC(self.minim.dyn.structure.unit_cell)
            BFGS.alpha = 1 / (3 * static_bulk_modulus * np.linalg.det(self.minim.dyn.structure.unit_cell))
        elif kind_minimizer == "PSD":
            BFGS = sscha.Optimizer.SD_PREC_UC(self.minim.dyn.structure.unit_cell, static_bulk_modulus)
        elif kind_minimizer == "BFGS":
            BFGS = sscha.Optimizer.BFGS_UC(self.minim.dyn.structure.unit_cell, static_bulk_modulus)

        # Initialize the bulk modulus
        # The gradient (stress) is in eV/A^3, we have the cell in Angstrom so the Hessian must be
        # in eV / A^6

        pop = start_pop
                
        running = True
        while running:
            # Compute the static bulk modulus if required
            if kind_minimizer == "RPSD":
                # Compute the static bulk modulus
                sbm = GetStaticBulkModulus(self.minim.dyn.structure, self.calc)
                print ("BM:")
                print (sbm)
                BFGS = sscha.Optimizer.SD_PREC_UC(self.minim.dyn.structure.unit_cell, sbm)

            # Generate the ensemble
            self.minim.ensemble.dyn_0 = self.minim.dyn.Copy()
            if pop != start_pop or not restart_from_ens:
                self.minim.ensemble.generate(self.N_configs)
            
                # Compute energies and forces
                self.minim.ensemble.compute_ensemble(self.calc, True, stress_numerical,
                                                 cluster = self.cluster)
                #self.minim.ensemble.get_energy_forces(self.calc, True, stress_numerical = stress_numerical)
            
                if ensemble_loc is not None and self.save_ensemble:
                    self.minim.ensemble.save_bin(ensemble_loc, pop)
            
            self.minim.population = pop
            self.minim.init(delete_previous_data = False)
        
            self.minim.run(custom_function_pre = self.__cfpre__,
                           custom_function_post = self.__cfpost__,
                           custom_function_gradient = self.__cfg__)
        
            
            self.minim.finalize()
            
            # Get the stress tensor [ev/A^3]
            stress_tensor, stress_err = self.minim.get_stress_tensor() 
            stress_tensor *= sscha.SchaMinimizer.__RyBohr3_to_evA3__
            stress_err *=  sscha.SchaMinimizer.__RyBohr3_to_evA3__

            # Get the pressure
            Press = np.trace(stress_tensor) / 3
            
            # Get the volume
            Vol = np.linalg.det(self.minim.dyn.structure.unit_cell)
            
            # Get the Helmoltz-Gibbs free energy
            helmoltz = self.minim.get_free_energy() * sscha.SchaMinimizer.__RyToev__
            gibbs = helmoltz + target_press_evA3 * Vol - self.minim.eq_energy
            
            # Prepare a mark to underline which quantity is actually minimized by the
            # Variable relaxation algorithm if the helmoltz free energy (in case of fixed volume)
            # Or the Gibbs free energy (in case of fixed pressure)
            mark_helmoltz = ""
            mark_gibbs = ""
            if fix_volume:
                mark_helmoltz = "<--"
            else:
                mark_helmoltz = "<--"

            # Print the enthalpic contribution
            message = """
 ====================== 
 ENTHALPIC CONTRIBUTION 
 ====================== 

 P = {:.4f} GPa   V = {:.4f} A^3

 P V = {:.8e} eV

 Helmoltz Free energy = {:.8e} eV {}
 Gibbs Free energy = {:.8e} eV {}
 Zero energy = {:.8e} eV

 """.format(target_press , Vol,target_press_evA3 * Vol, helmoltz, mark_helmoltz, gibbs, mark_gibbs, self.minim.eq_energy)
            print(message)
            # print " ====================== "
            # print " ENTHALPIC CONTRIBUTION "
            # print " ====================== "
            # print ""
            # print "  P = %.4f GPa    V = %.4f A^3" % (target_press , Vol)
            # print ""
            # print "  P V = %.8e eV " % (target_press_evA3 * Vol)
            # print ""
            # print " Helmoltz Free energy = %.8e eV " % helmoltz,
            # if fix_volume:
            #     print "  <-- "
            # else:
            #     print ""
            # print " Gibbs Free energy = %.8e eV " % gibbs,
            # if fix_volume:
            #     print ""
            # else:
            #     print "  <-- "
            # print " (Zero energy = %.8e eV) " % self.minim.eq_energy
            # print ""
        
            # Perform the cell step
            cell_gradient = (stress_tensor - I *target_press_evA3)
            
            new_uc = self.minim.dyn.structure.unit_cell.copy()
            BFGS.UpdateCell(new_uc,  cell_gradient, fix_volume)
            
            # Strain the structure and the q points preserving the symmetries
            self.minim.dyn.AdjustToNewCell(new_uc)
            #self.minim.dyn.structure.change_unit_cell(new_uc)
            

            print (" New unit cell:")
            print (" v1 [A] = (%16.8f %16.8f %16.8f)" % (new_uc[0,0], new_uc[0,1], new_uc[0,2]))
            print (" v2 [A] = (%16.8f %16.8f %16.8f)" % (new_uc[1,0], new_uc[1,1], new_uc[1,2]))
            print (" v3 [A] = (%16.8f %16.8f %16.8f)" % (new_uc[2,0], new_uc[2,1], new_uc[2,2]))
            
            print ()
            print ("Check the symmetries in the new cell:")
            qe_sym = CC.symmetries.QE_Symmetry(self.minim.dyn.structure)
            qe_sym.SetupQPoint(verbose = True)
            
            print ("Forcing the symmetries in the dynamical matrix.")
            fcq = np.array(self.minim.dyn.dynmats, dtype = np.complex128)
            qe_sym.SymmetrizeFCQ(fcq, self.minim.dyn.q_stars, asr = "custom")
            for iq,q in enumerate(self.minim.dyn.q_tot):
                self.minim.dyn.dynmats[iq] = fcq[iq, :, :]

            # Save the dynamical matrix
            self.minim.dyn.save_qe("dyn_pop%d_" % pop)

            # Check if the constant volume calculation is converged
            running1 = not self.minim.is_converged()

            # Check if the cell variation is converged
            running2 = True
            not_zero_mask = stress_err != 0
            if not fix_volume:
                if np.max(np.abs(cell_gradient[not_zero_mask]) / stress_err[not_zero_mask]) <= 1:
                    running2 = False
            else:
                if np.max(np.abs((stress_tensor - I * Press)[not_zero_mask] / 
                                 stress_err[not_zero_mask])) <= 1:
                    running2 = False
                

            running = running1 or running2

            pop += 1
            
            if pop > self.max_pop:
                running = False
                
        return (not running1) and (not running2)


def GetStaticBulkModulus(structure, ase_calculator, eps = 1e-3):
    """
    GET STATIC BULK MODULUS
    =======================

    This method uses finite differences on the cell to compute
    the static bulk modulus. The cell is strained into several volumes,
    and the stress tensor is computed in orther to obtain the bulk modulus.
    Only the symmmetry relevant terms are computed.

    Parameters
    ----------
        structure : CC.Structure.Structure()
            The structure on which you want to compute the static bulk modulus
        ase_calculator : ase.calculators.calculator.Calculator()
            One of the ase calculators to get the stress tensor in several strained
            cells.
        eps : float
            The strain module

    Results
    -------
        bk_mod : ndarray (9x9)
            The bulk modulus as a 9x9 matrix, expressed in eV / A^3
    """

    # Initialize the symmetries
    qe_sym = CC.symmetries.QE_Symmetry(structure)

    # Perform the firts calculation
    atm_center = structure.get_ase_atoms()
    atm_center.set_calculator(ase_calculator)
    V_0 = np.linalg.det(structure.unit_cell)
    stress_0 = atm_center.get_stress(False)
    I = np.eye(3, dtype = np.float64)
    

    qe_sym.ApplySymmetryToMatrix(stress_0)

    bk_mod = np.zeros( (9,9), dtype = np.float64)
    # Get the non zero elements
    for i in range(3):
        for j in range(i, 3):
            trials = np.zeros((3,3), dtype = np.float64)
            trials[i,j] = 1
            qe_sym.ApplySymmetryToMatrix(trials)


            if trials[i,j] == 0:
                continue


            mask = (np.abs(trials) < __EPSILON__).astype(np.float64)

            strain = eps * mask
            
            dstrain = eps * np.sqrt(np.sum(mask))

            # Strain the cell and perform the calculation
            new_structure = structure.copy()
            new_structure.change_unit_cell(structure.unit_cell.dot(I + strain))

            atm = new_structure.get_ase_atoms()
            atm.set_calculator(ase_calculator)
            stress_1 = atm.get_stress(False)
            V_1 = np.linalg.det(new_structure.unit_cell)


            qe_sym.ApplySymmetryToMatrix(stress_1)

            bk_mod[3*i + j, :] = -(stress_1 - stress_0).ravel() / dstrain
            if j!= i:
                bk_mod[3*j + i, :] = -(stress_1 - stress_0).ravel() / dstrain
            
            # Apply hermitianity
            bk_mod = 0.5 * (bk_mod.transpose() + bk_mod)

    return bk_mod





