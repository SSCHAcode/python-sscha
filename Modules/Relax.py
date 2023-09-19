# -*- coding: utf-8 -*-
from __future__ import print_function, annotations

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
from sscha.aiida_ensemble import AiiDAEnsemble

import sys, os

from sscha.Parallel import pprint as print

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
__RELAX_SOBOL__ = "sobol_sampling"
__RELAX_SOBOL_SCATTER__ = "sobol_scatter"

__TYPE_SINGLE__ = "sscha"
__TYPE_RELAX__ = "relax"
__TYPE_VCRELAX__ = "vc-relax"
__ALLOWED_RELAX_TYPES__ = [__TYPE_RELAX__, __TYPE_VCRELAX__]

__REQ_KEYS__ = [__RELAX_TYPE__, __RELAX_NCONFIGS__]
__ALLOWED_KEYS__ = [__RELAX_TYPE__, __RELAX_NCONFIGS__, __RELAX_MAX_POP__,
                    __RELAX_START_POP__, __RELAX_SAVE_ENSEMBLE__,
                    __RELAX_FIXVOLUME__, __RELAX_TARGET_PRESSURE__,
                    __RELAX_BULK_MODULUS__, __RELAX_GENERATE_FIRST_ENSEMBLE__,
                    __RELAX_SOBOL__, __RELAX_SOBOL_SCATTER__]

class SSCHA(object):

    def __init__(
        self,
        minimizer = None,
        ase_calculator=None, 
        aiida_inputs: dict | None = None,
        N_configs=1, 
        max_pop = 20,
        save_ensemble = False, 
        cluster = None,
        **kwargs
    ):
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
            aiida_input: dict
                Dictionary containing the input data for the 
                :class:`~sscha.aiida_ensemble.AiiDAEnsemble.compute_ensmble`
                method; namely:
                    * pw_code: str,
                    * protocol: str = 'moderate',
                    * options: dict = None,
                    * overrides: dict = None,
                    * group_label: str = None,
                    * kwargs
            N_configs : int
                The number of configuration to be used for each population
            max_pop: int, optional
                The maximum number of iteration (The code will stop)
            save_ensemble : bool, optional
                If True (default false) the ensemble is saved after each energy and forces calculation.
            cluster : Cluster.Cluster, optional
                If different from None, the ensemble force and energy calculations
                will be runned in the provided cluster.
            **kwargs : any other keyword that matches an object of this structure
        """

        if minimizer == None:
            self.minim = sscha.SchaMinimizer.SSCHA_Minimizer()
        else:
            self.minim = minimizer

        self.calc = ase_calculator
        self.N_configs = N_configs
        self.max_pop = max_pop
        self.cluster = cluster
        self.aiida_inputs = aiida_inputs
        self.start_pop = 1

        # If the ensemble must be saved at each iteration.
        #
        self.save_ensemble = save_ensemble
        self.data_dir = "data"



        self.__cfpre__ = None
        self.__cfpost__ = None
        self.__cfg__ = None

        # The variable cell attributes
        self.bulk_modulus = 15
        self.target_pressure = 0
        self.fix_volume = False

        # Options for the cell relaxation
        # If true the cell shape is fixed in a variable-cell relaxation
        # even if the symmetries allows for more degrees of freedom in the cell shape.
        # Usefull (for example) if you want to enfoce a cubic cell even if the structure brakes the symmetries
        self.fix_cell_shape = False


        # Set the Sobol Parameters by default (aka no Sobol)
        self.sobol_sampling = False
        self.sobol_scatter = 0.0


        # Setup the attribute control
        self.__total_attributes__ = [item for item in self.__dict__.keys()]
        self.fixed_attributes = True # This must be the last attribute to be setted


        # Setup any other keyword given in input (raising the error if not already defined)
        for key in kwargs:
            self.__setattr__(key, kwargs[key])


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

# ****Diegom_test****
        if __RELAX_SOBOL__ in keys:
            self.sobol_sampling = bool(c_info[__RELAX_SOBOL__])

        if __RELAX_SOBOL_SCATTER__ in keys:
            self.sobol_scatter = np.float64(c_info[__RELAX_SOBOL_SCATTER__])

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
                           self.start_pop, sobol = self.sobol_sampling, sobol_scatter = self.sobol_scatter)
            elif rtype == __TYPE_VCRELAX__:
                self.vc_relax(self.target_pressure, self.bulk_modulus,
                              restart_from_ens, self.data_dir, self.start_pop, fix_volume=self.fix_volume,
                               sobol = self.sobol_sampling, sobol_scatter = self.sobol_scatter)

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
              ensemble_loc = None, start_pop = None, sobol = False,
               sobol_scramble = False, sobol_scatter = 0.0):
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
                use the content of self.data_dir. It is just a way to override the variable self.data_dir
            start_pop : int, optional
                The starting index for the population, used only for saving the ensemble and the dynamical
                matrix. If None, the content of self.start_pop will be used.
            sobol : bool, optional (Default = False)
                 Defines if the calculation uses random Gaussian generator or Sobol Gaussian generator.
            sobol_scramble : bool, optional (Default = False)
                Set the optional scrambling of the generated numbers taken from the Sobol sequence.
            sobol_scatter : real (0.0 to 1) (Deafault = 0.0)
                Set the scatter parameter to displace the Sobol positions randommly.

        Returns
        -------
            status : bool
                True if the minimization converged, False if the maximum number of
                populations has been reached.
        """

        if ensemble_loc is None:
            ensemble_loc = self.data_dir

        if (not ensemble_loc) and self.save_ensemble:
            ERR_MSG = """
Error, you must specify where to save the ensembles.
       this can be done either passing ensemble_loc = "path/to/dir"
       for the ensemble, or by setting the data_dir attribute of this object.
"""
            raise IOError(ERR_MSG)

        if self.save_ensemble:
            if not os.path.exists(ensemble_loc):
                os.makedirs(ensemble_loc)
            else:
                if not os.path.isdir(ensemble_loc):
                    ERR_MSG = """
Error, the specified location to save the ensemble:
       '{}'
       already exists and it is not a directory.
""".format(ensemble_loc)
                    raise IOError(ERR_MSG)


        if start_pop is None:
            start_pop = self.start_pop

        pop = start_pop

        running = True
        while running:
            # Generate the ensemble
            self.minim.ensemble.dyn_0 = self.minim.dyn.Copy()

            if pop != start_pop or not restart_from_ens:
                self.minim.ensemble.generate(self.N_configs, sobol = sobol, sobol_scramble = sobol_scramble, sobol_scatter = sobol_scatter)

                # Compute energies and forces
                if isinstance(self.minim.ensemble, AiiDAEnsemble):
                    self.minim.ensemble.compute_ensemble(**self.aiida_inputs)
                else:
                    self.minim.ensemble.compute_ensemble(
                        self.calc, get_stress, cluster = self.cluster)
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


        self.start_pop = pop
        print('Population = ',pop) #**** Diegom_test ****
        return self.minim.is_converged()


    def vc_relax(self, target_press = 0, static_bulk_modulus = 100,
                 restart_from_ens = False,
                 ensemble_loc = None, start_pop = None, stress_numerical = False,
                 cell_relax_algorithm = "sd", fix_volume = False, sobol = False,
                  sobol_scramble = False, sobol_scatter = 0.0):
        """
        VARIABLE CELL RELAX
        ====================

        This function performs a variable cell SCHA relaxation at constant pressure,
        It is similar to the relax calculation, but the unit cell is updated according
        to the anharmonic stress tensor at each new population.

        By default, all the degrees of freedom compatible with the symmetry group are relaxed in the cell.
        You can constrain the cell to keep the same shape by setting fix_cell_shape = True.


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
                cg (conjugate-gradient) or bfgs (Quasi-newton).
                The most robust one is SD. Do not change if you are not sure what you are doing.
            fix_volume : bool, optional
                If true (default False) the volume is fixed, therefore only the cell shape is relaxed.
            sobol : bool, optional (Default = False)
                 Defines if the calculation uses random Gaussian generator or Sobol Gaussian generator.
            sobol_scramble : bool, optional (Default = False)
                Set the optional scrambling of the generated numbers taken from the Sobol sequence.
            sobol_scatter : real (0.0 to 1) (Deafault = 0.0)
                Set the scatter parameter to displace the Sobol positions randommly.

        Returns
        -------
            status : bool
                True if the minimization converged, False if the maximum number of
                populations has been reached.
        """

        # Prepare the saving directory
        if ensemble_loc is None:
            ensemble_loc = self.data_dir

        if (not ensemble_loc) and self.save_ensemble:
            ERR_MSG = """
Error, you must specify where to save the ensembles.
       this can be done either passing ensemble_loc = "path/to/dir"
       for the ensemble, or by setting the data_dir attribute of this object.
"""
            raise IOError(ERR_MSG)

        if self.save_ensemble:
            if not os.path.exists(ensemble_loc):
                os.makedirs(ensemble_loc)
            else:
                if not os.path.isdir(ensemble_loc):
                    ERR_MSG = """
Error, the specified location to save the ensemble:
       '{}'
       already exists and it is not a directory.
""".format(ensemble_loc)
                    raise IOError(ERR_MSG)



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




        if static_bulk_modulus != "recalc":
            # Rescale the static bulk modulus in eV / A^3
            static_bulk_modulus /= sscha.SchaMinimizer.__evA3_to_GPa__

        # initilaize the cell minimizer
        #BFGS = sscha.Optimizer.BFGS_UC(self.minim.dyn.structure.unit_cell, static_bulk_modulus)
        if kind_minimizer in ["SD", "CG"] :
            BFGS = sscha.Optimizer.UC_OPTIMIZER(self.minim.dyn.structure.unit_cell)
            BFGS.alpha = 1 / (3 * static_bulk_modulus * self.minim.dyn.structure.get_volume())
            BFGS.algorithm = kind_minimizer.lower()
        elif kind_minimizer == "PSD":
            BFGS = sscha.Optimizer.SD_PREC_UC(self.minim.dyn.structure.unit_cell, static_bulk_modulus)
        elif kind_minimizer == "BFGS":
            BFGS = sscha.Optimizer.BFGS_UC(self.minim.dyn.structure.unit_cell, static_bulk_modulus)

        # Initialize the bulk modulus
        # The gradient (stress) is in eV/A^3, we have the cell in Angstrom so the Hessian must be
        # in eV / A^6
        if start_pop is not None:
            pop = start_pop
        else:
            pop = self.start_pop

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
                self.minim.ensemble.generate(self.N_configs, sobol=sobol, sobol_scramble = sobol_scramble, sobol_scatter = sobol_scatter)

                # Save also the generation
                #if ensemble_loc is not None and self.save_ensemble:
                #    self.minim.ensemble.save_bin(ensemble_loc, pop)

                # Compute energies and forces
                if isinstance(self.minim.ensemble, AiiDAEnsemble):
                    self.minim.ensemble.compute_ensemble(**self.aiida_inputs)
                else:
                    self.minim.ensemble.compute_ensemble(
                        self.calc, True, stress_numerical, cluster = self.cluster)
                #self.minim.ensemble.get_energy_forces(self.calc, True, stress_numerical = stress_numerical)

                print("RELAX force length:", len(self.minim.ensemble.force_computed))
                
                if ensemble_loc is not None and self.save_ensemble:
                    self.minim.ensemble.save_bin(ensemble_loc, pop)
                print("RELAX force length:", len(self.minim.ensemble.force_computed))

            self.minim.population = pop
            self.minim.init(delete_previous_data = False)

            print("RELAX force length:", len(self.minim.ensemble.force_computed))
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
            Vol = self.minim.dyn.structure.get_volume()

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
                mark_gibbs = "<--"

            # Extract the bulk modulus from the cell minimization
            new_bulkmodulus = 1 / (3 * BFGS.alpha * self.minim.dyn.structure.get_volume())
            new_bulkmodulus *= sscha.SchaMinimizer.__evA3_to_GPa__

            # Print the enthalpic contribution
            message = """
 ======================
 ENTHALPIC CONTRIBUTION
 ======================

 P = {:.4f} GPa   V = {:.4f} A^3

 P V = {:.8e} eV

 Helmoltz Free energy = {:.10e} eV {}
 Gibbs Free energy = {:.10e} eV {}
 Zero energy = {:.10e} eV

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
            if self.fix_cell_shape:
                # Use a isotropic stress tensor to keep the same cell shape
                cell_gradient = I * (Press - target_press_evA3)
            else:
                cell_gradient = (stress_tensor - I *target_press_evA3)

            new_uc = self.minim.dyn.structure.unit_cell.copy()
            BFGS.UpdateCell(new_uc,  cell_gradient, fix_volume)

            # Strain the structure and the q points preserving the symmetries
            self.minim.dyn.AdjustToNewCell(new_uc)
            #self.minim.dyn.structure.change_unit_cell(new_uc)

            message = """
 Currently estimated bulk modulus = {:8.3f} GPa
 (Note: this is just indicative, do not use it for computing bulk modulus)

 """.format(new_bulkmodulus)
            print(message)


            print (" New unit cell:")
            print (" v1 [A] = (%16.8f %16.8f %16.8f)" % (new_uc[0,0], new_uc[0,1], new_uc[0,2]))
            print (" v2 [A] = (%16.8f %16.8f %16.8f)" % (new_uc[1,0], new_uc[1,1], new_uc[1,2]))
            print (" v3 [A] = (%16.8f %16.8f %16.8f)" % (new_uc[2,0], new_uc[2,1], new_uc[2,2]))

            print ()
            print ("Check the symmetries in the new cell:")
            sys.stdout.flush()
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

        self.start_pop = pop
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





































import ase
from ase.calculators.lammpsrun import LAMMPS
from ase.units import Bohr, Ry
from datetime import datetime
from cellconstructor.calculators import Espresso


class SSCHA_MTP(SSCHA):
    def __init__(self, minimizer = None, ase_calculator=None, N_configs=1, max_pop = 20,
                 save_ensemble = False, cluster = None, 
                 specorder = None,
                 pot_name = 'fitted.mtp', 
                 mlip_run_command = 'mpirun -np 1', 
                 path_to_mlip = 'mlp',
                 ab_initio_calculator = 'QE',
                 ab_initio_parameters = None,
                 ab_initio_modules_load = None,
                 ab_initio_kresol = 0.25,
                 ab_initio_pseudos = None,
                 ab_initio_cluster = None,
                 iteration_limit = 500,
                 energy_weight = 1.0,
                 force_weight = 0.01,
                 stress_weight = 0.001,
                 include_stress = False,
                 train_on_every_ensemble = False,
                 train_local_mtps = False,
                 retrain = False,
                 **kwargs):

        # super().__init__(minimizer=minimizer, ase_calculator=ase_calculator, N_configs=N_configs, max_pop=max_pop,
        #          save_ensemble = save_ensemble, cluster = cluster)

        # self.minimizer = minimizer
        # self.ase_calculator = ase_calculator
        # self.N_configs = N_configs
        # self.max_pop = max_pop
        # self.save_ensemble = save_ensemble
        # self.cluster = cluster

        self.specorder = specorder
        self.mlip_run_command = mlip_run_command
        self.path_to_mlip = path_to_mlip
        self.pot_name = pot_name 
        self.ab_initio_calculator = ab_initio_calculator
        self.ab_initio_parameters = ab_initio_parameters 
        self.ab_initio_modules_load = ab_initio_modules_load
        self.ab_initio_kresol = ab_initio_kresol 
        self.ab_initio_pseudos = ab_initio_pseudos 
        self.ab_initio_cluster = ab_initio_cluster
        self.iteration_limit = iteration_limit 
        self.energy_weight = energy_weight 
        self.force_weight = force_weight
        self.stress_weight = stress_weight
        self.include_stress = include_stress
        self.train_on_every_ensemble = train_on_every_ensemble # if True the MTP is trained every time the new ensemble is generated  
        self.train_local_mtps = train_local_mtps # if True the new MTP is trained from scratch every time the training new set is generated (e.g. from ensemble)  
        self.retrain = retrain # IMPORTANT TAG!!! if we wish to retrain the MTP on structures produced with extrapolation control

        # super().__init__(minimizer=self.minimizer, ase_calculator=self.ase_calculator, N_configs=self.N_configs, max_pop=self.max_pop,
        #          save_ensemble = self.save_ensemble, cluster = self.cluster)

        super().__init__(minimizer, ase_calculator, N_configs, max_pop,
                 save_ensemble, cluster)


    def relax(self, restart_from_ens = False, get_stress = False,
              ensemble_loc = None, start_pop = None, sobol = False,
              sobol_scramble = False, sobol_scatter = 0.0):#,
            #   train_on_every_ensemble = False,
            #   train_local_mtps = False,
            #   retrain = False):
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
                use the content of self.data_dir. It is just a way to override the variable self.data_dir
            start_pop : int, optional
                The starting index for the population, used only for saving the ensemble and the dynamical
                matrix. If None, the content of self.start_pop will be used.
            sobol : bool, optional (Default = False)
                 Defines if the calculation uses random Gaussian generator or Sobol Gaussian generator.
            sobol_scramble : bool, optional (Default = False)
                Set the optional scrambling of the generated numbers taken from the Sobol sequence.
            sobol_scatter : real (0.0 to 1) (Deafault = 0.0)
                Set the scatter parameter to displace the Sobol positions randommly.

        Returns
        -------
            status : bool
                True if the minimization converged, False if the maximum number of
                populations has been reached.
        """

        if ensemble_loc is None:
            ensemble_loc = self.data_dir

        if (not ensemble_loc) and self.save_ensemble:
            ERR_MSG = """
Error, you must specify where to save the ensembles.
       this can be done either passing ensemble_loc = "path/to/dir"
       for the ensemble, or by setting the data_dir attribute of this object.
"""
            raise IOError(ERR_MSG)

        if self.save_ensemble:
            if not os.path.exists(ensemble_loc):
                os.makedirs(ensemble_loc)
            else:
                if not os.path.isdir(ensemble_loc):
                    ERR_MSG = """
Error, the specified location to save the ensemble:
       '{}'
       already exists and it is not a directory.
""".format(ensemble_loc)
                    raise IOError(ERR_MSG)


        if start_pop is None:
            start_pop = self.start_pop

        pop = start_pop

        running = True
        while running:
            # Generate the ensemble
            self.minim.ensemble.dyn_0 = self.minim.dyn.Copy()

            if pop != start_pop or not restart_from_ens:
                self.minim.ensemble.generate(self.N_configs, sobol = sobol, sobol_scramble = sobol_scramble, sobol_scatter = sobol_scatter)
                
                # Train MTP on generated ensemble
                if self.train_on_every_ensemble:
                    ensemble_structures = self.minim.ensemble.structures
                    ase_structures_list_to_cfg(ensemble_structures,'preselected.cfg')
                    os.system('touch set.cfg')
                    train_mtp_on_cfg(self.specorder,self.mlip_run_command, self.path_to_mlip, self.pot_name, 
                        self.ab_initio_calculator, self.ab_initio_parameters, self.ab_initio_modules_load, self.ab_initio_kresol, self.ab_initio_pseudos, self.ab_initio_cluster, 
                        self.iteration_limit, self.energy_weight, self.force_weight, self.stress_weight, self.include_stress, self.train_local_mtps,pop)
                elif not self.train_on_every_ensemble:
                    if pop == start_pop:
                        ensemble_structures = self.minim.ensemble.structures
                        ase_structures_list_to_cfg(ensemble_structures,'preselected.cfg')
                        os.system('touch set.cfg')
                        train_mtp_on_cfg(self.specorder,self.mlip_run_command, self.path_to_mlip, self.pot_name, 
                            self.ab_initio_calculator, self.ab_initio_parameters, self.ab_initio_modules_load, self.ab_initio_kresol, self.ab_initio_pseudos, self.ab_initio_cluster, 
                            self.iteration_limit, self.energy_weight, self.force_weight, self.stress_weight, self.include_stress, self.train_local_mtps,pop)
                    else:
                        pass

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


        self.start_pop = pop
        print('Population = ',pop) #**** Diegom_test ****
        return self.minim.is_converged()


    def vc_relax(self, target_press = 0, static_bulk_modulus = 100,
                 restart_from_ens = False,
                 ensemble_loc = None, start_pop = None, stress_numerical = False,
                 cell_relax_algorithm = "sd", fix_volume = False, sobol = False,
                 sobol_scramble = False, sobol_scatter = 0.0):#, 
                #  train_on_every_ensemble = False,
                #  train_local_mtps = False,
                #  retrain = False):
        """
        VARIABLE CELL RELAX
        ====================

        This function performs a variable cell SCHA relaxation at constant pressure,
        It is similar to the relax calculation, but the unit cell is updated according
        to the anharmonic stress tensor at each new population.

        By default, all the degrees of freedom compatible with the symmetry group are relaxed in the cell.
        You can constrain the cell to keep the same shape by setting fix_cell_shape = True.


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
                cg (conjugate-gradient) or bfgs (Quasi-newton).
                The most robust one is SD. Do not change if you are not sure what you are doing.
            fix_volume : bool, optional
                If true (default False) the volume is fixed, therefore only the cell shape is relaxed.
            sobol : bool, optional (Default = False)
                 Defines if the calculation uses random Gaussian generator or Sobol Gaussian generator.
            sobol_scramble : bool, optional (Default = False)
                Set the optional scrambling of the generated numbers taken from the Sobol sequence.
            sobol_scatter : real (0.0 to 1) (Deafault = 0.0)
                Set the scatter parameter to displace the Sobol positions randommly.

        Returns
        -------
            status : bool
                True if the minimization converged, False if the maximum number of
                populations has been reached.
        """

        # Prepare the saving directory
        if ensemble_loc is None:
            ensemble_loc = self.data_dir

        if (not ensemble_loc) and self.save_ensemble:
            ERR_MSG = """
Error, you must specify where to save the ensembles.
       this can be done either passing ensemble_loc = "path/to/dir"
       for the ensemble, or by setting the data_dir attribute of this object.
"""
            raise IOError(ERR_MSG)

        if self.save_ensemble:
            if not os.path.exists(ensemble_loc):
                os.makedirs(ensemble_loc)
            else:
                if not os.path.isdir(ensemble_loc):
                    ERR_MSG = """
Error, the specified location to save the ensemble:
       '{}'
       already exists and it is not a directory.
""".format(ensemble_loc)
                    raise IOError(ERR_MSG)



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




        if static_bulk_modulus != "recalc":
            # Rescale the static bulk modulus in eV / A^3
            static_bulk_modulus /= sscha.SchaMinimizer.__evA3_to_GPa__

        # initilaize the cell minimizer
        #BFGS = sscha.Optimizer.BFGS_UC(self.minim.dyn.structure.unit_cell, static_bulk_modulus)
        if kind_minimizer in ["SD", "CG"] :
            BFGS = sscha.Optimizer.UC_OPTIMIZER(self.minim.dyn.structure.unit_cell)
            BFGS.alpha = 1 / (3 * static_bulk_modulus * self.minim.dyn.structure.get_volume())
            BFGS.algorithm = kind_minimizer.lower()
        elif kind_minimizer == "PSD":
            BFGS = sscha.Optimizer.SD_PREC_UC(self.minim.dyn.structure.unit_cell, static_bulk_modulus)
        elif kind_minimizer == "BFGS":
            BFGS = sscha.Optimizer.BFGS_UC(self.minim.dyn.structure.unit_cell, static_bulk_modulus)

        # Initialize the bulk modulus
        # The gradient (stress) is in eV/A^3, we have the cell in Angstrom so the Hessian must be
        # in eV / A^6
        if start_pop is not None:
            pop = start_pop
        else:
            pop = self.start_pop

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
                self.minim.ensemble.generate(self.N_configs, sobol=sobol, sobol_scramble = sobol_scramble, sobol_scatter = sobol_scatter)

                # Save also the generation
                #if ensemble_loc is not None and self.save_ensemble:
                #    self.minim.ensemble.save_bin(ensemble_loc, pop)
                # Train MTP on generated ensemble
                if self.train_on_every_ensemble:
                    ensemble_structures = self.minim.ensemble.structures
                    ase_structures_list_to_cfg(ensemble_structures,'preselected.cfg')
                    os.system('touch set.cfg')
                    train_mtp_on_cfg(self.specorder,self.mlip_run_command, self.path_to_mlip, self.pot_name, 
                        self.ab_initio_calculator, self.ab_initio_parameters, self.ab_initio_modules_load, self.ab_initio_kresol, self.ab_initio_pseudos, self.ab_initio_cluster, 
                        self.iteration_limit, self.energy_weight, self.force_weight, self.stress_weight, self.include_stress, self.train_local_mtps,pop)
                elif not self.train_on_every_ensemble:
                    if pop == start_pop:
                        ensemble_structures = self.minim.ensemble.structures
                        ase_structures_list_to_cfg(ensemble_structures,'preselected.cfg')
                        os.system('touch set.cfg')
                        train_mtp_on_cfg(self.specorder,self.mlip_run_command, self.path_to_mlip, self.pot_name, 
                            self.ab_initio_calculator, self.ab_initio_parameters, self.ab_initio_modules_load, self.ab_initio_kresol, self.ab_initio_pseudos, self.ab_initio_cluster, 
                            self.iteration_limit, self.energy_weight, self.force_weight, self.stress_weight, self.include_stress, self.train_local_mtps,pop)
                    else:
                        pass

                # Compute energies and forces
                self.minim.ensemble.compute_ensemble(self.calc, True, stress_numerical,
                                                 cluster = self.cluster)
                #self.minim.ensemble.get_energy_forces(self.calc, True, stress_numerical = stress_numerical)

                print("RELAX force length:", len(self.minim.ensemble.force_computed))
                
                if ensemble_loc is not None and self.save_ensemble:
                    self.minim.ensemble.save_bin(ensemble_loc, pop)
                print("RELAX force length:", len(self.minim.ensemble.force_computed))

            self.minim.population = pop
            self.minim.init(delete_previous_data = False)

            print("RELAX force length:", len(self.minim.ensemble.force_computed))
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
            Vol = self.minim.dyn.structure.get_volume()

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
                mark_gibbs = "<--"

            # Extract the bulk modulus from the cell minimization
            new_bulkmodulus = 1 / (3 * BFGS.alpha * self.minim.dyn.structure.get_volume())
            new_bulkmodulus *= sscha.SchaMinimizer.__evA3_to_GPa__

            # Print the enthalpic contribution
            message = """
 ======================
 ENTHALPIC CONTRIBUTION
 ======================

 P = {:.4f} GPa   V = {:.4f} A^3

 P V = {:.8e} eV

 Helmoltz Free energy = {:.10e} eV {}
 Gibbs Free energy = {:.10e} eV {}
 Zero energy = {:.10e} eV

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
            if self.fix_cell_shape:
                # Use a isotropic stress tensor to keep the same cell shape
                cell_gradient = I * (Press - target_press_evA3)
            else:
                cell_gradient = (stress_tensor - I *target_press_evA3)

            new_uc = self.minim.dyn.structure.unit_cell.copy()
            BFGS.UpdateCell(new_uc,  cell_gradient, fix_volume)

            # Strain the structure and the q points preserving the symmetries
            self.minim.dyn.AdjustToNewCell(new_uc)
            #self.minim.dyn.structure.change_unit_cell(new_uc)

            message = """
 Currently estimated bulk modulus = {:8.3f} GPa
 (Note: this is just indicative, do not use it for computing bulk modulus)

 """.format(new_bulkmodulus)
            print(message)


            print (" New unit cell:")
            print (" v1 [A] = (%16.8f %16.8f %16.8f)" % (new_uc[0,0], new_uc[0,1], new_uc[0,2]))
            print (" v2 [A] = (%16.8f %16.8f %16.8f)" % (new_uc[1,0], new_uc[1,1], new_uc[1,2]))
            print (" v3 [A] = (%16.8f %16.8f %16.8f)" % (new_uc[2,0], new_uc[2,1], new_uc[2,2]))

            print ()
            print ("Check the symmetries in the new cell:")
            sys.stdout.flush()
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

        self.start_pop = pop
        return (not running1) and (not running2)

def train_mtp_on_cfg(specorder, mlip_run_command, path_to_mlip, pot_name, 
                    ab_initio_calculator, ab_initio_parameters, ab_initio_modules_load,
                    ab_initio_kresol, ab_initio_pseudos, ab_initio_cluster,
                    iteration_limit = 500, 
                    energy_weight = 1.0, 
                    force_weight = 0.01,
                    stress_weight = 0.001,
                    include_stress = False,
                    train_local_mtps = False,
                    pop = 0,
                    np_ab_initio = 1,
                    np_mlp_train = 1):

    print(f"Preparing files for (re)training MTP {pot_name}")
    # cmd_gain_cfg   = 'cat preselected.cfg* >> preselected.cfg; rm -f preselected.cfg.*'
    cmd_gain_cfg   = 'cat preselected.cfg* >> preselected.cfg'
    cmd_select_add = f'{mlip_run_command} {path_to_mlip} \
                        select_add {pot_name} set.cfg preselected.cfg selected.cfg'
    cmd_mlip_train = f'{mlip_run_command} {path_to_mlip} \
                        train {pot_name} set.cfg \
                        --tolerance=0.01 \
                        --energy_weight={energy_weight} \
                        --force_weight={force_weight} \
                        --stress_weight={stress_weight} \
                        --iteration_limit={iteration_limit}'

    if train_local_mtps:
        try:
            os.system(f'cp {pot_name}.bak {pot_name}')
        except:
            pass


    try:
        os.system(cmd_gain_cfg)
    except:
        pass
    os.system(cmd_select_add)
    # n_cfg_cmd = 'ls preselected.cfg.* | wc -l'
    n_cfg_cmd = 'grep "BEGIN" selected.cfg | wc -l'
    # n_cfg_cmd = 'grep "BEGIN" preselected.cfg | wc -l'
    n_cfg = int(os.popen(n_cfg_cmd).read().split()[0])
    print(f"There are {n_cfg} configurations that need to be added to the training set") 

    try:
        split_cfg('selected.cfg')
    except FileNotFoundError:
        print('There are no new structures that need to be added to the train set!')
        return

    for i in range(n_cfg):
        print(f"Calculating ab initio energies and forces for configuration selected.cfg.{i}")
        # cmd_convert  = f'{mlip_run_command} {path_to_mlip} convert --output_format=poscar sampled.cfg.{i} {i}.POSCAR' 
        # os.system(cmd_convert)
        ab_initio_dir = f'ab_initio_dir_{i}'+ '_' + str(datetime.now()).replace(' ','_').replace(':','-')
        cmd_mkdir_abinitio = f'mkdir {ab_initio_dir}'
        if not os.path.exists(ab_initio_dir): 
            os.system(cmd_mkdir_abinitio)
        else: 
            print(f'Ab initio directory {ab_initio_dir} already exists!')
        
        ab_initio_ase_atoms = one_cfg_to_atoms(f'selected.cfg.{i}',specorder) 

        if ab_initio_calculator == "QE":
            
            k_points = calc_ngkpt(ab_initio_ase_atoms.cell.reciprocal(),ab_initio_kresol)
            print(f'KPOINTS: {k_points}')
            if ab_initio_modules_load != None:
                ab_initio_run_command = f"module load {ab_initio_modules_load}; mpirun -np {np_ab_initio} pw.x -in ESP.pwi > ESP.pwo"
            else:
                ab_initio_run_command = f"mpirun -np {np_ab_initio} pw.x -in ESP.pwi > ESP.pwo" 
            ab_initio_cc_calc = Espresso(pseudopotentials = ab_initio_pseudos,
                                        input_data = ab_initio_parameters,
                                        kpts = k_points,
                                        command = ab_initio_run_command,
                                        koffset  = (0,0,0))
            ab_initio_cc_calc.set_directory(ab_initio_dir)
            in_ext  = ".pwi"
            out_ext = ".pwo"   
            print(ab_initio_cc_calc.command)

        ### Starting Ab initio calculations ###
        if ab_initio_cluster != None:
            print("Using HPC for ab initio calculations")
            if ab_initio_calculator == "QE":
                # Here we need to improve run_atoms method to allow it to wait 
                ab_initio_cluster.run_atoms(ab_initio_cc_calc, ab_initio_ase_atoms, in_extension = in_ext, out_extension = out_ext)

            # my_hpc.read_results(cc_calc, my_hpc.label)
        elif ab_initio_cluster == None:

            if ab_initio_calculator == "QE":
                cc_struct = CC.Structure.Structure()
                cc_struct.generate_from_ase_atoms(ab_initio_ase_atoms)
                ab_initio_cc_calc.set_label("ESP")
                ab_initio_cc_calc.calculate(cc_struct)

        ### Ending Ab initio calculations ###

        # print(dir(cc_calc.structure))
        if ab_initio_calculator == "QE":
            # using finction for CellConstructor calculator
            cc_calc_to_cfg(ab_initio_cc_calc,f'input.cfg.{i}',specorder, include_stress)

    os.system(f'cp set.cfg set.cfg.bak')

    os.system('cat input.cfg* >> set.cfg')
    os.system(f'cp {pot_name} {pot_name}.bak')
    print('Start training MLIP')
    os.system(cmd_mlip_train)
    print('End training MLIP')
    if train_local_mtps:
        os.system(f'cat input.cfg* >> all_input_pop_{str(pop)}.cfg')
        os.system(f'cat selected.cfg* >> all_selected_pop_{str(pop)}.cfg')
        os.system(f'mv set.cfg set_pop_{str(pop)}.cfg')
        os.system(f'cp {pot_name} {pot_name}.pop_{str(pop)}')
        # os.system(f'cp {pot_name}.bak {pot_name}')
        os.system('rm -f input.cfg*')
        os.system('rm -f selected.cfg*')
        os.system('rm -f preselected.cfg*')
    else:
        os.system(f'cat input.cfg* >> all_input.cfg')
        os.system(f'cat selected.cfg* >> all_selected.cfg')
        os.system('rm -f input.cfg*')
        os.system('rm -f selected.cfg*')
        os.system('rm -f preselected.cfg*')

    return  

def ase_structures_list_to_cfg(ase_structures_list,path_to_cfg):

    """
    Function for converting list of ase stuctures to one cfg file for MLIP package.

    ase_structures_list - list of ase structures
    path_to_cfg - path where the cfg file will be saved

    """

    with open(path_to_cfg, 'w') as f:
        for ase_structure in ase_structures_list:
            try:
                # if structure is ase structure
                structure = CC.Structure.Structure()
                structure.generate_from_ase_atoms(ase_structure)
            except AttributeError:
                # if structure is CC structure
                structure = ase_structure
                
            # f.write('\n')
            f.write('BEGIN_CFG\n')
            f.write(' Size\n')
            f.write(f'    {structure.N_atoms}\n')
            f.write(' Supercell\n')
            cell = structure.unit_cell
            f.write(f'         {cell[0][0]:.6f}      {cell[0][1]:.6f}      {cell[0][2]:.6f}\n')
            f.write(f'         {cell[1][0]:.6f}      {cell[1][1]:.6f}      {cell[1][2]:.6f}\n')
            f.write(f'         {cell[2][0]:.6f}      {cell[2][1]:.6f}      {cell[2][2]:.6f}\n')
            f.write(' AtomData:  id type       cartes_x      cartes_y      cartes_z\n')
            # coords = structure.get_xcoords()
            coords = structure.coords
            typat  = structure.get_ityp()
            for i in range(structure.N_atoms):
                f.write(f'             {i+1:3}    {typat[i]}       {coords[i][0]:.6f}      {coords[i][1]:.6f}      {coords[i][2]:.6f} \n')
            f.write(f'END_CFG\n')
            f.write(f'\n')


    return

def calc_ngkpt(recip, kspacing):
    to_ang_local = 1
    
    N_from_kspacing = []
    for i in 0, 1, 2:
        N_from_kspacing.append( int(np.ceil( (np.linalg.norm(recip[i]) / to_ang_local) / kspacing)) )

    return N_from_kspacing

def one_cfg_to_atoms(path_to_cfg,specorder):
    """
    Function for converting cfg file with ONE structure to ase_atoms object

    path_to_cfg - a path to the cfg file th
    specorder - a list with species order

    """

    with open(path_to_cfg,'r') as f:
        lines = f.readlines()
        num_at = int(lines[2].split('/n')[0])
        vec1 = [float(lines[4].split()[0]), float(lines[4].split()[1]), float(lines[4].split()[2])]
        vec2 = [float(lines[5].split()[0]), float(lines[5].split()[1]), float(lines[5].split()[2])]
        vec3 = [float(lines[6].split()[0]), float(lines[6].split()[1]), float(lines[6].split()[2])]
        xcart = []
        symbols = []
        for i in range(num_at):
            xcart.append([])
            t = int(lines[8+i].split()[1])
            x = float(lines[8+i].split()[2])
            y = float(lines[8+i].split()[3])
            z = float(lines[8+i].split()[4])
            xcart[i].append(x)
            xcart[i].append(y)
            xcart[i].append(z)
            symbols.append(specorder[t])

        cell = [vec1,vec2,vec3]

    # print(f'num_at = {num_at}')
    # print(f'cell = {cell}')
    # print(f'xcart = {xcart}')
    # print(f'symbols = {symbols}')

    ase_atoms_cfg = ase.Atoms(symbols = symbols,
                              positions = xcart,
                              cell = cell,
                            #   atoms = symbols,
                             )

    return ase_atoms_cfg

def split_cfg(path_to_cfg):
    """
    Function for splitting one cfg with multiple configurations to several cfg's with one configuration in each
    """
    n_config = 0
    configs = []
    with open(path_to_cfg, 'r') as f:
        lines = f.readlines()
        for line in lines:
        # while True:
            # line = f.readline()
            if 'BEGIN_CFG' in line:
                configs.append([])
                configs[n_config].append(line)
            elif 'END_CFG' in line:
                configs[n_config].append(line)
                n_config += 1
            else:
                try:
                    configs[n_config].append(line)
                except IndexError:
                    continue

    for i, config in enumerate(configs):
        path_to_cfg_i = f'{path_to_cfg}.{i}'
        with open(path_to_cfg_i, 'w') as f:
            f.writelines(config)

    return
                

def cc_calc_to_cfg(cc_calc,path_to_cfg,specorder,include_stress = False):
    """
    Function for converting cellconstructor calculator object to cfg file for MTP training.
    The cfg is written 

    cc_calc - cellconstructor calculator object
    path_to_cfg - path where the cfg file will be saved
    specorder - a list with species order

    """

    with open(path_to_cfg, 'w') as f:
        f.write('\n')
        f.write('BEGIN_CFG\n')
        f.write(' Size\n')
        f.write(f'    {cc_calc.structure.N_atoms}\n')
        f.write(' Supercell\n')
        cell = cc_calc.structure.unit_cell
        f.write(f'         {cell[0][0]:.6f}      {cell[0][1]:.6f}      {cell[0][2]:.6f}\n')
        f.write(f'         {cell[1][0]:.6f}      {cell[1][1]:.6f}      {cell[1][2]:.6f}\n')
        f.write(f'         {cell[2][0]:.6f}      {cell[2][1]:.6f}      {cell[2][2]:.6f}\n')
        f.write(' AtomData:  id type       cartes_x      cartes_y      cartes_z           fx          fy          fz\n')
        # coords = cc_calc.structure.get_xcoords()
        coords = cc_calc.structure.coords

        # count = 1
        # for i, atm in enumerate(structure.atoms):
        #     if not atm in dictionary:
        #         dictionary[atm] = count
        #         count += 1
        
        # ityp = [dictionary[x] for x in self.atoms]

        typat  = cc_calc.structure.get_atomic_types()

        forces = cc_calc.results['forces']
        for i in range(cc_calc.structure.N_atoms):
            f.write(f'             {i+1:3}    {typat[i]-1}       {coords[i][0]:.6f}      {coords[i][1]:.6f}      {coords[i][2]:.6f}      {forces[i][0]:.6f} {forces[i][1]:.6f} {forces[i][2]:.6f}\n')
        f.write(f' Energy\n')
        f.write(f' {cc_calc.results["energy"]}\n')
        if "stress" in cc_calc.results and include_stress:
            pxx = -cc_calc.results["stress"][0]*cc_calc.structure.get_volume()
            pyy = -cc_calc.results["stress"][1]*cc_calc.structure.get_volume()
            pzz = -cc_calc.results["stress"][2]*cc_calc.structure.get_volume()
            pyz = -cc_calc.results["stress"][3]*cc_calc.structure.get_volume()
            pxz = -cc_calc.results["stress"][4]*cc_calc.structure.get_volume()
            pxy = -cc_calc.results["stress"][5]*cc_calc.structure.get_volume()
            f.write(f' PlusStress:  xx          yy          zz          yz          xz          xy\n')
            f.write(f'        {pxx}    {pyy}    {pzz}    {pyz}    {pxz}    {pxy}\n')
        f.write(f'END_CFG\n')
        f.write(f'\n')


    return

def ase_calc_to_cfg(ase_calc,path_to_cfg,specorder):
    return