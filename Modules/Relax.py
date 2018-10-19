# -*- coding: utf-8 -*-

"""
This module performs the relax over more
SCHA minimization. It can be both a constant temperature and a 
constant pressure relax. 
"""
import numpy as np
import sscha, sscha.Ensemble, sscha.SchaMinimizer

class SSCHA:
    def __init__(self, minimizer, ase_calculator, N_configs, max_pop = 20):
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
        """
        
        self.minim = minimizer
        self.calc = ase_calculator
        self.N_configs = N_configs
        self.max_pop = max_pop
        
        self.__cfpre__ = None
        self.__cfpost__ = None
        self.__cfg__ = None
        
    def setup_custom_functions(self, custom_function_pre = None,
                               custom_function_post = None,
                               custom_function_gradient = None):
        """
        This subroutine setup which custom functions should be called during the minimization.
        Look for the SCHA_Minimizer.run() method for other details.
        """
        
        self.__cfpre__ = custom_function_pre
        self.__cfpost__ = custom_function_post
        self.__cfg__ = custom_function_gradient
        
        
    def relax(self, restart_from_ens = False, get_stress = False,
              ensemble_loc = ".", start_pop = 1):
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
                Where the ensemble of each population is saved on the disk. You can specify None
                if you do not want to save the ensemble (useful to avoid disk I/O for force fields)
            start_pop : int, optional
                The starting index for the population, used only for saving the ensemble and the dynamical 
                matrix.
            
        Returns
        -------
            status : bool
                True if the minimization converged, False if the maximum number of 
                populations has been reached.
        """
        
        pop = start_pop
                
        running = True
        while running:
            # Generate the ensemble
            self.minim.ensemble.generate(self.N_configs)
            
            # Compute energies and forces
            self.minim.ensemble.get_energy_forces(self.calc, get_stress)
            
            if ensemble_loc is not None:
                ens.save_bin(ensemble_loc, pop)
            
            self.minim.population = pop
            self.minim.init()
        
            self.minim.run(custom_function_pre = self.__cfpre__,
                           custom_function_post = self.__cfpost__,
                           custom_function_gradient = self.__cfg__)
        
            
            self.minim.finalize()
        
            # Save the dynamical matrix
            self.minim.dyn.save_qe("dyn_pop%d_" % pop)
        
            # Check if it is converged
            running = not self.minim.is_converged()
            pop += 1
            
            
            if pop > self.max_pop:
                running = False
                
        return self.minim.is_converged()
    
    
    def vc_relax(self, target_press = 0, 
                 restart_from_ens = False,
                 ensemble_loc = ".", start_pop = 1):
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
            restart_from_ens : bool, optional
                If True the ensemble is used to start the first population, without recomputing
                energies and forces. If False (default) the first ensemble is overwritten with
                a new one, and the minimization starts.
            get_stress : bool, optional
                If true the stress tensor is calculated. This may increase the computational
                cost, as it will be computed for each ab-initio configuration (it may be not available
                with some ase calculator)
            ensemble_loc : string
                Where the ensemble of each population is saved on the disk. You can specify None
                if you do not want to save the ensemble (useful to avoid disk I/O for force fields)
            start_pop : int, optional
                The starting index for the population, used only for saving the ensemble and the dynamical 
                matrix.
            
        Returns
        -------
            status : bool
                True if the minimization converged, False if the maximum number of 
                populations has been reached.
        """
        
        
        
        pop = start_pop
                
        running = True
        while running:
            # Generate the ensemble
            self.minim.ensemble.generate(self.N_configs)
            
            # Compute energies and forces
            self.minim.ensemble.get_energy_forces(self.calc, True)
            
            if ensemble_loc is not None:
                ens.save_bin(ensemble_loc, pop)
            
            self.minim.population = pop
            self.minim.init()
        
            self.minim.run(custom_function_pre = self.__cfpre__,
                           custom_function_post = self.__cfpost__,
                           custom_function_gradient = self.__cfg__)
        
            
            self.minim.finalize()
            
            # Get the stress tensor [ev/A^3]
            stress_tensor = self.minim.get_stress_tensor() * sscha.SchaMinimizer.__RyBohr3_to_evA3__
            
            # Get the pressure
            Press = np.trace(stress_tensor) / 3
            
            # Get the volume
            Vol = np.linalg.det(self.minim.dyn.structure.unit_cell)
            
            # Get the Gibbs free energy
            gibbs = self.minim.get_free_energy() * sscha.SchaMinimizer.__RyToev__ - Press * Vol - self.minim.eq_energy
            
            
            # Print the enthalpic contribution
            print ""
            print " ENTHALPIC CONTRIBUTION "
            print " ====================== "
            print ""
            print "  P V = %.8e eV " % (Press * Vol)
            print ""
            print " Gibbs Free energy = %.8e eV " % gibbs
            print " (Zero energy = %.8e eV) " % self.minim.eq_energy
        
        
            # Save the dynamical matrix
            self.minim.dyn.save_qe("dyn_pop%d_" % pop)
        
            # Check if it is converged
            running1 = not self.minim.is_converged()
            pop += 1
            
            
            if pop > self.max_pop:
                running = False
                
        return self.minim.is_converged()