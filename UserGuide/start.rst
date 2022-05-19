Quick start
===========

In this chapter we provide ready to use examples to setup your first SSCHA calculation.

The free energy of gold: a simulation in the NVT ensemble
---------------------------------------------------------

This simple tutorial explains how to setup a SSCHA calculation starting just from the structure, in this case a cif file we downloaded from the [Materials Project](https://materialsproject.org/materials/mp-81/) database.

You can find there a lot of structures ready to use for your SSCHA runs.

For the purpouse of this tutorial, we are going to use the EMT force field, so that the calculation can be run in a laptop without the need of a supercomputer.
We explain in a later section how to couple the SSCHA with a cluster to submit the same calculation fully ab-initio.

Starting from the Gold structure in the primitive cell, to run the SSCHA we need:
 - Compute the harmonic phonons (dynamical matrix)
 - Remove imaginary frequencies (if any)
 - Run the SSCHA

We prepared an input file in the form of a python script (tested with python-sscha version 1.2) which makes all these passages automatically.

.. code-block:: python

   # Import the sscha code
   import sscha, sscha.Ensemble, sscha.SchaMinimizer, sscha.Relax, sscha.Utilities

   # Import the cellconstructor library to manage phonons
   import cellconstructor as CC, cellconstructor.Phonons
   import cellconstructor.Structure, cellconstructor.calculators

   # Import the force field of Gold
   import ase, ase.calculators
   from ase.calculators.emt import EMT

   # Import numerical and general pourpouse libraries
   import numpy as np, matplotlib.pyplot as plt
   import sys, os


   """
   Here we load the primitive cell of Gold from a cif file.
   And we use CellConstructor to compute phonons from finite differences.
   The phonons are computed on a q-mesh 4x4x4
   """

   gold_structure = CC.Structure.Structure()
   gold_structure.read_generic_file("Au.cif")

   # Get the force field for gold
   calculator = EMT()

   # Relax the gold structure (useless since for symmetries it is already relaxed)
   relax = CC.calculators.Relax(gold_structure, calculator)
   gold_structure_relaxed = relax.static_relax()

   # Compute the harmonic phonons
   # NOTE: if the code is run with mpirun, the calculation goes in parallel
   gold_harmonic_dyn = CC.Phonons.compute_phonons_finite_displacements(gold_structure_relaxed, calculator, supercell = (4,4,4))

   # Impose the symmetries and 
   # save the dynamical matrix in the quantum espresso format
   gold_harmonic_dyn.Symmetrize()
   gold_harmonic_dyn.save_qe("harmonic_dyn")


   # If the dynamical matrix has imaginary frequencies, remove them
   gold_harmonic_dyn.ForcePositiveDefinite()

   """
   gold_harmonic_dyn is ready to start the SSCHA calculation.

   Now let us initialize the ensemble, and the calculation at 300 K.
   We will run a NVT calculation, using 100 configurations at each step
   """

   TEMPERATURE = 300
   N_CONFIGS = 50
   MAX_ITERATIONS = 20

   # Initialize the random ionic ensemble
   ensemble = sscha.Ensemble.Ensemble(gold_harmonic_dyn, TEMPERATURE)

   # Initialize the free energy minimizer
   minim = sscha.SchaMinimizer.SSCHA_Minimizer(ensemble)
   minim.set_minimization_step(0.01) 

   # Initialize the NVT simulation
   relax = sscha.Relax.SSCHA(minim, calculator, N_configs = N_CONFIGS,
   max_pop = MAX_ITERATIONS)

   # Define the I/O operations
   # To save info about the free energy minimization after each step
   ioinfo = sscha.Utilities.IOInfo()
   ioinfo.SetupSaving("minim_info")
   relax.setup_custom_functions(custom_function_post = ioinfo.CFP_SaveAll)


   # Run the NVT simulation (save the stress to compute the pressure)
   relax.relax(get_stress = True)

   # If instead you want to run a NPT simulation, use
   # The target pressure is given in GPa.
   #relax.vc_relax(target_press = 0)

   # You can also run a mixed simulation (NVT) but with variable lattice parameters
   #relax.vc_relax(fix_volume = True)

   # Now we can save the final dynamical matrix
   # And print in stdout the info about the minimization
   relax.minim.finalize()
   relax.minim.dyn.save_qe("sscha_T{}_dyn".format(TEMPERATURE))


Now save the file as `sscha_gold.py` and execute it with:

.. code-block:: bash

   $ python sscha_gold.py > output.log

And that's it.
**Congratulations!**


While the input may seem long, it is heavily commented, but lets go through it step by step.
At the very beginning, we simply import the sscha libraries, cellconstructor, the math libraries and the force field. This is done in python with the `import` statemets. 

The first real part of the code is:

.. code-block:: python
   
   gold_structure = CC.Structure.Structure()
   gold_structure.read_generic_file("Au.cif")

   # Get the force field for gold
   calculator = EMT()

   # Relax the gold structure (useless since for symmetries it is already relaxed)
   relax = CC.calculators.Relax(gold_structure, calculator)
   gold_structure_relaxed = relax.static_relax()

Here we initialize a cellconstructor structure from the cif file downloaded from the material database (*Au.cif*). We initialize the EMT calculator from ASE, and relax the structure.

In the case of Gold the relaxation is useless, as it is a FCC structure with Fm-3m symmetry group and 1 atom per primitive cell. This means the atomic positions have no degrees of freedom, thus the relaxation will end before even start.

In the next part of the code, we perform the harmonic phonon calculation using cellconstructor and a finite displacement approach:

.. code-block:: python

   gold_harmonic_dyn = CC.Phonons.compute_phonons_finite_displacements(gold_structure_relaxed, calculator, supercell = (4,4,4))

   # Impose the symmetries and 
   # save the dynamical matrix in the quantum espresso format
   gold_harmonic_dyn.Symmetrize()
   gold_harmonic_dyn.save_qe("harmonic_dyn")


   # If the dynamical matrix has imaginary frequencies, remove them
   gold_harmonic_dyn.ForcePositiveDefinite()


The method `compute_phonons_finite_displacements` is documented in the CellConstructor guide. It requires the structure (in this case `gold_structure_relaxed`), the force-field (`calculator`) and the supercell for the calculation. In this case we use a 4x4x4 (equivalent to 64 atoms). This may not be sufficient to converge all the properties, especially at very high temperature, but it is just a start.

Note that  `compute_phonons_finite_displacements`  works in parallel with MPI, therefore, if the script is executed with `mpirun -np 16 python myscript.py` it will split the calculations of the finite displacements across 16 processors. You need to have mpi4py installed.

After computing the harmonic phonons in gold_harmonic_dyn, we impose the correct symmetrization and the acousitic sum rule with the `Symmetrize` method, and save the result in the quantum ESPRESSO format with `save_qe`.
This should not be the case for Gold, however, if we have a structure which has imaginary phonon frequencies, we need to get rid of them before starting the SSCHA. This is achieved with `ForcePositiveDefinite` (see CellConstructor documentation for more details on how these methods work).


**Now we are ready to submit the SSCHA calculation in the NVT ensemble!**.
The important parameters are:
  - The temperature
  - The number of random configurations in the ensemble
  - The maximum number of iterations


These parameters are almost self-explaining. However, we give a brief overview of how the SSCHA works to help you understand which are the best one for your case.
While MD or MC calculation represent the equilibrium probability distribution over time of the system by updating a single structure, the SSCHA encodes the whole probability distribution as an analytical function. Therefore, to compute properties, we can generate on the fly the ionic configurations that represent the equilibrium distributions.
The number of random configuration is exactly how many ionic configuration we generate to compute the properties (Free energy and Stress tensors)

The code that sets up and perform the SSCHA is the following:

.. code-block:: python

   TEMPERATURE = 300
   N_CONFIGS = 50
   MAX_ITERATIONS = 20

   # Initialize the random ionic ensemble
   ensemble = sscha.Ensemble.Ensemble(gold_harmonic_dyn, TEMPERATURE)

   # Initialize the free energy minimizer
   minim = sscha.SchaMinimizer.SSCHA_Minimizer(ensemble)
   minim.set_minimization_step(0.01) 

   # Initialize the NVT simulation
   relax = sscha.Relax.SSCHA(minim, calculator, N_configs = N_CONFIGS,
   max_pop = MAX_ITERATIONS)

   # Define the I/O operations
   # To save info about the free energy minimization after each step
   ioinfo = sscha.Utilities.IOInfo()
   ioinfo.SetupSaving("minim_info")
   relax.setup_custom_functions(custom_function_post = ioinfo.CFP_SaveAll)

   # Run the NVT simulation
   relax.relax(get_stress = True)
   



So you see many classes. `ensemble` represent the ensemble of ionic configurations. We initialize it with the dynamical matrix (which represent how much atoms fluctuate around the centroids) and the temperature.
`minim` is a `SSCHA_Minimizer` object, which performs the free energy minimization. It contains all the info regarding the minimization algorithm, as the initial timestep (that here we set to 0.01). You can avoid setting the time-step, as the code will automatically guess the best value.
The `relax` is a `SSCHA` object: the class that takes care about the simulation and automatizes all the steps to perform a NVT or NPT calculation.
We pass the minimizer (which contains the ensemble with the temperature), the force-field (`calculator`), the number of configurations `N_configs` and the maximum number of iterations.

In this example, most of the time is spent in the minimization, however, if we replace the force-field with ab-initio DFT, the time tu run the minimization is negligible with respect to the time to compute energies and forces on the ensemble configurations.
The total (maximum) number of energy/forces calculations is equal to the number of configurations times the number of iterations (passed through the `max_pop` argument).

The calculation is submitted with `relax.relax()`. However, before running the calculation we introduce another object, the `IOInfo`.
This tells the `relax` to save information of the free energy, its gradient and the anharmonic phonon frequencies during the minimization in the files *minim_info.dat° and *minim_info.freqs*. It is not mandatory to introduce them, but it is very usefull as it allows to visualize the minimization while it is running.


















To quickly start using the code, we recommend using the jupyter notebooks with examples we provide in the Tutorials directory of the source code.

Tutorials are organized as follows:

1. Setup from the structure and manual submission: PbTe tutorial. Here you learn how to set up a SSCHA calculation starting just with the structure (we provide a .cif file of the PbTe at high temperature). The tutorial will guide you step by step. You will learn how to: prepare the starting data needed for the SSCHA calculation, generate a random ensemble, save the ensemble and prepare input files for your favorite ab-initio code, read back the energies and the forces inside SSCHA, run a SSCHA minimization. You will also learn how to use ASE and the Cluster module to automatize the calculation of the ensemble and submit it to a HPC system.
2. Automatic relaxation with a force field: SnTe_ToyModel. Here, we show how to use a force-field for a SSCHA calculation, running everything on your computer. We also will explain how to calculate the free energy hessian for second-order phase transitions, and study a phase transition as a function of temperature.
3. Variable cell relaxation: LaH10 tutorial. Here you learn how to perform an automatic calculation with a variable cell. You will exploit the quantum effect to search the high-temperature superconductive phase (Fm-3m) of LaH10 below 200 GPa, starting from a distorted structure. 
4. Hessian matrix calculation for second-order phase transitions: H3S tutorial. Here you reproduce the calculation of the Hessian of the free energy to assert the stability of the H3S phase.
5. Spectral properties: Spectral_Properties. In this tutorial, we explain how to use the post-processing utilities of the SSCHA to calculate the phonon spectral function, and computing phonon lifetimes, and plotting interacting phonon dispersion. We provide an ensemble for PbTe already computed ab-initio.


The jupyter notebooks are interactive, to quickly start with your simulation, pick the tutorial that resembles the kind of calculation you want to run, and simply edit it directly in the notebook. 

