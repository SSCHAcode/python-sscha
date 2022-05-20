Frequently Asked Questions (FAQs)
=================================

.. role:: raw-html(raw)
    :format: html


Here we answer to most common question we received.

Setup the calculation
---------------------

How do I start a calculation if the Dynamical matrices have imaginary frequencies?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A good starting point for a sscha minimization are the dynamical matrix obtained from a harmonic calculation. However, they can have imaginary frequencies. This may be related to both instabilities (the structure is a saddle-point of the Born-Oppenheimer energy landscape) or to a not well-converged choice of the parameters for computing the harmonic frequencies.
In both cases, it is very easy to get a new dynamical matrix that is positive definite and can be used as a starting point. An example is made in Turorial on H3S.
Assuming your not positive definite dynamical matrix is in Quantum Espresso format "harm1" ... "harmN" (with N the number of irreducible q points), you can generate a positive definite dynamical matrix "positive1" ... "positiveN" with the following python script that uses CellConstructor.

.. code-block:: python

   # Load the cellconstructor library
   import cellconstructor as CC
   import cellconstructor.Phonons

   # Load the harmonic not-positive definite dynamical matrix
   # We are reading 6 dynamical matrices
   harm = CC.Phonons.Phonons("harm", nqirr = 6) 

   # Apply the acoustic sum rule and the symmetries
   harm.Symmetrize() 

   # Force the frequencies to be positive definite
   harm.ForcePositiveDefinite() 

   # Save the final dynamical matrix, ready to be used in a sscha run
   harm.save_qe("positive")
       

The previous script (that we can save into *script.py*) will generate the positive definite matrix ready for the sscha run. It may be executed with

.. code-block:: console
		    
   $ python script.py

  

What are the reasonable values for the steps?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Starting from version 1.2, the line minimization is implemented. This means that there is no need to specify the value of the minimization step as the code will automatically find it.

However, if the code takes too long to get a good timestep at the begining of a calculation (especially at the very first iteration or if few configurations are employed), you could speedup the calculation providing a smaller initial guess than the default one (1).
This is done in the python script by calling the function

.. code-block:: python
   
   minim.set_minimization_step(0.1)


Where minim is the sscha.SchaMinimizer.SSCHA_Minimizer class. You can select the step also in the namespace input by setting the following variables in the inputscha namespace

**lambda_w** is the step in the atomic positions (stand-alone program input).
   
**lambda_a** is the step in the dynamical matrix (stand-alone program input).


   
In a NPT or NVT with variable lattice, what is a reasonable value for the bulk modulus?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The bulk modulus is just an indicative parameter used to guess the optimal step of the lattice parameters to converge as quickly as possible.
It is expressed in GPa. You can find online the bulk modulus for many materials. Find a material similar to the one you are studying and look if there is in literature a bulk modulus.

The default value is good for most case (equal to 100), but it could be too low for very hard materials (like diamond, which is 500 GPa, or high-pressure stuff). If you are not sure, it is safer to choose an higher value of the bulk modulus, as the code is going to optimize it during the simulation anyway.

If you have no idea on the bulk modulus, you can easily compute them by doing two static *ab initio* calculations at very close volumes (by varying the cell size), and then computing the differences between the pressure:

.. math::

   B = - \Omega \frac{dP}{d\Omega}

where :math:`\Omega` is the unit-cell volume and :math:`P` is the pressure (in GPa).



It is always good to run NVT before any NPT simulation?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In general, it is good to have a reasonable dynamical matrix before starting with a relaxation with variable cell (vc_relax).
Therefore, to avoid mooving the volume upward and backward, always start with a NVT simulation with fixed lattice (the relax method of SSCHA class) and then run a NPT or a NVT with variable lattice (vc_relax method), starting from the static lattice solution.

    
How may I run a calculation neglecting symmetries?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can tell the code to neglect symmetries with the :code:`neglect_symmetries = .true.` flag.
In the python script, this is done setting the attribute *neglect_symmetries* of sscha.SchaMinimizer.SSCHA_Minimizer to False.



In which units are the lattice vectors, the atomic positions, and the mass of the atoms in the dynamical matrix file?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The dynamical matrix follows the quantum espresso units. They are Rydberg atomic units (unit of mass is 1/2  the electron mass, energy is Ry, positions are in Bohr. However, espresso may have an ibrav not equal to zero (the third number in the header of the dynamical matrix). In this case, please, refer to the espresso ibrav guide in the `PW.x input description <https://www.quantum-espresso.org/Doc/INPUT_PW.html#idm199>`
  

What is the difference between different kinds of minimization (preconditioning and root_representation)?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You do not need to worry about these parameters, as starting from version 1.2 the code has a robust workflow that should avoid bothering you with these details.
However, if you are curious and want to know a bit more on the details here it is the explanation:
We provide three different advanced algorithms called in **root_representation**, that can be either **normal**, or **sqrt**, or **root4** (inside &inputscha namespace or the SSCHA_Minimizer object)
In this way, instead of minimizing the :math:`\Phi` matrix, we minimize with respect to :math:`\sqrt{\Phi}` or :math:`\sqrt[4]{\Phi}`.
Therefore the new dynamical matrix is constrained in a space that is positive definite. Moreover, it has been proved that :math:`\sqrt[4]{\Phi}`
minimization has a better condition number than the original one and thus should reach the minimum faster.

Alternatively, a similar effect to the speedup in the minimization obtained with **root4** is possible using the preconditioning (by setting **preconditioning** or **precond_dyn** to True in the input file or the python script, respectively). This way also the single minimization step runs faster, as it avoids passing in the root space of the dynamical matrix (but indeed, you can have imaginary frequencies).

Since the gradient computation is much slower (especially for a system with more than 80 atoms in the supercell) without the preconditioning,
it is possible to combine the preconditioning with the root representation to have a faster gradient computation and to be guaranteed that
the dynamical matrix is positive definite by construction at each step.
However, in this way the good condition number obtained by the preconditioning (or the root4 representation) is spoiled. For this reason, when using the preconditioning, avoid using **root4**, and chose instead **sqrt** as root_representation.

The default values are:

.. code-block:: console
		    
   &inputscha
      root_representation = "normal"
      preconditioning = .true.
   &end

or in python
    
.. code-block:: python
		    
   # The ensemble has been loaded as ens
   minim = sscha.SchaMinimizer.SSCHA_Minimizer(ens)
   minim.root_representation = "normal"
   minim.precond_dyn = True



How do I fix the random number generator seed to make a calculation reproducible?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
As for version 1.2, this can be achieved only by using the python script.
Since python uses NumPy for random numbers generation, you can, at the beginning of the script that generates the ensemble, use the following:

.. code-block:: python
		    
   import numpy as np

   X = 0
   np.random.seed(seed = X)

where :code:`X` is the integer used as a seed. By default, if not specified, it is initialized with None that it is equivalent to initializing with the current local time.

   

On error and convergence of the free energy minimization
--------------------------------------------------------

    
The code stops saying it has found imaginary frequencies, how do I fix it?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
**Update python-sscha to version 1.2 (at least)!** This should be fixed.

If you do not want to update the code, set

.. code-block:: python

   minim.root_representation = 'root2'
       
This way the minimization strategy changes and it is mathematically impossible to get imaginary frequencies.
The same option can be activated within the namespace input

.. code-block:: console

   &inputscha
      root_representation = 'root2'
   &end

       
Why the gradient sometimes increases during a minimization?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Noting in principle assures that a gradient should always go down. It is possible at the beginning of the calculation when we are far from the solution that one of the gradients increases.
However, when we get closer to the solution, indeed the gradient must decrease.
If this does not happen it could be due to the ensemble that has fewer configurations than necessary. In this case, the good choice is to increase the number of effective sample size (the Kong-Liu ratio), to stop the minimization when the gradient starts increasing, or to increase the number of configurations in the ensemble.

In any case, what must decrease is free energy. If you see that the gradient is increasing but the free energy decreases, then the minimization is correct. However, if both the gradient and free energy are increasing, something is wrong, and you may require more configurations in each iteration.
This is especially true for system with few symmetries (or big primitive cells).

    

How do I check if my calculations are well converged?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
In general, if the gradient goes to zero and the Kong Liu ratio is above 0.5 probably your calculation converged very well. This means that when your calculation stops because it converged (not because it runs out of iterations), then it should be well converged.

There are some cases (especially in systems with many atoms) in which it is difficult to have an ensemble sufficiently big to reach this condition.
In these cases, you can look at the history of the frequencies in the last populations (there is a drift or random fluctuations?)


What is the final error on the structure or the dynamical matrix of a SCHA minimization?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To test the error, you can split the ensemble in two half and repeat the last minimization.
Then check at the difference between the result to have a rough estimation of the fluctuations.

To split the ensemble, refer to the FAQ :ref:`FAQ split`.


How do I understand if the free energy hessian calculation is converged?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The free energy hessian requires much more configurations than the SCHA minimization. First of all, to run the free energy Hessian, the SSCHA minimization must end with a gradient that can be decreased indefinitely without decreasing the KL below 0.7 /0.8.
Then you can estimate the error by repeating the hessian calculation with half of the ensemble and check how the frequencies of the hessian changes. This is also a good check for the final error on the frequencies.
    
You can split your ensemble in two by using the split function.

To split the ensemble, refer to the FAQ :ref:`FAQ split`.


.. _FAQ split:

How do I split the ensemble?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After you load or compute an ensemble you can split it and select only a portion of it to run the code.

.. code-block:: python

   # Assuming you loaded or computed the ensemble inside
   # ensemble

   # Let us create a mask that selects only the first half of the ensemble
   first_half_mask = np.zeros(ensemble.N, dtype = bool)
   first_half_mask[:ensemble.N//2] = True

   # Now we pass the mask to the ensemble to extract a new one
   # Containing only the configurations that correspond to the True
   # values of the mask
   first_half_ensemble = ensemble.split(first_half_mask)


After this code, the varialbe first_half_ensemble is a sscha.Ensemble.Ensemble that
can be used for any caluclation. 



How can I add more configurations to an existing ensemble?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
You can use the split and merge functions of the Ensemble class.
First of all you generate a new ensemble, you compute the energy and force for that ensemble,
then you merge it inside another one.

.. code-block:: python

   # Load the original ensemble (first population with 1000 configurations)
   ens = sscha.Ensemble.Ensemble(dynmat, T, dynmat.GetSupercell())
   ens.load("data_dir", population = 1, N = 1000)

   # Generate a new ensemble with other 1000 configurations
   new_ensemble = sscha.Ensemble.Ensemble(dynmat, T, dynmat.GetSupercell())
   new_ensemble.generate(1000)

   # Compute the energy and forces for the new ensemble
   # For example in this case we assume to have initialized 'calc' as an ASE calculator.
   # But you can also save it with a different population,
   # manually compute energy and forces, and then load again the ensemble.
   new_ensemble.get_energy_forces(calc)

   # Merge the two ensembles
   ens.merge(new_ensemble)

   # Now ens contains the two ensembles. You can save it or directly use it for a SSCHA calculation
   ens.save("data_dir", population = 2)

Indeed, to avoid mistakes, when merging the ensemble you must be carefull that the dynamical matrix and the temperature used to generate both ensembles are the same.


    
How does the error over the gradients scale with the number of configurations?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    The error scales as any stochastic method, with the inverse of the square root of the number of configurations. So to double the accuracy, the number of configurations must be multiplied by 4. 


I cannot remove the pressure anisotropy after relaxing the cell, what is happening?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Variable cell calculation is a tricky algorithm. It could be that your bulk modulus is strongly anisotropic, so the algorithm has difficulties in optimizing well.
In general, the stress tensor is also affected by the stochastic error, so it is impossible to completely remove anisotropy. However, a converged result is one in which the residual anisotropy in the stress tensor is comparable to the stochastic error on the stress tensor.
If you are not able to converge, you can either increase the number of configurations, modify the bulk_modulus parameter (increase it if the stress change too much between two populations, decrease it if it does not changes enough) or fix the overall volume (by using the fix_volume flag in the &relax namespace or the vc_relax method if you are using the python script).

Fixing the volume improves the convergence of the variable cell algorithm (using the fix_volume = True argument of the vc_relax method).


How do I choose the appropriate value of Kong-Liu effective sample size or ratio?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The Kong-Liu (KL) effective sample size is an estimation of how good is the extracted set of configurations to describe the BO landscape around the current values of the dynamical matrix and the centroid position. After the ensemble is generated, the KL sample size matches with the actual number of configurations, however, as the minimization goes, the KL sample size is reduced. The code stops when the KL sample size is below a certain threshold.

The default value for the Kong-Liu threshold ratio is 0.5 (effective sample size = 0.5 the original number of configurations). This is a good and safe value for most situations. However, if you are very far from the minimum and the gradient is big, you can trust it even if it is very noisy. For this reason, you can lower the Kong-Liu ratio to 0.2 or 0.1. However, notice that by construction the KL effective sample size is always bigger than 2.  Therefore, if you use 10 configurations, and you set a threshold ratio below 0.2, you will never reach the threshold, and your minimization will continue forever (going into a very bad regime where you are minimizing something completely random). On the other side, on some very complex systems close to the minimum, it could be safe to increase the KL ratio even at 0.6.




       
Post-processing the output
--------------------------

How do I plot the phonon dispersion after the calculation?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

See :ref:`Plot the phonon dispersion` section.


How do I plot the frequencies of the dynamical matrix during the optimization?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To check if the SSCHA is converging, you should plot the dynamical matrix's frequencies during the minimization.
In particular, you should look if, between different populations, the evolution of each frequency is consistent. If it seems that frequencies are evolving randomly from a population to the next one, you should increase the number of configurations, otherwise, you can keep the number fixed.


The code can print the frequencies at each step.
If you run the code with an input script, you should provide in the &utils tag the filename for the frequencies:

.. code-block:: fortran

   &utils
       save_frequencies = "minim_info"
   &utils
		
You can use the same function from the python script by calling a custom function that saves the frequencies after each optimization step. The Utilities module of the SSCHA offers this function:

.. code-block:: python
		
   IO_freq = sscha.Utilities.IOInfo()
   IO_freq.SetupSaving("minim_info")

   # Initialize the minimizer as minim [...]
   minim.run(custom_function_post = IO_freq.CFP_SaveAll)

Then, while running you can plot all the information about the minimization with:

.. code-block:: console

   $ sscha-plot-data.py minim_info

   And you will see both frequencies, free energy, gradients and everything how it evolves during the
   minimization.

   If you are using a version older than 1.2, the previous command should be replaced with:


.. code-block:: console

   $ plot_frequencies.py minim_info

If you restart the calculation and save it in multiple files, you can concatenate the results with:

.. code-block:: console

   $ sscha-plot-data.py minim_info1 minim_info2 ...

       

Constrains and custom minimization
----------------------------------
 

How do I lock modes from m to n in the minimization?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Constrains to the minimization within the mode space may be added in both the input file (for the stand-alone execution) and in the python script.
In the input script, inside the namespace **&utils**, you should add:

**mu_free_start = 30** and **mu_free_end = 36** : optimize only between mode 30 and 36 (for each q point).

You can also use the keywords **mu_lock_start** and **mu_lock_end** to freeze only a subset of modes.

You can also choose if you want to freeze only the dynamical matrix or also the structure relaxation along with those directions, by picking:

**project_dyn = .true.** and **project_structure = .false.**. In this way, I freeze only the dynamical matrix along with the specified modes, but not the structure.

Modes may be also locked within the python scripting. Look at the LockModes example in the Examples directory.

TODO: Add the same guide for the python code


How do I lock a special atom in the minimization?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    More complex constraints may be activated in the minimization, but their use is limited within the python scripting.
    You can write your constraining function that will be applied to the structure gradient or the dynamical matrix gradient.
    This function should take as input the two gradients (dynamical matrix and structure) and operate directly on them.
    Then it can be passed to the minimization engine as *custom_function_gradient*.

    .. code-block:: python
		    
       LIST_OF_ATOMS_TO_FIX = [0, 2, 3]
       def fix_atoms(gradient_dyn, gradient_struct):
           # Fix the atoms in the list
	   gradient_struct[LIST_OF_ATOMS_TO_FIX, :] = 0

       minim.run( custom_function_gradient = fix_atoms )

    Here, :code:`minim` is the :code:`SSCHA_Minimizer` class. In this case, we only fix the structure gradient. However, the overall gradient will have a translation (acoustic sum rule is violated). Be very careful when doing this kind of constrains, and check if it is really what you want.

    A more detailed and working example that fixes also the degrees of freedom of the dynamical matrix is reported in the FixAtoms example.



    
