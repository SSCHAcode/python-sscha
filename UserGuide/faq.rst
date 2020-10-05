Frequently Asked Questions (FAQs)
=================================

.. role:: raw-html(raw)
    :format: html

How do I start a calculation if the Dynamical matrices have imaginary frequencies?
    :raw-html:`<br />`
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
		    
       python script.py

  

What are the reasonable values for the steps (lambda_a, lambda_w, min_step_dyn, and min_step_struc)?
    :raw-html:`<br />`
    The code minimizes using a Newton method: preconditioned gradient descend. Thanks to an analytical evaluation of the Hessian matrix, the step is rescaled so that the theoretical best step is close to 1.
    In other words: **one is theoretically the best (and the default) choice for the steps**. However, the SSCHA is a stochastic algorithm, therefore, if the ensemble is too small, or the gradient is very big, this step could bring you outside the region in which the ensemble is describing well the physics very soon.
    Since SSCHA can exploit the reweighting, and the most computationally expensive part of the algorithm is the computation of forces and energies, it is often much better using a small step (smaller than the optimal one). **Good values of the steps are usually around 0.01 and 0.1**. Rule of thumbs: the minimization should not end because it went outside the stochastic regime before that at least 10 steps have been made. This will depend on the system, the number of configurations, and how far from the correct solution you are.

    **lambda_w** is the step in the atomic positions (stand-alone program input).
    
    **lambda_a** is the step in the dynamical matrix (stand-alone program input).

    If you are using the python script, the equivalent variables are the attributes of the sscha.SchaMinimizer.SSCHA_Minimizer class.
  
    **min_step_struc** is the step in the atomic positions (stand-alone program input).
    
    **min_step_dyn** is the step in the dynamical matrix (stand-alone program input).  
  

In a variable cell optimization, what is a reasonable value for the bulk modulus?
    :raw-html:`<br />`
    The bulk modulus is just an indicative parameter used to guess the optimal step of the lattice parameters to converge as quickly as possible.
    It is expressed in GPa. You can find online the bulk modulus for many materials. Find a material similar to the one you are studying and look if there is in literature a bulk modulus.

    Usual values are between 10 GPa and 100 GPa for a system at ambient conditions. Diamond has a bulk modulus of about 500 GPa. High-pressure hydrates have a bulk modulus of around 300 GPa as well.

    If you have no idea on the bulk modulus, you can easily compute them by doing two static *ab initio* calculations at very close volumes (by varying the cell size), and then computing the differences between the pressure:

    .. math::

       B = - \Omega \frac{dP}{d\Omega}

    where :math:`\Omega` is the unit-cell volume and :math:`P` is the pressure (in GPa).


The code stops saying it has found imaginary frequencies, how do I fix it?
    :raw-html:`<br />`
    This means that your step is too large. You can reduce the step of the minimization. An alternative (often more efficient) is to switch to the root representation.
    In this way, the square root of the dynamical matrix is minimized, and the total dynamical matrix is positive definite in the whole minimization by construction.

    In the input name-list, you activate this minimization with the following keywords inside the &inputscha name-list

    .. code-block:: fortran
		    
       preconditioning = .false.
       root_representation = "root4"

    Or, in the python script, you set up the attributes of the sscha.SchaMinimizer.SSCHA_Minimizer class

    .. code-block:: python
		    
       minim.preconditioning = False
       minim.root_representation = "root4"

    The optimal step size for the root_representation may be different than the other one.
     

    
Why the gradient sometimes increases during a minimization?
    :raw-html:`<br />`
    Noting in principle assures that a gradient should always go down. It is possible at the beginning of the calculation when we are far from the solution that one of the gradients increases.
    However, when we get closer to the solution, indeed the gradient must decrease.
    If this does not happen it could be due to the ensemble that has fewer configurations than necessary. In this case, the good choice is to increase the number of effective sample size (the Kong-Liu ratio), to stop the minimization when the gradient starts increasing, or to increase the number of configurations in the ensemble.

    In any case, what must decrease is free energy. If you see that the gradient is increasing but the free energy decreases, then the minimization is correct. However, if both the gradient and free energy are increasing, something is wrong. This could be due to a step size too big, then try to reduce the value of **lambda_a** and **lambda_w** (in the input file) or **min_step_dyn** and **min_step_struc** (in the python script). It could also be due to a wasted ensemble, in this case, check the value of the Kong-Liu effective sample size, if it is below or around 0.5, then try to increase the threshold at which stop the calculation, **kong_liu_ratio** (in the python script) or **N_random_eff** (in the input file), or increase the number of configurations for the next population.


The gradients on my simulations are increasing a lot, why is this happening?
    :raw-html:`<br />`
    See the previous question. 


How do I check if my calculations are well converged?
    :raw-html:`<br />`
    In general, if the gradient goes to zero and the Kong Liu ratio is above 0.5 probably your calculation converged very well.
    There are some cases (especially in systems with many atoms) in which it is difficult to have an ensemble sufficiently big to reach this condition.
    In these cases, you can look at the history of the frequencies in the last populations.

    If the code is provided with a &utils namespace, on which the code

    .. code-block:: console
		    
       &utils
          save_freq_filename = "frequencies_popX.dat"
       &end

    You can after the minimization use the plotting program to see the frequencies as they evolve during the minimizations:

    .. code-block:: console
		    
       plot_frequencies_new.pyx frequencies_pop*.dat

       
    This will plot all the files *frequencies_popX.dat* in the directory. You can see all the history of the frequency minimization.
    If between different populations (that you will distinguish by a kink in the frequency evolutions) the frequencies will fluctuate due to the stochastic nature of the algorithm, with no general drift, then the algorithm reached its maximum accuracy with the given number of configurations.
    You may either stop the minimization or increase the ensemble to improve the accuracy.


What is the final error on the structure or the dynamical matrix of a SCHA minimization?
    :raw-html:`<br />`
    This is a difficult question. The best way to estimate the error is to generate a new ensemble with the same number of configurations at the end of the minimization and check how the final optimized solution changes with this new ensemble. This is also a good way to test if the solution is converged to the correct solution. The magnitude of the changes in the dynamical matrix's frequencies and structure is an accurate estimation of the stochastic error.

    You can always split the ensemble in two and run two minimizations with the two half of the ensemble to get a hint on the error on the structure or the dynamical matrix.
    To split the ensemble, refer to the *FAQ* about the error on the hessian matrix.


How does the error over the gradients scale with the number of configurations?
    :raw-html:`<br />`
    The error scales as any stochastic method, with the inverse of the square root of the number of configurations. So to double the accuracy, the number of configurations must be multiplied by 4. 


When I relax the cell, is it necessary for the gradients to reach zero before making a step with the new cell?
    :raw-html:`<br />`
    In general, it is good to have a reasonable dynamical matrix before starting with a variable cell relaxation. The best strategy is to perform a fixed cell relaxation with few configurations until you are close to the final solution (the gradients are comparable with their errors). Then you can start a variable cell relaxation and submit new populations in the suggested new cell even if the previous one was not perfectly converged.

I cannot remove the pressure anisotropy after relaxing the cell, what is happening?
    :raw-html:`<br />`
    Variable cell calculation is a tricky algorithm. It could be that your bulk modulus is strongly anisotropic, so the algorithm has difficulties in optimizing well.
    In general, the stress tensor is also affected by the stochastic error, so it is impossible to completely remove anisotropy. However, a converged result is one in which the residual anisotropy in the stress tensor is comparable to the stochastic error on the stress tensor.
    If you are not able to converge, you can either increase the number of configurations, modify the bulk_modulus parameter (increase it if the stress change too much between two populations, decrease it if it does not changes enough) or fix the overall volume (by using the fix_volume flag in the &relax namespace or the vc_relax method if you are using the python script).
    Fixing the volume can improve the convergence of the variable cell algorithm by a lot.

    
How may I run a calculation neglecting symmetries?
    :raw-html:`<br />`
    You can tell the code to neglect symmetries with the :code:`neglect_symmetries = .true.` flag.
    In the python script, this is done setting the attribute *neglect_symmetries* of sscha.SchaMinimizer.SSCHA_Minimizer to False.


In which units are the lattice vectors, the atomic positions, and the mass of the atoms in the dynamical matrix file?
    :raw-html:`<br />`
    The dynamical matrix follows the quantum espresso units. They are Rydberg atomic units (unit of mass is 1/2  the electron mass, energy is Ry, positions are in Bohr. However, espresso may have an ibrav not equal to zero (the third number in the header of the dynamical matrix). In this case, please, refer to the espresso ibrav guide in the `PW.x input description <https://www.quantum-espresso.org/Doc/INPUT_PW.html#idm199>`
  

What is the difference between different kinds of minimization (preconditioning and root_representation)?
    :raw-html:`<br />`
    The target of a SSCHA minimization is to find the ionic density matrix :math:`\rho(\Phi, \vec {\mathcal R})` that minimizes the total free energy. It may happen, if we are using a too big step for the dynamical matrix :math:`\Phi` that it becomes not positive definite.
    This may be due to the stochastic noise during the minimization.
    For avoid this to happen, you may set **root_representation** to either **sqrt** or **root4** (inside &inputscha namespace or the SSCHA_Minimizer object)
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

    
 

How do I lock modes from m to n in the minimization?
    :raw-html:`<br />`
    Constrains to the minimization within the mode space may be added in both the input file (for the stand-alone execution) and in the python script.
    In the input script, inside the namespace **&utils**, you should add:

    **mu_free_start = 30** and **mu_free_end = 36** : optimize only between mode 30 and 36 (for each q point).

    You can also use the keywords **mu_lock_start** and **mu_lock_end** to freeze only a subset of modes.

    You can also choose if you want to freeze only the dynamical matrix or also the structure relaxation along with those directions, by picking:

    **project_dyn = .true.** and **project_structure = .false.**. In this way, I freeze only the dynamical matrix along with the specified modes, but not the structure.

    Modes may be also locked within the python scripting. Look at the LockModes example in the Examples directory.


How do I lock a special atom in the minimization?
    :raw-html:`<br />`
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



How do I understand if I have to generate a new population or the minimization converged?
    :raw-html:`<br />`
    In general, if the code stops because the gradient is much below the error (less then 1\%), then it is converged (with a Kong-Liu threshold ratio of at least 0.5). If the code ends the minimization because it went outside the stochastic criteria, a new population is required.
    There are cases in which you use to few configurations to reach a small gradient before wasting the ensemble. If this is the case, print the frequencies during the minimizations (using the &utils card with :code:`save_freq_filename` attribute). You may compare subsequent minimizations, if the frequencies are randomly moving between different minimization (and you cannot identify a trend in any of them), then you reach the limit of accuracy of the ensemble. Frequencies are a much better parameter to control for convergence than free energy, as the free energy close to the minimum is quadratic.



How do I choose the appropriate value of Kong-Liu effective sample size or ratio?
    :raw-html:`<br />`
    The Kong-Liu (KL) effective sample size is an estimation of how good is the extracted set of configurations to describe the BO landscape around the current values of the dynamical matrix and the centroid position. After the ensemble is generated, the KL sample size matches with the actual number of configurations, however, as the minimization goes, the KL sample size is reduced. The code stops when the KL sample size is below a certain threshold.
    
    The default value for the Kong-Liu threshold ratio is 0.5 (effective sample size = 0.5 the original number of configurations). This is a good and safe value for most situations. However, if you are very far from the minimum and the gradient is big, you can trust it even if it is very noisy. For this reason, you can lower the Kong-Liu ratio to 0.2 or 0.1. However, notice that by construction the KL effective sample size is always bigger than 2.  Therefore, if you use 10 configurations, and you set a threshold ratio below 0.2, you will never reach the threshold, and your minimization will continue forever (going into a very bad regime where you are minimizing something completely random). On the other side, on some very complex systems close to the minimum, it could be safe to increase the KL ratio even at 0.6.


How do I understand if the free energy hessian calculation is converged?
    :raw-html:`<br />`
    The free energy hessian requires much more configurations than the SCHA minimization. First of all, to run the free energy Hessian, the SSCHA minimization must end with a gradient that can be decreased indefinitely without decreasing the KL below 0.7 /0.8.
    Then you can estimate the error by repeating the hessian calculation with half of the ensemble and check how the frequencies of the hessian changes. This is also a good check for the final error on the frequencies.
    
    You can split your ensemble in two by using the split function.

    .. code-block:: python

       import sscha, sscha.Ensemble

       # Load the dynamical matrix as dyn
       # [...]
       
       # ens is the Ensemble() class correctly initialized.
       # We can for example load it
       # Assuming it is stored in 'data_dir' with population 1 and 1000 configurations
       # We assume to have loaded the original dynamical matrix dyn and to know the generating temperature T
       ens = sscha.Ensemble.Ensemble(dyn, T, dyn.GetSupercell())
       ens.load("data_dir", population = 1, N = 1000)

       # We create a mask that selects which configurations to take
       first_half_mask = np.zeros(ens.N, dtype = bool)
       first_half_mask[: ens.N // 2] = True

       # We create also the mask for the second half
       # by taking the not operation on the first_half_mask
       second_half_mask = ~first_half_mask

       # Now we split the ensemble
       ens_first_half = ens.split(first_half_mask)
       ens_second_half = ens.split(second_half_mask)

       # We can save the two half ensembles as population 2 and 3.
       ens_first_half.save("data_dir", population = 2)
       ens_second_half.save("data_dir", population = 3)

    This simple script will generate two ensembles inside :code:`data_dir` directory with population 2 and 3, each one containing the first and the second half of the ensemble with population 1 respectively. You can perform then your calculation of the free energy Hessian with both the ensemble to estimate the error on the frequencies and the polarization vectors.
       


How can I add more configurations to an existing ensemble?
    :raw-html:`<br />`
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

    Indeed, to avoid mistakes, when merging the ensemble you must be carefull that the dynamical matrix and the temperature
    used to generate both ensembles are the same.
    

How do I fix the random number generator seed to make a calculation reproducible?
    :raw-html:`<br />`
    As for version 1.0, this can be achieved only by using the python script.
    Since python uses NumPy for random numbers generation, you can, at the beginning of the script that generates the ensemble, use the following:

    .. code-block:: python
		    
       import numpy as np

       X = 0
       np.random.seed(seed = X)

    where :code:`X` is the integer used as a seed. By default, if not specified, it is initialized with None that it is equivalent to initializing with the current local time.
