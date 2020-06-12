Frequently Asked Questions (FAQs)
=================================

* Q: How do I start a calculation if the Dynamical matrices have imaginary frequencies?

  A: Good starting point for a sscha minimization are the dynamical matrix obtained from a harmonic calculation. However, they can have imaginary frequencies. This may be related to both instabilities (the structure is a saddle-point of the Born-Oppenheimer energy landscape) or to a not well converged choice of the parameters for computing the harmonic frequencies.

  In both cases, it is very easy to get a new dynamical matrix that is positive definite and can be used as starting point. An example is made in Turorial on H3S.
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
  ```bash
  >>> python script.py
  ```
  

* Q: What are the reasonable values for the steps (lambda_a, lambda_w, min_step_dyn and min_step_struc)?

  A: The code minimizes using a Newton method: preconditioned gradient descend. Thanks to an analytical evaluation of the hessian matrix, the step is rescaled so that the theoretical best step is close to 1.
  In other words: **one is theoretically the  best (and the default) choice for the steps**. However, the SSCHA is a stochastic algorithm, therefore, if the ensemble is too small, or the gradient is very big, this step could bring you outside the region in which the ensemble is describing well the physics very soon.
  Since SSCHA can exploit the reweighting, and the most computational expensive part of the algorithm is the computation of forces and energies, it is often much better using a small step (smaller than the optimal one). **Good values of the steps are usually around 0.01 and 0.1**. Rule of thumbs: the minimization should not end because it whent outside the stochastic regime before that at least 10 steps have been made. This will depend on the system, the number of configurations and how far from the correct solution you are.

  **lambda_w** is the step in the atomic positions (stand-alone program input)
  **lambda_a** is the step in the dynamical matrix (stand-alone program input)

  If you are using the python script, the equivalent variables are the attributes of the sscha.SchaMinimizer.SSCHA_Minimizer class.
  
  **min_step_struc** is the step in the atomic positions (stand-alone program input) 
  **min_step_dyn** is the step in the dynamical matrix (stand-alone program input)  
  
* Q: In a variable cell optimization, what is a reasonable value for the bulk modulus?

  The bulk modulus is just an indicative parameter used to guess the optimal step of the lattice parameters in order to converge as quickly as possible.
  It is expressed in GPa. You can find online the bulk modulus for many materials. Find a material similar to the one you are studying and look if there is in letterature a bulk modulus.

  Usual values are between 10 GPa and 100 GPa for system at ambient conditions. Diamond has a bulk modulus about 500 GPa. High pressure hydrates have a bulk modulus around 500 GPa as well.

  If you have no idea on the bulk modulus, you can easily compute them by doing two static *ab initio* calculations at very close volumes (by varying the cell size), and then computing the differences between the pressure:

  .. math::

     B = - \Omega \frac{dP}{d\Omega}

  where :math:`\Omega` is the unit-cell volume and :math:`P` is the pressure (in GPa).

* Q: The code stops saying it has found imaginary frequencies, how do I fix it?

  This means that you step is too large. You can reduce the step of the minimization. An alternative (often more efficient) is to switch to the root representation.
  In this way the square root of the dynamical matrix is minimized, and the total dynamical matrix is positive definite in the whole minimization by construction.

  In the namelist input you activate this minimization with the following keywords inside the &inputscha namelist

  .. code-block:: fortran
     preconditioning = .false.
     root_representation = "root4"

  Or, in the python script, you setup the attributes of the sscha.SchaMinimizer.SSCHA_Minimizer class

  .. code-block:: python
     minim.preconditioning = False
     minim.root_representation = "root4"

  It is possible that the optimal step size for the root_representation is different than the other one.
     
* Q: While the gradient on the atomic positions goes down, the gradient on the dynamical matrix is going up (or vice versa), what is happening?

  Noting in principle assures that a gradient should always go down. It is possible at the beginning of the calculation when we are far from the solution that one of the gradients increases.
  However, when we get closer to the solution, indeed the gradient must decrease.
  If this does not happen it could be due to the ensemble that has too less configurations. In this case, the good choice is to increase the number of effective sample size (the kong-liu ratio), in order to stop the minimization when the gradient start increasing, or to increase the number of configurations in the ensemble.

* Q: The gradients on my simulations are increasing a lot, why is this happening?

 See the previous question.

* Q: How do I check if my calculations are well converged?

  In general, if the gradient goes to zero and the Kong Liu ratio is above 0.5 probably your calculation converged very well.
  There are some cases (especially in systems with many atoms) in which it is difficult to have an ensemble sufficiently big to reach this condition.
  In these cases, you can look at the history of the frequencies in the last populations.

  If the code is provided with a &utils namespace, on which the code
  .. code-block:: fortran
     &utils
        save_freq_filename = "frequencies_popX.dat"
     &end

  You can after the minimization use the plotting program to see the frequencies as they evolve during the minimizations:
  ```bash
  >>> plot_frequencies_new.pyx frequencies_pop*.dat
  ```
  This will plot all the files *frequencies_popX.dat* in the directory. You can see all the history of the frequency minimization.
  If between different populations (that you will distinguish by kink in the frequency evolutions) the frequencies will fluctuate due to the stochastic nature of the algorithm, with no general drift, then the algorithm reached its maximum accuracy with the given number of configurations.
  You may either stop the minimization, or increase the ensemble to improve the accuracy.

* Q: What is the final error on the structure or the dynamical matrix of a SCHA minimization?
  That is a difficult question... XXX TO BE ANSWERED  XXX

* Q: How does the error over the gradients scale with the number of configurations?

  The error scales as any stochastic method, with the inverse of the square root of the number of configurations. 

* Q: When I relax the cell, is it necessary for the gradients to reach zero at the end of the run before making a step with the new cell?

  In general it is good to have a reasonable dynamical matrix before starting with a variable cell relaxation. The best strategy is to perform a fixed cell relaxation with few configurations until you are close to the final solution (the gradients are comparable with their errors). Then you can start a variable cell relaxation and submit new populations in the suggested new cell even if the previous one was not perfectly converged.

* Q: I cannot remove the pressure anisotropy after relaxing the cell, what is happening?

  Variable cell calculation is a tricky algorithm. It could be that your bulk modulus is stronlgy anisotropic, so the algorithm has difficulties in optimizing well.
  In general the stress tensor is also affected by stochastic error, so it is impossible to completely remove anisotropy. However, a converged result is one in which the residual anisotropy in the stress tensor is comparable to the stochastic error on the stress tensor.
  If you are not able to converge, you can either increase the number of configurations, modify the bulk_modulus parameter (increase it if the stress change too much between two populations, decrease it if it does not changes enough) or fix the overall volume (by using the fix_volume flag in the &relax namespace or in the vc_relax method if you are using the python script).
  Fixing the volume can improve the convergence of the variable cell algorithm by a lot.

* Q: How may I run a calculation neglecting symmetries?

  You can tell the code to neglect symmetries with the neglect_symmetries = .true. flag.
  In the python script this is done setting the attribute *neglect_symmetries* of sscha.SchaMinimizer.SSCHA_Minimizer to False.

* Q: I am looking at the dynamical matrices obtained after running the SSCHA, in which units are the vectors of the new obtained cell, the atomic positions and the mass of the atoms?

  The dynamical matrix follows the quantum espresso units. They are Rydberg atomic units (unit of mass is 1/2  the electron mass, energy is Ry, positions are in Bohr. However, espresso may have an ibrav not equal to zero (the third number in the header of the dynamical matrix). In this case, please, refer to the espresso ibrav guide in the `PW.x input description <https://www.quantum-espresso.org/Doc/INPUT_PW.html#idm199>`
  

* Q: What is the difference between the different kind of minimization (preconditioning and root_representation)?
  


* Q: How do I lock modes from m to n in the minimization ?


* Q: How do I lock a special atom in the minimization ?


* Q: How do I understand if I have to generate a new population or the minimization converged?


* Q: How do I choose the appropriate value of Kong-Liu effective sample size or ratio?

* Q: How do I understand if the free energy hessian calculation is converged ?

* Q: How can I add more configurations to an existing ensembe?

* Q: How do I fix the random number generator seed to make a calculation reproducible?

