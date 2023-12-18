Advanced Features
=================


The python-sscha code can be runned both as a stand-alone application with an input file and as a python library, writing a python script.

We will cover the python scripting as it is more general.

The SSCHA calculation is divided into 3 main steps:
 - The generation of a random ensemble of ionic configurations
 - Calculations of energies and forces on the ensemble
 - The SSCHA free energy minimization

Then this steps are iterated until convergence is achieved.

In this chapter, we cover some advanced features of the SSCHA code, as the manual submission, the use of constrains on modes and atoms or the configurations of cluster with a DFT code different from quantum ESPRESSO.

Manual submission
-----------------

The manual submission allows the user to take full controll over any steps in the simulation. It also means that the code perform just one iteration, and the user must interact with it to provide the forces and energies of the ensemble at each iterations.

It is usefull if you want to have full control on the number of configurations required to converge, or if you simply do not want to configure the automatic submission through a cluster because you have limited resources and are scared that the code could burn too much computer time without you realizing.

Indeed, it is strongly discuraged in variable cell simulations, as the code exploits the results from previous iterations to optimize the cell in a clever way.

The manual submission means that the user manually computes the energies and forces of all the configurations in the ensemble. The code stops after generating the random ensemble, and the user is requested to provide data files that contain the forces and total energies for each configuration.

Thus, the code works in two steps.
In the first step we generate the ensemble. Here is the code

.. code-block:: python

   import sscha, sscha.Ensemble, sscha.SchaMinimizer
   import cellconstructor as CC, cellconstructor.Phonons
   
   # Load the harmonic dynamical matrix
   dyn = CC.Phonons.Phonons('dyn', nqirr = 4)
   
   # If the dynamical matrix contains imaginary frequencies
   # we get rid of them with
   dyn.ForcePositiveDefinite()

   # Now we initialize the ensemble at the target temperature
   TEMPERATURE = 300 # Kelvin
   ensemble = sscha.Ensemble.Ensemble(dyn, TEMPERATURE)

   # We generate the random ensemble with N_CONFIGS configurations
   N_CONFIGS = 1000
   ensemble.generate(N_CONFIGS)

   # We save the ensemble on the disk (inside directory data)
   # We specify an integer 'population' which distinguish several ensembles
   # inside the same directory
   ensemble.save('data', population = 1)

   

To start, we need an initial guess of the dynamical matrix (the dyn file).
The default format is the one of Quantum ESPRESSO, but also phonopy and
ASE formats are supported (refer to the CellConstructor documentation to load these formats). Here we assume that the dynamical matrices are 4 (4 irreducible q points) called 'dyn1', 'dyn2', 'dyn3' and 'dyn4', as the standard quantum espresso format.

The dynamical matrix contain both the information about the atomic structure
and the ionic fluctuations. These can be obtained with a linear response
calculation from DFT.

The previous code generates the ensemble which is stored in the disk.
Inside the data directory you will find a lot of files

The files named 'scf_population1_X.dat' with X going over all the configurations contain the atomic structure in cartesian coordinates. It uses the standard espresso formalism.

You need to compute total energies and forces of each configuration, with your favourite code.
The total energies are written in column inside the file 'total_energies_population1.dat', in Rydberg atomic units and ordered with the index of the configurations.
The forces for each configuration should be inside 'forces_population1_X.dat' in Ry/Borh (Rydberg atomic units).

When you compute energies and forces, you can load them and run the SSCHA minimization:

.. code-block:: python

   import sscha, sscha.Ensemble, sscha.SchaMinimizer
   import cellconstructor as CC, cellconstructor.Phonons
   
   # Load the harmonic dynamical matrix
   dyn = CC.Phonons.Phonons('dyn', nqirr = 4)
   
   # If the dynamical matrix contains imaginary frequencies
   # we get rid of them with
   dyn.ForcePositiveDefinite()

   # Now we initialize the ensemble at the target temperature
   TEMPERATURE = 300 # Kelvin
   ensemble = sscha.Ensemble.Ensemble(dyn, TEMPERATURE)

   # We load the ensemble
   N_CONFIGS = 1000
   ensemble.load('data', population = 1, N = N_CONFIGS)

   # Now we can run the sscha minimization
   minim = sscha.SchaMinimizer.SSCHA_Minimizer(ensemble)
   minim.init()
   minim.run()

   # Print on stdout the final results
   minim.finalize()

   # Save the output dynamical matrix
   minim.dyn.save_qe('final_dyn')
   

And that's it. You run your first manual calculation.


Keep track of free energy, gradients and frequencies during minimization
------------------------------------------------------------------------

It is convenient to store on the file the information during the minimization, as the Free Energy, its gradient values and the frequencies.

To do this, we need to tell the code to save them into a file.

Let us replace the 'minim.run()' line in the previous example with the following code:

.. code-block:: python

   import sscha.Utilities
   IO = sscha.Utilities.IOinfo()
   IO.SetupSaving('minim_data')

   minim.run(custom_function_post = IO.CFP_SaveAll)


If you run it again, the code produces (starting from verison 1.2) two data files: minim_data.dat and minim_data.freqs.
You can plot all the minimization path (frequencies, free energy, gradients) calling the program:

.. code-block:: bash

   $ sscha-plot-data.py minim_data

The sscha-plot-data.py script is automatically installed within the SSCHA code.


Load from the output files
--------------------------

It is possible to load an ensemble directly from a list of output files from a specific program, 
like quantum ESPRESSO.
This is usefull when a calculation ended with an error on the cluster, 
and the ensemble has not been saved on the disk (but you have the output files of the 
configurations that have been already computed succesfully).

In this case, you can restart the minimization using an ensemble with the following code:

.. code-block:: python

    # Load the dynamical matrix
    dyn = CC.Phonons.Phonons("dyn", 4)

    # Load the ensemble
    ens = sscha.Ensemble.Ensemble(dyn, 300)

    # Load the ensemble from the output of the calculator
    # In this case, the pwo files are output of the quantum espresso program.
    # Any output file that ASE is able to read can be used to load the ensemble.
    ens.load_from_calculator_output(directory="data", out_ext=".pwo")


    # Run the minimization
    minim = sscha.SchaMinimizer.SSCHA_Minimizer(ens)
    minim.init()
    minim.set_minimization_step(0.01)
    minim.run()
    minim.finalize()

An example of this procedure is provided in the ``tests/test_load_ensemble_from_calculator_output`` directory.

Alternatively, you can also use the ensemble to restart a full relax procedure.
In this case, you need to provide the key ``restart_from_ens = True`` to the
``relax`` or ``vc_relax`` methods of the ``SSCHA`` class in the ``Relax`` module.

.. code-block:: python

   # Initialize the minimizer minim (see example above)
   relax = sscha.Relax.SSCHA(minim, ase_calculator=calculator,
                           N_configs = 100, max_pop = 20)
   relax.relax(restart_from_ens = True)

   # Or, if you want to perform a variable cell relaxation
   relax.vc_relax(restart_from_ens = True)  



Cluster configuration with a code different from Quantum ESPRESSO
-----------------------------------------------------------------

TODO


Employ a custom function
------------------------

An interesting feature provided by the SSCHA code is the customization of the algorithm. The user has access to all the variables at each iteration of the minimization. 
In this way, the user can print on files additional info or introduce constraints on the structure or on the dynamical matrix.
The interaction between the user and the SSCHA minimization occurs through three functions, that are defined by the user and passed to the **run** method of the **SSCHA_Minimizer** class (in the **SchaMinimizer** module): 
 - custom_function_pre
 - custom_function_gradient
 - custom_function_post

These functions are called by the code before, during, and after each iteration.

The **Utilities** module already provides some basic functions, that can be used for standard purpouses.
For example, the following code employs *custom_function_post* to print on a file the auxiliary dynamical matrix's frequencies at each step.

.. code-block:: python
	
	IO = sscha.Utilities.IOinfo()
	IO.SetupSaving("freqs.dat")
	# .... initialize minim as SSCHA_Minimizer class
	minim.run( custom_function_post = IO.CFP_SaveAll)

In this case *IO.CFP_SaveAll* is the *custom_function_post*. It is a standard python method, that takes one argument (the SSCHA_Minimizer).
*IO.CFP_SaveAll*  prints the frequencies of the current dynamical matrix (stored in minim.dyn) in the filename defined by *IO.SetupSaving("freqs.dat")*.

The following example, we define a *custom_function_post* not provided by the Utilities module. The following code generate a file with the full dynamical matrix for each iteration of the minimization algorithm.

.. code-block:: python
	
	def print_dyn(current_minim):
		# Get the current step id checking the lenght of the __fe__ variable (the free energy)
		step_id = len(current_minim.__fe__)

		# Save the dynamical matrix
		minim.dyn.save_qe("dyn_at_step_{}_".format(step_id))

Here, *print_dyn* is the *custom_function_post*. We must pass it to the *run* method of the *SSCHA_Minimizer* class (minim in the following case).

.. code-block:: python
	
	minim.run(custom_function_post = print_dyn)

In this way, you can interact with the code, getting access to all the variables of the minimization after each step. This could be exploited, for example, to print atomic positions, bond lenght distances or angles during the minimization, or to setup a live self-updating plot of the free energy and its gradient, that automatically refreshes at each step.


Constraints
-----------

Another important case in which you want to interact with the code is to constrain the minimization. 
A standard constraint is the locking of modes, in which you only optimize a subset of phonon branches defined from the beginning. Let us have a look at the code to constrain the modes:

.. code-block:: python

	# [...] Load the initial dynamical matrix as dyn
	ModeLock = sscha.Utilities.ModeProjection(dyn)
	
	# Setup the constrain on phonon branches from 4 to 8 (ascending energy)
	ModeLock.SetupFreeModes(4, 8)
	
	# [...] Define the SSCHA_Minimizer as minim
	minim.run(custom_function_gradient = ModeLock.CFG_ProjectOnModes)

The function *ModeLock.CFG_ProjectOnModes* is the *custom_function_gradient*. It takes two numpy array as input: the gradient of the dynamical matrix and the gradient on the structure.
Since numpy array are pointers to memory allocations, the content of the array can be modified by the function.
The *SSCHA_Minimizer* calls *custom_function_gradient* immediately before emplying the gradient to generate the dyanmical matrix and the structure for the next iteration.
Therefore, *custom_function_gradient* is employed to apply costraints, projecting the gradients in the desidered subspace.

In particular, *CFG_ProjectOnModes* projects the gradient of the dynamical matrix into the subspace defined only by the mode branches selected with *ModeLock.SetupFreeModes*. As done for *custom_function_post*, also here we can define a custom function instead of using the predefined one provided by the *Utilities* module.

The following code limit the projection on the subspace of modes only on the fourth q-point of the dynamical matrix.

.. code-block:: python
		
	iq = 4
	def my_constrain(dyn_gradient, structure_gradient):
		# Let us apply the standard constrain on modes
		ModeLock.CFG_ProjectOnModes(dyn_gradient, structure_gradient)

		# Now we set to zero the gradient of the dynamical matrix if it does not belong to the iq-th q point (ordered as they appear in the dynamical matrix used to initialize the minimization).
		
		nq, nat3, nat3_ = dyn_gradient.shape
		for i in range(nq):
			if i != iq:
				dyn_gradient[i, :, :] = 0

	
	# [...] define minim as the SSCHA_Minimizer 
	minim.run(custom_function_gradient = my_constrain)

The two arguments taken by custom_function_gradient are the gradient of the dynamical matrix of size (nq, 3*nat, 3*nat) and the gradient of the structure of size (nat, 3).
Notice also how, inside *my_constrain*, we call *ModeLock.CFG_ProjectOnModes*. You can concatenate many different custom functions following this approach.

Remember that the gradients are numpy arrays; **you must modify their content accessing their memory using the slices** [x,y,z] as we did.
In fact, if you overwrite the pointer to the memory (defining a new array), the content of the gradient will not be modified outside the function.
In the following code we show an example of correct and wrong.

.. code-block:: python
	
	# This puts the gradient to zero 
	dyn_gradient[:,:,:] = 0  # CORRECT

	# This does not put to zero the gradient
	dyn_gradient = np.zeros( (nq, 3*nat, 3*nat))  # WRONG

In particular, the second expression redefines the name *dyn_gradient* only inside the function, allocating new memory on a different position, and overwriting the name *dyn_gradient* only inside the function to point to this new memory location.  It **does not** write in the memory where *dyn_gradient* is stored: the gradient outside the function is unchanged. 

Indeed, you can also constrain the structure gradient. The ModeLocking class provides a function also to constrain the atomic displacement to follow the lattice vibrations identified by the selected branches at gamma.
This is *ModeLock.CFG_ProjectStructure*. If you want to constrain both the dynamical matrix and the structure, you can simply concatenate them as:

.. code-block:: python
	
	def my_constrain(dyn_grad, structure_grad):
		ModeLock.CFG_ProjectOnModes(dyn_grad, structure_grad)
		ModeLock.CFG_ProjectStructure(dyn_grad, structure_grad)

	# [...]
	minim.run(custom_function_gradient = my_constrain)

Resuming, *custom functions* can be used to inject your personal code inside each SSCHA iteration. Proper use of this function gives you full control over the minimization and allows you to personalize the SSCHA without editing the source code.
