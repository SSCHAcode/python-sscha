Customizing the SSCHA
=====================

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
