Customizing the SSCHA
=====================

One of the most useful features of the SSCHA code is the possibility to interact with the minimization procedure. The user has access to all the variables during the minimization, before, during, and after the step. 
In this way, the user can print in output additional information, introduce constraints, or even change the minimization algorithm with a custom one.

The interaction between the user and the SSCHA minimization occurs through three functions: 
custom_function_pre
custom_function_gradient
custom_function_post

These functions are defined by the user, passed to the SSCHA_Minimizer, and called before, during, and after each step. If you practiced a bit with the code, probably you already employed at least one of them.

For example, *custom_function_post* is used to print out the auxiliary dynamical matrix's frequencies at each step. The method that performs this task is standard and provided in the **Utilities** module.

.. code-block:: python
	
	IO = sscha.Utilities.IOinfo()
	IO.SetupSaving("freqs.dat")
	# .... initialize minim as SSCHA_Minimizer class
	minim.run( custom_function_post = IO.CFP_SaveAll)

In this case *IO.CFP_SaveAll* is a function that takes minim as an argument, and prints the frequencies of the current dynamical matrix (stored in minim.dyn) The filename is set by the *IO.SetupSaving("freqs.dat")*.
This is an example of the interaction user-code.
Let us imagine we want to print out the full dynamical matrix at each step. In this case, we must define a custom function:

.. code-block:: python
	
	def print_dyn(current_minim):
		# Get the current step id checking the lenght of the __fe__ variable (the free energy)
		step_id = len(current_minim.__fe__)

		# Save the dynamical matrix
		minim.dyn.save_qe("dyn_at_step_{}_".format(step_id))

Then, we pass this function to the SSCHA_Minimizer class (minim in the following case) during the run:

.. code-block:: python
	
	minim.run(custom_function_post = print_dyn)

In this way, you can interact with the code, getting access to all the variables of the minimization after each step.

Another important case in which you want to interact with the code is to constrain the minimization. 
A standard constraint is the locking of modes, in which you only optimize a subset of phonon branches defined from the beginning. Let us have a look at the code to constrain the modes:

.. code-block:: python

	# [...] Load the initial dynamical matrix as dyn
	ModeLock = sscha.Utilities.ModeProjection(dyn)
	
	# Setup the constrain on phonon branches from 4 to 8 (ascending energy)
	ModeLock.SetupFreeModes(4, 8)
	
	# [...] Define the SSCHA_Minimizer as minim
	minim.run(custom_function_gradient = ModeLock.CFG_ProjectOnModes)

The function *ModeLock.CFG_ProjectOnModes* is the one that interacts with the code. It takes as input the gradient of the dynamical matrix and the gradient on the structure. They can be modified by the function, in this way, before applying the step, you can set to zero some components of the gradient, or project the gradients in a subspace. This allows you to constrain the minimization. In particular, *CFG_ProjectOnModes* projects the gradient of the dynamical matrix into the subspace defined only by the mode branches selected with *ModeLock.SetupFreeModes*. 
They are done for each q point at which the dynamical matrix is defined. Let us imagine that you want only to minimize some phonon branches at a particular q point. We can modify this function as:

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

As you see, here we are calling inside *my_constrain* the *ModeLock.CFG_ProjectOnModes*. You can concatenate many different custom functions following this approach. In particular, the two arguments taken by custom_function_gradient are the gradient of the dynamical matrix of size (nq, 3*nat, 3*nat) and the gradient of the structure of size (nat, 3). Remember that they are array, so you can access their memory using the slices [x,y,z] as we did. However, if you redefine them inside the function, the gradient values will not be modified:

.. code-block:: python
	
	# This puts the gradient to zero correctly
	dyn_gradient[:,:,:] = 0

	# This does not put to zero the gradient
	dyn_gradient = np.zeros( (nq, 3*nat, 3*nat))

In particular, the second expression will redefine the name dyn_gradient inside the custom_function, allocating new memory. Thus, when the code exit from the custom_function_gradient and returns to the minimization, it will not have overwritten the memory part of the gradient, that will result unaffected. Instead, the first method (the one with slices) access directly the memory of dyn_gradient overwriting it, and the result will be visible also outside the function. This is how Python works, but take care, as if you are not aware of it, it can take you some headache to figure out why the code is not imposing your constrain correctly.

Indeed, you can also constrain the structure gradient. the ModeLocking class provides a function also to constrain the atomic displacement to follow the lattice vibrations identified by the selected branches at gamma.
This is *ModeLock.CFG_ProjectStructure*. If you want to constrain both the dynamical matrix and the structure, you can simply concatenate them as:

.. code-block:: python
	
	def my_constrain(dyn_grad, structure_grad):
		ModeLock.CFG_ProjectOnModes(dyn_grad, structure_grad)
		ModeLock.CFG_ProjectStructure(dyn_grad, structure_grad)

	# [...]
	minim.run(custom_function_gradient = my_constrain)

Resuming, *custom functions* can be used to inject your personal code inside each SSCHA iteration. Proper use of this function gives you full control over the minimization and allows you to personalize the SSCHA without editing the source code.
