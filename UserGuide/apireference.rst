THE API
=======

This chapter contains the documentation for the main methods of the python-sscha code.
It can be used both by advanced users, that wants to exploit python-sscha as a library,
or developers, willing to add new features to the code (or adapt existing ones for their purposes).

The API is divided into Modules.


The Ensemble Module
-------------------

This module deals with the ensembles of configurations.
It is used to generate random configurations from the dynamical matrix, to compute observables on the ensemble used in the SSCHA optimization.
These include the average force on atoms, the gradient of the SSCHA minimization, the quantum-thermal stress tensor, as well as properties of
the ensemble, like reweighting.

.. autoclass:: sscha.Ensemble.Ensemble
   :members:

The SchaMinimizer Module
------------------------

This module is the main SSCHA minimizer. It allows us to set up a single (one population) minimization.
In this module, the minimization algorithm is introduced, as well as stopping conditions and all the parameters
usually located in the &inputscha name list are read.
      
.. autoclass:: sscha.SchaMinimizer.SSCHA_Minimizer
   :members: 
      

The Relax Module
----------------

This module deals with relaxations that are iterated over more populations. It includes the variable cell optimization algorithm.
Here the parameters read in the &relax name list are read and setup.

.. autoclass:: sscha.Relax.SSCHA
   :members:


      
The Utilities Module
----------------

This module both provides the constrains and the IOinput

IOInfo class
^^^^^^^^^^^^

Use this class to make the python-sscha print information during the minimization

.. autoclass:: sscha.Utilities.IOInfo
   :members:


Constraints
^^^^^^^^^^^

The constrains are a member of the Utilities module.
To implement constrains on phonon modes, use the ModeProjection class


.. autoclass:: sscha.Utilities.ModeProjection
   :members:


      

The Cluster Module
------------------

The Cluster module provides the interface between python-sscha and remote servers to which you submit the energy and forces calculations.
The input in &cluster namespace is interpreted in this module

.. autoclass:: sscha.Cluster.Cluster
   :members:

      
