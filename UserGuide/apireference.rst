THE API
=======

This chapter contains the documentations for the main methods of the python-sscha code.
It can be used both from advanced users, that wants to exploit python-sscha as a library,
or developers, willing to add new features to the code (or adapt existing ones for their own purpouses).

The API is divided in Modules.


The Ensemble Module
-------------------

This module deals with ensembles of configurations.
It is used to generate random configurations from the dynamical matrix, to compute observables on the ensemble used in the SSCHA optimization.
These includes the average force on atoms, the gradient of the SSCHA minimization, the quantum-thermal stress tensor, as well as properties of
the ensemble, like reweighting.

.. autoclass:: sscha.Ensemble
   :members:

The SchaMinimizer Module
------------------------

This module is the main SSCHA minimizer. It allows to setup a single (one population) minimization.
In this module the minimization algorithm is introduced, as well as stopping conditions and all the parameters
usually located in the &inputscha namelist are read.
      
.. automodule:: sscha.SchaMinimizer
   :members: 
      

The Relax Module
----------------

This module deals with relaxations that are iterated over more populations. It includes the variable cell optimization algorithm.
Here the parameters readed in the &relax namelist are read and setup.

.. automodule:: sscha.Relax
   :members:


The Cluster Module
------------------

The Cluster module provide the interface between python-sscha and remote servers to witch you submit the energy and forces calculations.
The input in &cluster namespace is interpreted in this module

.. automodule:: sscha.Cluster
   :members:

      
