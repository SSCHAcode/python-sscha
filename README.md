Introduction
============

Website for more information [sscha.eu](http://sscha.eu/).

The full documentation of the package is in the [python-sscha.pdf](https://github.com/SSCHAcode/python-sscha/blob/master/python-sscha.pdf) file, in the root directory of this repository

What is python-sscha?
---------------------

python-sscha is both a python library and a stand-alone program to simulate quantum and thermal fluctuations in solid systems.




Why do I need python-sscha?
---------------------------


If you are simulating transport or thermal properties of materials, phase diagrams, or phonon-related properties of materials, you need python-sscha.
It is a package that enables you to include the effect of both thermal and quantum phonon fluctuations into your *ab initio* simulations.

The method used by this package is the  Stochastic self-consistent Harmonic Approximation (SSCHA). The SSCHA is a full-quantum method that optimizes the nuclear wave-function (or density matrix at finite temperature) to minimize the free energy.
In this way, you can simulate highly anharmonic systems, like those close to a second-order phase transition (as charge density waves and thermoelectric materials). 
Despite the full quantum and thermal nature of the algorithm, the overall computational cost is comparable to standard classical molecular dynamics. Since the algorithm correctly exploits the symmetries of the crystal, it is also much cheaper. 

python-sscha comes both as a python library that can be run inside your workflows and as stand-alone software, initialized by input scripts with the same syntax as Quantum ESPRESSO.

You can couple it with any *ab initio* engine for force and energy calculations. It can interact through the Atomic Simulation Environment (ASE) and has an implemented interface for automatic submission of jobs in a remote cluster.

Moreover, it is easy to use, with short input files highly human-readable.
What are you waiting for? Download and install python-sscha, and start enjoying the Tutorials!


How to install
==============

The SSCHA code is a collection of 2 python packages: CellConstructor and python-sscha.
In this guide, we refer to the installation of python-sscha.


Requirements
------------

To install python-sscha you need:
1. python (either 2.7 or 3.*)
2. numpy
3. scipy
4. matplotlib
5. Lapack
6. Blas
7. gfortran (or any fortran compiler)
8. CellConstructor

For python, we strongly recommend using the anaconda distribution, that already comes with numerical packages correctly compiled to exploit multithreading.

The numpy, scipy and matplotlib are python packages. These are usually provided with a new installation
of python distributions like anaconda. Lapack and Blas are needed for compiling the FORTRAN code (together with a FORTRAN compiler like gfortran).
In many Linux distributions like ubuntu they can be installed as 

.. code-block:: console

   sudo apt-get install libblas-dev liblapack-dev liblapacke-dev gfortran



Note that this specific command may change in time. 


Together with these mandatory requirements (otherwise, the code will not compile correctly or raise an exception at the startup), we
strongly recommend installing also the following libraries:
1. Atomic Simulation Environment (ASE)
2. SPGLIB

If these packages are available, they will enable the automatic cluster/local calculation of forces (ASE) and the symmetry recognition in the supercell (SPGLIB).


To install all the python dependencies (and recommended) automatically, you may just run:

.. code-block:: console
   
   pip install -r requirements.txt




Installation from pip
---------------------

The easiest way to install python-sscha (and CellConstructor) is through the python package manager:

.. code-block:: console
   
   pip install python-sscha 



Eventually, you can append the --user option to install the package only for the user (without requiring administrator powers).
Pip will check for requirements automatically and install them. This method only works if pip is already installed with python.



Installation from source
------------------------

Once all the dependences of the codes are satisfied, you can unzip the source code downloaded from the website.
Then run, inside the directory that contains the setup.py script, the following command:

.. code-block:: console

   python setup.py install


As for the pip installation, you may append the --user option to install the package only for the user (without requiring administrator powers).


Install with Intel FORTRAN compiler
-----------------------------------

The setup.py script works automatically with the GNU FORTRAN compiler. However, due to some differences in linking lapack,
to use the intel compiler you need to edit a bit the setup.py script:

In this case, you need to delete the lapack linking from the
setup.py and include the -mkl as linker option.
Note that you must force to use the same liker compiler as the one used for the compilation. 

Install with a specific compiler path
-------------------------------------

This can be achieved by specifying the environment variables on which setup.py relies:

1. CC (C compiler)
2. FC (Fortran compiler)
3. LDSHARED (linking)

If we want to use a custom compiler in /path/to/fcompiler we may run the setup as:

.. code-block:: console

   FC=/path/to/fcompiler LDSHARED=/path/to/fcompiler python setup.py install



A specific setup.py script is provided to install it easily in FOSS clusters.


Quick start
===========


To quickly start using the code, we recommend using the jupyter notebooks with examples we provide in the Tutorials directory of the source code.

Tutorials are organized as follows:

1. Setup the first calculation: PbTe tutorial. Here you learn how to set up a SSCHA calculation starting just with the structure (we provide a .cif file of the PbTe at high temperature). The tutorial will guide you step by step. You will learn how to: prepare the starting data needed for the SSCHA calculation, generate a random ensemble, save the ensemble and prepare input files for your favorite ab-initio code, read back the energies and the forces inside SSCHA, run a SSCHA minimization. You will also learn how to use ASE and the Cluster module to automatize the calculation of the ensemble and submit it to a HPC system.
2. Automatic relaxation with a force field: SnTe_ToyModel. Here, we show how to use a force-field for a SSCHA calculation, running everything on your computer. We also will explain how to calculate the free energy hessian for second-order phase transitions, and study a phase transition as a function of temperature.
3. Variable cell relaxation: LaH10 tutorial. Here you learn how to perform an automatic calculation with a variable cell. You will exploit the quantum effect to search the high-temperature superconductive phase (Fm-3m) of LaH10 below 200 GPa, starting from a distorted structure. 
4. Hessian matrix calculation for second-order phase transitions: H3S tutorial. Here you reproduce the calculation of the Hessian of the free energy to assert the stability of the H3S phase.
5. Spectral properties: Spectral_Properties. In this tutorial, we explain how to use the post-processing utilities of the SSCHA to calculate the phonon spectral function, and computing phonon lifetimes, and plotting interacting phonon dispersion. We provide an ensemble for PbTe already computed ab-initio.


The jupyter notebooks are interactive, to quickly start with your simulation, pick the tutorial that resembles the kind of calculation you want to run, and simply edit it directly in the notebook. 

Stop
