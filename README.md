# python-sscha

Python implementation of the stochastic self-consistent harmonic approximation (SSCHA).
This program, when istalled, contains both the python libraries to perform a 
custom minimization, and manipulate the simulation at your will.

## Why do I need python-sscha?

If you are simulating transport or thermal properties of materials, phase diagrams, or phonon-related properties of materials, than you need python-sscha.
This is a package that enables you to include the effect of both thermal and quantum phonon fluctuations into your *ab initio* simulations.

The method used by this package is the  Stochastic self-consistent Harmonic Approximation (SSCHA). This is a full-quantum method that optimizes the nuclear wave-function (or density matrix at finite temperature) to minimize the free energy.
In this way you can simulatie highly anharmonic systems, like those close to a second order phase transition (as charge density waves and thermoelectric materials). 
Despite the full quantum and thermal nature of the algorithm, the overall computational cost is the comparable to standard classical molecular dynamics. Since the algorithm correctly exploits the symmetries of the crystal, in highly symmetric cases it is also much cheaper. 

python-sscha comes both as a python library that can be runned inside your own workflows and as a stand-alone software, initialized by input scripts with the same syntax as Quantum ESPRESSO.

It can be coupled with any *ab initio* engine for force and energy calculations. It can interact through the Atomic Simulation Environment (ASE), and has an implemented interface for automatic submission of jobs in a remote cluster.

Moreover, it is quite easy to use, a standard input is highly human readable and less than 10 lines!
So, what are you waiting? Download and install python-sscha, and start enjoing the Tutorials!


## Requirements

The requirements of the python-sscha package are:
1. python >= 2.7 
2. numpy
3. matplotlib
3. Lapack
4. Blas
5. gfortran (or any fortran compiler)
6. CellConstructor

These packages are fundamental. In particular CellConstructor is the shoulder on
which python-sscha is builded on. You can find the last development version on
GitHub

The ASE (Atomic Simulation Environment) is another dependency, since CellConstructor
relies on it.


## Installation

To install the package it is recommended the last version of anaconda-python2,
that cames with all the updated numpy and matplotlib, already compiled to work
in parallel. 
Moreover the last version of matplotlib will allow the user to modify the plots 
after they are produced.

It can be simply installed from command line:

```bash
python setup.py install
```

This command must be executed in the same directory as the setup.py script.
Note, if you installed python in a system directory, administration rights may be
requested (add a sudo before the command). 
Please, consider adopting anaconda-python to install the software on clusters where you do not have 
administration rights.

Installation on clusters:
It is suggested to install the package with the anaconda distribution.
For MPI enable parallelism (used by the Lanczos algorithm), you need to specify the MPICC 
environment flag. 
Please, remember that if you use the intel compiler, you need to delete the lapack linking from the
setup.py and include the -mkl (as done for cellconstructor).
Note that you must force to use the same liker compiler as the one used for the compilation.
it may be necessary to specify

MPICC=mpiicc LDSHARED="mpiicc -shared" python setup.py install

A specific setup.py script is provided to install it easily in FOSS clusters.

Remember to use the same compiler that you used to compile the pypar library, 
otherwise when calling MPI_Init() on the initialization, the code will complain and crash.

## Usage

The installer will build and install both the python library and the executable.
The executable can then be runned with:

```bash
sscha -i input_file
```
Where the "input_file" is a valid input for the code. Please look at the SimpleExample
inside the Examples directory to see a template and how to create it.

You are strongly encouraged to follow the Tutorials (in the Tutorial directory).

Alternative the code can be used directly as a python library (expert mode).
This will allow much more customization on the input and the behaviour of the
code, and the easy interface with calculators different from Quantum ESPRESSO.

## ENJOY!
