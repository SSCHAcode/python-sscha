# python-sscha

Python implementation of the stochastic self-consistent harmonic approximation (SSCHA).
This program, when istalled, contains both the python libraries to perform a 
custom minimization, and manipulate the simulation at your will.

## Requirements

The requirements of the python-sscha package are:
1. python >= 2.7 and < 3
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

It is recommen

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

## Usage

The installer will build and install both the python library and the executable.
The executable can then be runned with:

```bash
sscha -i input_file>
```
Where the "input_file" is a valid input for the code. Please look at the SimpleExample
inside the Examples directory to see a template and how to create it.

Alternative the code can be used directly as a python library (expert mode).
This will allow much more customization on the input and the behaviour of the
code, and the easy interface with calculators different from Quantum ESPRESSO.

## ENJOY!
