How to install
==============


Requirements
------------

The requirements of the python-sscha package are:
1. python >= 2.7 and < 3
2. numpy
3. matplotlib
3. Lapack
4. Blas
5. gfortran (or any fortran compiler)
6. CellConstructor

Recommended packages:
1. Atomic Simulation Environment (ASE)
2. SPGLIB
   
These packages are fundamental. In particular CellConstructor is the shoulder on
which python-sscha is builded on. You can find the last development version on
GitHub

The ASE (Atomic Simulation Environment) is another dependency, since CellConstructor
relies on it.

To install all the dependences, simply run

```bash
pip install -r requirements.txt
```



Installation
------------

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

Please, remember that if you use the intel compiler, you need to delete the lapack linking from the
setup.py and include the -mkl (as done for cellconstructor).
Note that you must force to use the same liker compiler as the one used for the compilation.

A specific setup.py script is provided to install it easily in FOSS clusters.


