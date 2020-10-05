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
3. LD_SHARED (linking)

If we want to use a custom compiler in /path/to/fcompiler we may run the setup as:

.. code-block:: console

   FC=/path/to/fcompiler LD_SHARED=/path/to/fcompiler python setup.py install



A specific setup.py script is provided to install it easily in FOSS clusters.


