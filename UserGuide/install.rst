How to install
==============

The SSCHA code is a collection of 2 python packages: CellConstructor and python-sscha.
In this guide, we refer to the installation of python-sscha.


Easy installation using Anaconda
--------------------------------

The easy way to install python-sscha is to use the anaconda distribution of python.

.. code-block:: console

    conda create -n sscha -c conda-forge python=3.11 gfortran libblas lapack openmpi julia openmpi-mpicc pip numpy scipy spglib
    conda activate sscha
    pip install ase julia mpi4py
    pip install cellconstructor python-sscha tdscha


This will create a new environment called sscha, install all the dependencies and the packages.
To use the code, you need to activate the environment:

.. code-block:: console

    conda activate sscha


The sscha code exploits the julia language to speed up the calculation.
To install the julia dependencies, you need to run the following command:

.. code-block:: console

   python -c 'import julia; julia.install()'


And that's it! You can now run the sscha code.

.. code-block:: console

   sscha -h


SSCHA also works as a library, so you can import it in your python scripts.
To check that everything is working, you can run the following script:

.. code-blocK:: python

   import sscha, sscha.Ensemble
   print("Hello world!")

If it runs without problems, the SSCHA code is working correctly.

Requirements
------------

To install python-sscha you need:
1. python<=3.11
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
Then run, inside the directory that contains the meson.build script, the following command:

.. code-block:: console

   pip install .


As for the pip installation, you may append the --user option to install the package only for the user (without requiring administrator powers).

An "editable" install is highly recommended for developers. It allows you to modify the source code and have the changes reflected immediately without needing to reinstall.
.. code-block:: console

   pip install -e .


Install with Intel FORTRAN compiler
-----------------------------------

Meson works automatically with the GNU FORTRAN compiler. However, due to some differences in linking lapack,
to use the intel compiler you need to:

Ensure MKL is installed in your Conda environment:
.. code-block:: console

   conda install mkl mkl-devel

You can pass Meson options through pip's \--config-settings flag.
.. code-block:: console

    pip install . --config-settings=--setup-args=-Duse_mkl=true

Or for an editable install:
.. code-block::
    pip install -e . --config-settings=--setup-args=-Duse_mkl=true

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
