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

Meson works automatically with several FORTRAN compilers, including Intel FORTRAN. However, due to some differences in linking lapack,
to use the intel compiler you need to:

Ensure MKL is installed in your Conda environment:
.. code-block:: console

   conda install mkl mkl-devel

You can pass Meson options through pip's \--config-settings flag.
.. code-block:: console

    pip install . --config-settings=--setup-args=-Duse_mkl=true

Or for an editable install:
.. code-block:: console
    pip install -e . --config-settings=--setup-args=-Duse_mkl=true

Install with a specific compiler path
-------------------------------------

Method 1: Using Environment Variables (Recommended for most cases)
------------------------------------------------------------------

Meson automatically detects compilers using standard environment variables. You can set these variables before running the installation command. This is the simplest way to specify a compiler for a single build.

The key variables are:

    CC: Specifies the C compiler executable.

    CXX: Specifies the C++ compiler executable.

    FC: Specifies the Fortran compiler executable.

Step-by-Step Instructions

1. Open your terminal. All commands must be run in the same session, as environment variables are typically not permanent.
2. Set the environment variables to point to your desired compilers.

Example for C (using a specific gcc):

.. code-block:: console
    export CC=/path/to/my/custom/gcc

Example for Fortran (using a specific gfortran):

.. code-block:: console
    export FC=/path/to/my/custom/gfortran

Example for C++ (if the project needed it):

.. code-block:: console
    export CXX=/path/to/my/custom/g++

3. Combine them as needed. For this project, you will likely need to set CC and FC.

# Example using compilers from a specific toolchain
.. code-block:: console
    export CC=/usr/local/bin/gcc-11
    export FC=/usr/local/bin/gfortran-11

4. Run the pip installation. With the variables set, run pip install from the project's root directory. pip will pass the environment variables down to Meson.

.. code-block:: console
    # Ensure you are in the project's root directory (where pyproject.toml is)
    pip install .

The build process will now use the compilers you specified instead of the system defaults.

Method 2: Using a Meson Cross File (Advanced & Reproducible)
------------------------------------------------------------

For a more permanent or reproducible setup (e.g., in CI/CD pipelines or complex environments), a Meson "cross file" is the best practice. This file explicitly defines the toolchain.
Step-by-Step Instructions

1. Create a cross file. In your project's root directory, create a file named native-toolchain.ini (the name can be anything).

2. Edit the file to specify the paths to your compilers in the [binaries] section.

Example native-toolchain.ini:

.. code-block:: console
    # native-toolchain.ini
    # Defines the compilers to use for the build.

    [binaries]
    c = '/path/to/my/custom/gcc'
    fortran = '/path/to/my/custom/gfortran'

    # If you also needed C++, you would add:
    # cpp = '/path/to/my/custom/g++'

Note: The keys are c, cpp, and fortran.

3. Run the pip installation with meson-args. You can instruct pip to pass configuration settings to meson-python, which in turn passes them to Meson. We use this to specify our cross file.

.. code-block:: console
    # The --native-file option tells Meson which toolchain to use.
    pip install . --config-settings=meson-args="--native-file=native-toolchain.ini"

This method is more explicit and less prone to errors from shell configurations. It ensures that anyone building the project can easily use the exact same toolchain by sharing the native-toolchain.ini file.
