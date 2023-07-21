# SSCHA

The stochastic self-consistent harmonic approximation (SSCHA) is a full computational python package that simulates thermodynamic and phononic properties of materials accounting for anharmonicity at a nonperturbative level, fully including quantum and thermal fluctuations.

See more info on the webpage:

[www.sscha.eu](https://sscha.eu)


## Easy installation through Anaconda

The SSCHA code comes as a python library, with computationally intense part speedup with C, Fortran and Julia. The easiest way to install is through Anaconda ([how to install anaconda](https://www.anaconda.com/download)) 


```
conda create -n sscha -c conda-forge python=3.10 gfortran libblas lapack openmpi julia openmpi-mpicc pip numpy scipy spglib
conda activate sscha
pip install ase julia mpi4py
pip install cellconstructor python-sscha tdscha
```

This is the safest and best way to install the SSCHA. The first line creates a new pristine python environment with all the required libraries to compile the source code. The second line activates the newly installed environment. Then, the thrid command installs the additional dependencies, the last line compiles and install the SSCHA code.

To use the SSCHA, you must activate the python environment with:

```
conda activate sscha
```

This installation method should work also on clusters and with computers with custom configurations. You must remember to activate the ``sscha`` environment even in your submission scripts on clusters.

To activate the julia speedup on the SSCHA minimization, you must ensure julia dependencies are correctly setup. To do this, run the following line:

```
python -c 'import julia; julia.install()'
```


## Installing without Anaconda 

If you do not have anaconda to handle your dependencies you need to manually compile the code.

Most of the codes require a fortran or C compiler and MPI configured. Here we install all the requirements to properly setup the SSCHA code. To properly compile and install the SSCHA code, you need a fortran compiler and LAPACK/BLAS available.

On Debian-based Linux distribution, all the software required is installed with (Tested on Quantum Mobile and ubuntu 20.04):
```
sudo apt update
sudo apt install libblas-dev liblapack-dev liblapacke-dev gfortran openmpi-bin
```
Note that some of the names of the libraries may change slightly in different linux versions or on MacOS.

### Python installation

Up to version 1.4 of SSCHA, it supports only python <= 3.10. If you are using the default python in the system, make sure to have installed the development header files. On ubuntu, they can be installed with:

```
sudo apt install python-dev
```

If you use anaconda, they are automatically installed.

### Prerequisites

#### Python libraries
We need to install the prerequisites with pip:
```
pip install ase spglib 
```

#### Julia speedup
The SSCHA benefits from julia being installed in the system. If present,
it will be automatically used to speedup the calculation.

To install julia, refer to the official website [julialang.org/downloads/](https://julialang.org/downloads/)
Alternatively, to install julia on linux we can employ juliaup:
```
curl -fsSL https://install.julialang.org | sh
```
Hit enter when asked to install julia.

Then, install the python bindings for julia with
```
pip install julia
```


The tdscha extension to compute Raman and IR requires some additional julia packages that can be installed within a julia terminal. Update your configuration to have access to the newly installed julia
```
source ~/.bashrc
```
Then, open a terminal and type `julia`. Inside the julia prompt, type `]`. The prompt should change color and display the julia version ending with `pkg>`

Install the required julia libraries
```
pkg> add SparseArrays, LinearAlgebra, InteractiveUtils, PyCall
```
This should install the required libraries. Press backspace to return to the standard julia prompt and exit with
```
julia> exit()
```
Now, you should be able to exploit the julia speedup in the TDSCHA calculations. It is not required to install julia before TDSCHA, it can also be done in a later moment.


### Compiling SSCHA

Once the prerequisites have been installed, python-sscha can be downloaded and installed with
```
pip install cellconstructor python-sscha
```

Alternatively, it is possible to use the most recent version from the [SSCHA GitHub](https://github.com/SSCHAcode) repository, under CellConstructor and python-sscha repositories. The installation is performed in this case with

```
python setup.py install
```

#### Personalize the compiler

If you have multiple compilers installed, and want to force pip to employ a specific fortran compiler, you can specify its path in the FC environment variable. Remember that the compiler employed to compile the code should match with the linker, indicated in the LDSHARED variable.

For example
```
FC=gfortran LDSHARED=gfortran pip install cellconstructor python-sscha
```
For the development version of the code, subtitute the pip call with the python setup.py install.

