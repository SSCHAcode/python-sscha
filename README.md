# SSCHA

The stochastic self-consistent harmonic approximation (SSCHA) is a full computational Python package that simulates thermodynamic and phononic properties of materials, accounting for anharmonicity at a nonperturbative level, fully including quantum and thermal fluctuations.

See more info on the webpage:

[www.sscha.eu](https://sscha.eu)

## Easy installation through Anaconda/Mamba

The SSCHA code is available as a Python library, with computationally intensive parts accelerated using C, Fortran, and Julia. 
A quick way to install all the requirements to compile the code in a separate environment is to use `Anaconda`.
Alternatively, when `anaconda` is too slow, we recommend using `mamba` or `micromamba`. Just replace `coda` with `mamba` or `micromamba`, respectively.

```
conda create -n sscha -c conda-forge python=3.12 gfortran libblas lapack openmpi julia openmpi-mpicc pip numpy scipy spglib=2.2
conda activate sscha
pip install ase julia mpi4py
pip install cellconstructor python-sscha tdscha
```

This is the safest and best way to install the SSCHA. The first line creates a new pristine Python environment with all the required libraries to compile the source code. The second line activates the newly installed environment. Then, the third command installs the additional dependencies, and the last line compiles and installs the SSCHA code.

Note, the first time you start up a sscha calculation, the code will try to download extra packages to set up the python
julia interface. This process may fail if you do not have an internet connection available or if the Julia installation failed.
Note that this is not mandatory, as the code will fall back to the old Fortran implementation (prior to version 1.4) and continue to run.
To achieve the Julia speedup in such a case, please take a look at the section on Manual Installation to preconfigure your system correctly.

## Video lessons from  the 2023 School are available

The full recordings, both of theoretical lectures, tutorials and Hands-on sessions can be found
in our youtube channel [SSCHAcode](https://www.youtube.com/@SSCHAcode>)

To use the SSCHA, you must activate the python environment with (replace `conda` with `mamba` or `micromamba` if you use those environment managers):

```
conda activate sscha
```

This installation method should also work on clusters and with computers with custom configurations. You must remember to activate the ``sscha`` environment even in your submission scripts on clusters.

To activate the Julia speedup for SSCHA minimization, ensure that Julia dependencies are correctly set up. To do this, run the following line:

```
python -c 'import julia; julia.install()'
```


## Installing without Anaconda

If you do not have Anaconda to handle your dependencies, you need to compile the code manually.

Most of the codes require a Fortran or C compiler and an MPI configured. Here we install all the requirements to set up the SSCHA code correctly. To properly compile and install the SSCHA code, you need a Fortran compiler and LAPACK/BLAS available.

On a Debian-based Linux distribution, all the software required is installed with (Tested on Quantum Mobile and Ubuntu 20.04):
```
sudo apt update
sudo apt install libblas-dev liblapack-dev liblapacke-dev gfortran openmpi-bin
```
Note that some of the names of the libraries may change slightly in different Linux versions or on macOS.

### Python installation

Up to version 1.4 of SSCHA, it supports only Python <= 3.10 and the old setuptools<=64 (we require the deprecated distutils). Starting from version 1.5, we dropped distutils and now support more recent version of python. 
If you are using the default Python in the system, make sure to have installed the development header files. On Ubuntu, they can be installed with:

```
sudo apt install python-dev
```

If you use anaconda (or mamba/micromamba), they are automatically installed.

### Prerequisites

The SSCHA code comprises three Python packages: CellConstructor, python-sscha, and tdscha.

- [CellConstructor](https://github.com/SSCHAcode/CellConstructor>): utility to manage phonon dispersions, atomic structures and crystal symmetries
- [sscha](https://github.com/SSCHAcode/python-sscha>) : This repository, relax the structure with anharmonicity and computes static linear response properties.
- [tdscha](<https://github.com/SSCHAcode/tdscha>) : Compute the dynamical linear response (Raman and IR, spectral functions)

More details about installations are on the official website [www.sscha.eu](https://sscha.eu/download>)

## Install with Anaconda/Mamba

The easiest way to install the code is through anaconda/mamba.

The following commands are sufficient to install the full sscha suite and its dependencies.
Replace `conda` with `mamba` or `micromamba` if that is the package manager you are using.

```
conda create -n sscha -c conda-forge python=3.12 gfortran libblas lapack openmpi julia openmpi-mpicc pip numpy scipy spglib=2.2
conda activate sscha
pip install ase julia mpi4py
pip install cellconstructor python-sscha tdscha
```

If you happen to get an error when using Julia, please try to install Julia from the official website and see the passages reported in the Manual installation.

To activate the environment and execute the SSCHA, run

```
   conda activate sscha
```


## Manual installation

The SSCHA benefits from Julia being installed in the system. If present,
it will be automatically used to speedup the calculation.

To install Julia, refer to the official website [julialang.org/downloads/](https://julialang.org/downloads/)
Alternatively, to install Julia on Linux we can employ JuliaUp:

```
  curl -fsSL https://install.julialang.org | sh
```

Hit enter when asked to install Julia.

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

This should install the required libraries. Press backspace to return to the standard Julia prompt and exit with

```
  julia> exit()
```

Now, you should be able to exploit the Julia speedup in the TDSCHA calculations. It is not required to install Julia before TDSCHA; it can also be done at a later time.


### Compiling SSCHA

Once the prerequisites have been installed, python-sscha can be downloaded and installed with

```
  pip install cellconstructor python-sscha
```

Alternatively, it is possible to use the most recent version from the [SSCHA GitHub](https://github.com/SSCHAcode) repository, under CellConstructor and python-sscha repositories. The installation is performed in this case with


```
  pip install .
```

### Personalize the compiler

If you have multiple compilers installed, and want to force pip to employ a specific fortran compiler, you can specify its path in the FC environment variable. Remember that the compiler employed to compile the code should match with the linker, indicated in the LDSHARED variable.

For example

```
  FC=gfortran LDSHARED=gfortran pip install cellconstructor python-sscha
```

For the development version of the code.

## Compiling with Meson

To compile and install SSCHA with Meson, follow these typical steps:

### 1. Change to the Source Directory

First, you can open a terminal and then navigate to the root directory of the project source code. This is where the `meson.build` file is located.

```bash
cd /path/to/source/root/python-sscha
```


### 2. Configure the Build Directory

Create and configure a build directory by running:

```bash
meson setup builddir
```

or if you are in a conda env (the best option for a local installation):
```bash
meson setup builddir --prefix=$CONDA_PREFIX
```

if you want to use Intel MKL:
```bash
setup builddir -Duse_mkl=true
```

This command sets up a separate build directory (`builddir`) where all compiled files and build artifacts will be placed, keeping the source directory clean. After this, change into the build directory:

```bash
cd builddir
```


### 3. Compile the Project

Once inside the build directory, compile the project using:

```bash
meson compile
```

This will compile the source code according to the configuration from the previous step.

### 4. Run Tests (Optional)

The project includes tests, you need to install pytest to work. You can run them with:

```bash
meson test
```

This step helps verify that the build works correctly.

### 5. Install the Project (Optional)

To install the compiled binaries, libraries, and other files system-wide (or to a custom location), run:

```bash
meson install
```

or

```bash
sudo meson install
```

You may need superuser privileges (hence `sudo`) to install to system directories.

***

Following these steps will help you successfully compile, test, and install SSCHA with Meson as their build system.
