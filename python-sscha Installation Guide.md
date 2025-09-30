# python-sscha Installation Guide

This guide provides step-by-step instructions to compile and install the *sscha* library. The project uses the Meson build system to compile C and Fortran extensions for Python.

## Recommended Installation: Conda / Mamba

The easiest and most reproducible way to install *python-sscha* is by using a Conda environment. We strongly recommend using *micromamba* or *mamba* as they are significantly faster. This method installs all the required compilers, libraries, and Python packages inside an isolated environment, avoiding the need to modify your base system.

### Step 1: Create and Activate the Environment

1. Install Conda/Mamba. If you don't have it, we recommend installing *micromamba*.

2. Create the environment. Open your terminal and run the following command. It will create a new environment named *sscha* with all the necessary dependencies, including *CellConstructor*.

```bash
micromamba create -n sscha python=3.12 gfortran libblas lapack openmpi openmpi-mpicc pkg-config pip numpy scipy spglib=2.2 cellconstructor matplotlib ase
```

* Python Version: We use Python 3.12. Newer versions are not yet supported due to a dependency constraint from *spglib <= 2.2*.

* pkg-config: This package is essential. Meson requires it to correctly detect the BLAS and LAPACK libraries provided by Conda.

3. Activate the environment. You must activate the environment before proceeding.

```bash
micromamba activate sscha
```

### Step 2: Clone the Repository (if not done)

If you don't have the source code locally, clone it from the repository.

```bash
git clone https://github.com/SSCHAcode/python-sscha.git
cd python-sscha
```

### Step 3: Install python-sscha

With the environment active, install the package using *pip*. *pip* will automatically use Meson to compile and install the project.

```bash
pip install .
```

The installation is now complete! You can verify it by running the tests or importing the modules in Python.

## Advanced Installation: Manual Configuration

This section is for users who cannot use Conda or need to configure the build with specific compilers or options.

### Prerequisites

* A C and Fortran compiler.

* Python 3.12 and pip.

* Ninja and Meson: *pip install meson ninja*.

* System libraries for BLAS, LAPACK, and MPI.

* All Python dependencies, including CellConstructor, must be installed in your environment.

### Method 1: Using Environment Variables

You can specify compilers by setting the *CC* (C compiler) and *FC* (Fortran compiler) environment variables before running *pip*.

```bash
# Example using compilers from a specific toolchain
export CC=/usr/local/bin/gcc-11
export FC=/usr/local/bin/gfortran-11

# Install the project
pip install .
```

### Method 2: Using a Meson Cross File

For a reproducible setup, define your compilers in a file (e.g., *native-toolchain.ini*).

Example *native-toolchain.ini*:

```bash
[binaries]
c = '/path/to/my/custom/gcc'
fortran = '/path/to/my/custom/gfortran'
```

Then, install by passing this file to Meson via *pip*:

```bash
pip install . --config-settings=meson-args="--native-file=native-toolchain.ini"
```

## Build Options

You can pass options to Meson to customize the build.

* Create a debug build:

```bash
pip install . --config-settings=meson-args="-Dbuildtype=debug"
```

* Enable Intel MKL (requires MKL to be installed and findable by your system):

```bash
pip install . --config-settings=meson-args="-Duse_mkl=true"
```

* Combine multiple options:

```bash
pip install . --config-settings=meson-args="-Dbuildtype=debug,-Duse_mkl=true"
```

## Reconfiguring a Build

If you need to change build options, it is best to perform a clean installation.

```bash
# 1. Uninstall the package
pip uninstall python-sscha

# 2. (Optional but recommended) Remove the build directory
rm -rf build/

# 3. Reinstall with the new options
pip install . --config-settings=meson-args="-Dnew_option=value"
```
