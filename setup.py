from numpy.distutils.core import setup, Extension
import os, sys 
import numpy 

extra_flags_c = ["-fopenmp"]
extra_link_args_c = ["-fopenmp"]
mpi_compile_args = []
mpi_link_args = []

# Get the MPI from environmental variables
parallel = False 
if os.environ.has_key("MPICC"):
        mpicc = os.environ["MPICC"]
        parallel = True
else:
        print("No MPI compiler found, please specify MPICC environmental variable")

# Setup the parallel environemnt
if parallel:
        # If we are here we can compile using MPI support
        mpi_compile_args = os.popen("%s --showme:compile" % mpicc).read().strip().split(' ')
        mpi_link_args    = os.popen("%s --showme:link" % mpicc).read().strip().split(' ')
        extra_flags_c += ["-D_MPI"]


# Compile the fortran SCHA modules
SCHAModules = Extension(name = "SCHAModules", 
                        sources = ["SCHAModules/module_stochastic.f90",
                                   "SCHAModules/module_new_thermodynamic.f90",
                                   "SCHAModules/module_anharmonic.f90",
                                   "SCHAModules/get_stress_tensor.f90",
                                   "SCHAModules/get_gradient_supercell.f90",
                                   "SCHAModules/get_upsilon_matrix.f90",
                                   "SCHAModules/multiply_lambda_tensor.f90",
                                   "SCHAModules/cell_force.f90",
                                   "SCHAModules/get_gradient_supercell_fast.f90"],
                        libraries = ["lapack", "blas"],
                        extra_f90_compile_args = ["-cpp"])



# Setup the HP performance module
odd_HP = Extension(name = "sscha_HP_odd",
                   sources= ["CModules/odd_corr_module.c"],
                   include_dirs=[numpy.get_include()],
                   extra_compile_args= extra_flags_c + mpi_compile_args,
                   extra_link_args = mpi_link_args + extra_link_args_c
                   )



# Prepare the compilation of the Python Conde
setup( name = "python-sscha",
       version = "0.1",
       description = "Python implementation of the sscha code",
       author = "Lorenzo Monacelli",
       url = "https://github.com/mesonepigreco/python-sscha",
       packages = ["sscha"],
       package_dir = {"sscha": "Modules"},
       install_requires = ["numpy", "ase", "scipy", "cellconstructor", "lapack", "blas"],
       ext_modules = [SCHAModules, odd_HP],
       scripts = ["scripts/sscha", "scripts/cluster_check.x", "scripts/plot_frequencies.pyx",
                  "scripts/static-vc-relax.pyx", "scripts/read_incomplete_ensemble.py",
                  "scripts/plot_lanczos_convergence.py"],
       license = "GPLv3"
       )

def readme():
    with open("README.md") as f:
        return f.read()
