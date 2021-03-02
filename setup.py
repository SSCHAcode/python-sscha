from __future__ import print_function

from numpy.distutils.core import setup, Extension
import os, sys 
import numpy 

extra_flags_c = ["-fopenmp"]
extra_link_args_c = ["-fopenmp"]
mpi_compile_args = []
mpi_link_args = []

# If true, do not check for parallel
avoid_parallel_test = True 

# Get the MPI from environmental variables
parallel = False 
if "MPICC" in  os.environ:
        mpicc = os.environ["MPICC"]
        parallel = True
        print ()
        print("Parallel compiler setted to:", mpicc)
        print()

# Check for the python parallel libraries
python_parallel = True
try:
        import pypar
except:
        try:
                import mpi4py
        except:
                #parallel = False
                python_parallel = False

# Setup the parallel environemnt
if parallel:
        # If we are here we can compile using MPI support
        mpi_compile_args = os.popen("%s -show" % mpicc).read().strip().split(' ')[1:]
        mpi_link_args    = os.popen("%s -show" % mpicc).read().strip().split(' ')[1:]
        extra_flags_c += ["-D_MPI"]

# Check if it is python2 or 3
if sys.version_info[0] < 3:
        print("Running on python 2, added the flag -D_PYTHON2")
        extra_flags_c += ["-D_PYTHON2"]




LIBRARIES = ["lapack", "blas"]
EXTRA_F90_FLAGS =  ["-cpp", "-fopenmp"]
EXTRA_LINK_ARGS = ["-fopenmp"]

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
                                   "SCHAModules/get_gradient_supercell_fast.f90",
                                   "SCHAModules/get_g.f90",
                                   "SCHAModules/get_emat.f90",
                                   "SCHAModules/get_v3.f90",
                                   "SCHAModules/get_odd_straight.f90",
                                   "SCHAModules/get_cmat.f90",
                                   "SCHAModules/get_v4.f90",
                                   "SCHAModules/get_odd_straight_with_v4.f90"],
                        libraries = LIBRARIES,
                        extra_f90_compile_args = EXTRA_F90_FLAGS,
                        extra_link_args= EXTRA_LINK_ARGS)




# Setup the HP performance module
#odd_HP = Extension(name = "sscha_HP_odd",
#                   sources= ["CModules/odd_corr_module.c", "CModules/LanczosFunctions.c"],
#                   include_dirs=[numpy.get_include()],
#                   extra_compile_args= extra_flags_c + mpi_compile_args,
#                   extra_link_args = mpi_link_args + extra_link_args_c
#                   )



# Prepare the compilation of the Python Conde
setup( name = "python-sscha",
       version = "1.0alpha4",
       description = "Python implementation of the sscha code",
       author = "Lorenzo Monacelli",
       url = "https://github.com/mesonepigreco/python-sscha",
       packages = ["sscha"],
       package_dir = {"sscha": "Modules"},
       install_requires = ["numpy", "ase", "scipy", "cellconstructor", "spglib", "matplotlib"],
       ext_modules = [SCHAModules], # odd_HP
       scripts = ["scripts/sscha", "scripts/cluster_check.x", "scripts/plot_frequencies.py",
                  "scripts/static-vc-relax.pyx", "scripts/read_incomplete_ensemble.py",
                  "scripts/plot_lanczos_convergence.py"],
       license = "GPLv3"
       )
                                                                                                                                                          
if not python_parallel and not parallel and not avoid_parallel_test:                                                                                                                      
        print()                                                                                                                                               
        print("======= WARNING =======")                                                                                                                      
        print("Nor python parallel neither MPI compiler found.")                                                                                              
        print("If you whish to activate MPI acceleration,")                                                                                                   
        print("Consider installing either pypar or mpi4py")                                                                                                   
        print("For example, try to run: ")                                                                                                                    
        print(" >>> MPICC=mpicc python " + " ".join(sys.argv))                                                                                                
        print("Note: clean the build directory if you whish to recompile the code.")                                                                          
        print("=======================")                                                                                                                      
        print()                                                                                                                                               
elif not parallel and not avoid_parallel_test:
        print()
        print("======= WARNING =======")
        print("No MPI compiler found, please specify MPICC environmental variable")
        print("For example, try to run: ")
        print(" >>> MPICC=mpicc python " + " ".join(sys.argv))
        print("Note: clean the build directory if you whish to recompile the code.")
        print("=======================")
        print()
elif not python_parallel and not avoid_parallel_test:
        print()
        print("======= WARNING =======")
        print("No Python MPI library found")
        print("Supported libraries:")
        print(" - pypar ")
        print(" - mpi4py ")
        print()
        print("Note: Fast MPI implemetation will crash if used")
        print("      consider to install one of these libraries.")
        print("      (No need to reinstall python-sscha)")
        print("=======================")
        print()
elif not avoid_parallel_test:
        print()
        print(" PARALLEL ENVIRONMENT DETECTED CORRECTLY ")
        print()



def readme():
    with open("README.md") as f:
        return f.read()
