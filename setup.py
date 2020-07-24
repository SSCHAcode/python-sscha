from __future__ import print_function
from numpy.distutils.core import setup, Extension




LIBRARIES = ["lapack", "blas"]
EXTRA_F90_FLAGS =  ["-cpp", "-fopenmp"]
EXTRA_LINK_ARGS = ["-fopenmp"]

# Compile the fortran SCHA modules
SCHAModules = Extension(name = "SCHAModules", 
                        sources = ["SCHAModules/module_stochastic.f90",
                                   "SCHAModules/module_new_thermodynamic.f90",
                                   "SCHAModules/module_anharmonic.f90",
				   "SCHAModules/module_harmonic.f90",
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


# Prepare the compilation of the Python Conde
setup( name = "python-sscha",
       version = "0.1",
       description = "Python implementation of the sscha code",
       author = "Lorenzo Monacelli",
       url = "https://github.com/mesonepigreco/python-sscha",
       packages = ["sscha"],
       package_dir = {"sscha": "Modules"},
       install_requires = ["numpy", "ase", "scipy", "cellconstructor", "lapack", "blas"],
       ext_modules = [SCHAModules],
       scripts = ["scripts/sscha", "scripts/cluster_check.x", "scripts/plot_frequencies.pyx",
                  "scripts/static-vc-relax.pyx"],
       license = "GPLv3"
       )

def readme():
    with open("README.md") as f:
        return f.read()
