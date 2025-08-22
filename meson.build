project('python-sscha', 'fortran', 'c',
  version : '1.4.2',
  license: 'GPL',
# Common compilation options for Fortran.
# -cpp is essential for the Fortran preprocessor.
# -O2 is a general optimization you had in setup.py.
  default_options : ['warning_level=1','fortran_args=-cpp', 'fortran_args=-O2'])

  # --- System and Python Dependencies ---
  # Imports the Meson Python module to handle extensions and Python installation.
  python = import('python')
  python_installation = python.find_installation()
  numpy_incdir = run_command(
    python_installation, '-c', 'import numpy; print(numpy.get_include())',
    check: true
  ).stdout().strip()

# --- BLAS/LAPACK Dependencies ---
# Meson will attempt to automatically detect LAPACK and BLAS.
# You don't need 'lapack_opt' as in numpy.distutils.
# If detection fails, you can pass configuration options to Meson
# (e.g., -Dblas_args='-L/path/to/blas -lblas').
lapack_dep = dependency('lapack', required : true)
blas_dep = dependency('blas', required : true)

# SSCHA needs CellConstructor
CellConstructor_dep = dependency('cellconstructor', required : true)

# Look for the OpenMP library if it is needed for parallelization.
openmp_dep = dependency('openmp', required: true)

# --- MPI Detection ---
# This overrides the logic in os.environ["MPICC"] and os.popen("%s -show" % mpicc).
# Meson has a built-in MPI module.
mpi_args = []
mpi_link_args = []
mpi_compile_args = []
has_mpi = false

# Attempts to find the MPI dependency.
# You can specify a specific MPI compiler with the 'mpi_compiler' parameter
# or, if you want, a specific Fortran compiler with 'mpi_fortran_compiler'.
# For OpenMPI, IntelMPI, MPICH, etc., Meson usually finds it automatically.
#mpi_dep = dependency('mpi', required: false, language: ['c', 'fortran'])
mpi_dep = dependency('mpi', required: false)

if mpi_dep.found()
  message('MPI environment detected correctly.')
  has_mpi = true
  # Meson handles adding appropriate flags. We just add the define.
  # If you need specific MPI flags beyond what Meson adds automatically,
  # you can get them via mpi_dep.get_compile_args() and mpi_dep.get_link_args()
  # and add them to extra_compile_args/extra_link_args.
  mpi_compile_args += ['-D_MPI']
else
  # Here you can add warning logic if MPI is not found.
  # Meson prints a warning if required: true and it is not found.
  # For required: false, you can print your own warning.
  warning('No MPI compiler found, please ensure MPI is installed and configured.')
  warning('If you wish to activate MPI acceleration, consider setting MPICC environment variable or providing Meson with appropriate flags.')
endif

# --- Python 2/3 Detection (not as relevant anymore, but translated) ---
# In Meson, this would be done based on the Python version Meson is using.
# Since Python 2 is EOL, this check is generally not needed.
# However, if you want to keep the flag:
#py2_flag = []
#if py_dep.version().version_compare('<3.0')
#  message('Running on Python 2, adding the flag -D_PYTHON2')
#  py2_flag += ['-D_PYTHON2']
#endif

# --- Common Flags ---
common_fortran_flags = ['-cpp', '-fopenmp']
common_link_flags = ['-fopenmp'] # For OpenMP

# --- SCHAModules Module (Fortran) ---
# This compiles and links your Fortran module with the Python, LAPACK, BLAS, and MPI dependencies
# Note: The Fortran extension will generate a .so or .pyd file that Python can import.
schamodules_sources = files(
  'SCHAModules/module_stochastic.f90',
  'SCHAModules/module_new_thermodynamic.f90',
  'SCHAModules/module_anharmonic.f90',
  'SCHAModules/get_stress_tensor.f90',
  'SCHAModules/get_gradient_supercell.f90',
  'SCHAModules/get_upsilon_matrix.f90',
  'SCHAModules/multiply_lambda_tensor.f90',
  'SCHAModules/cell_force.f90',
  'SCHAModules/get_gradient_supercell_fast.f90',
  'SCHAModules/get_g.f90',
  'SCHAModules/get_emat.f90',
  'SCHAModules/get_v3.f90',
  'SCHAModules/get_odd_straight.f90',
  'SCHAModules/get_cmat.f90',
  'SCHAModules/get_v4.f90',
  'SCHAModules/get_odd_straight_with_v4.f90',
  'SCHAModules/module_polarization.f90'
)

# Meson's `extension_module` function for Python handles much of the magic.
# 'SCHAModules' will be the name of the resulting Python module.
py_mod_schamodules = python_installation.extension_module(
  'SCHAModules',
  schamodules_sources,
  # Make sure BLAS/LAPACK dependencies are linked correctly
  dependencies : [lapack_dep, blas_dep, mpi_dep, CellConstructor_dep],
  # Fortran-specific flags
  fortran_args : common_fortran_flags,
  link_args : common_link_flags,
  # Add -fPIC for Linux if needed (Meson usually handles this)
  # This is the equivalent of extra_f90_compile_args and extra_link_args
  # in setup.py for the Fortran module.
  include_directories: include_directories(numpy_incdir),
)

# --- Module odd_HP (C) ---
# Uncomment and adapt if you want to include this module.
# This would be done similarly to the Fortran module.
# odd_hp_sources = files(
#   'CModules/odd_corr_module.c',
#   'CModules/LanczosFunctions.c'
# )
# py_mod_odd_hp = python.extension_module(
#   'sscha_HP_odd', # Make sure this is the desired import name
#   odd_hp_sources,
#   dependencies : [py_dep, mpi_dep], # assumes it only depends on Python and MPI
#   c_args : py2_flag + mpi_compile_args + ['-O3'], # Extra compile args for C
#   link_args : mpi_link_args # Extra link args for C
# )


# --- Installing Python Packages ---
# The 'Modules' directory will be installed as the 'sscha' package.
#python.install_sources(
#  'Modules/', # Path to your Python package directory
# The sources in this directory will be copied.
# You can also use install_dir: python.get_install_dir() / 'sscha'
# If you want to install it in a specific subdirectory of site-packages.
# By default, Meson-Python will install it correctly.
#  subdir: 'sscha'
#)
install_data(
  'Modules/__init__.py',
  'Modules/AdvancedCalculations.py',
  'Modules/aiida_ensemble.py',
  'Modules/Calculator.py',
  'Modules/Cluster.py',
  'Modules/Dynamical.py',
  'Modules/Ensemble.py',
#  'Modules/fourier_gradient.jl',
  'Modules/LocalCluster.py',
  'Modules/Minimizer.py',
  'Modules/Optimizer.py',
  'Modules/Parallel.py',
  'Modules/Relax.py',
  'Modules/SchaMinimizer.py',
  'Modules/Tools.py',
  'Modules/Utilities.py',
  install_dir: python_installation.get_install_dir() / 'sscha',
)

# --- Installing Scripts ---
# Meson is great for installing scripts and making them executable.
# Scripts will be installed in the `bin` directory of the Python environment (e.g., venv/bin/).
#executable(
#  'scripts/sscha',
#  'scripts/cluster_check.x',
#  'scripts/plot_frequencies.py',
#  'scripts/sscha-plot-data.py',
#  'scripts/static-vc-relax.pyx', # Si esto es un script ejecutable directamente
#  'scripts/read_incomplete_ensemble.py'
#, install : true)
#install_data = (
#  [
#  'scripts/sscha',
#  'scripts/cluster_check.x',
#  'scripts/plot_frequencies.py',
#  'scripts/sscha-plot-data.py',
#  'scripts/static-vc-relax.pyx')#, # Si esto es un script ejecutable directamente
#  'scripts/read_incomplete_ensemble.py'
#  ],
#  install_dir: get_option('bindir')
#  )

#foreach s : install_scripts
#  install_data(s, install_dir : python.get_install_subdir() / '..' / 'bin') # Installs in bin
#endforeach

# --- Copy package data (e.g. *.jl) ---
# This emulates package_data={"": ["*.jl"]}
# Create a 'sscha' subdirectory within the Python installation directory
# and copy the .jl files there.
install_data(
  ['Modules/fourier_gradient.jl'], # List the .jl files you need
  install_dir : python_installation.get_install_dir() / 'sscha'
)
# If there are many .jl files in multiple subdirectories,
# you would need to use 'install_subdir' or list them all explicitly.
# For a more general `package_data` approach for non-Python files,
# it is often recommended to use the `package_data` from `pyproject.toml` or
# an `sdist` that includes those files and then the `wheel` packages them.
# Meson focuses more on building and compiling.

# You can use configure_file to generate files if needed,
# for example for dynamic versions or generated data.
# configure_file(input: 'src/my_template.py.in',
# output: 'my_module/version.py',
# configuration: conf)

# Define el intérprete de Python
python = find_program('python3')

# Crea un array con los nombres de todos los scripts que quieres instalar
python_scripts = [
#  'scripts/sscha',
#  'scripts/cluster_check.x',
#  'scripts/plot_frequencies.py',
#  'scripts/sscha-plot-data.py',
#  'scripts/static-vc-relax.pyx',
#  'scripts/read_incomplete_ensemble.py'
]

# Itera sobre el array y crea un ejecutable para cada script
foreach script : python_scripts
  # Obtiene el nombre del archivo de la ruta, por ejemplo 'cluster_check.x'
  script_basename = script.split('/')[-1]

  # Divide el nombre por el punto y toma el primer elemento para quitar la extensión
  # Si el nombre no tiene punto, esto lo maneja correctamente
  script_name_without_ext = script_basename.split('.')[0]

  executable(
    script_name_without_ext,
#    python,
    script,
    install: true
  )
endforeach
# executable('sscha','scripts/sscha',install: true)
# executable('cluster_check','scripts/cluster_check.x',install: true)
# executable('plot_frequencies','scripts/plot_frequencies.py',install: true)
# executable('sscha-plot-data','scripts/sscha-plot-data.py',install: true)
# executable('static-vc-relax','scripts/static-vc-relax.pyx', install: true)
# executable('read_incomplete_ensemble','scripts/read_incomplete_ensemble.py',install: true)
install_data([
  'scripts/sscha',
  'scripts/cluster_check.x',
  'scripts/plot_frequencies.py',
  'scripts/sscha-plot-data.py',
  'scripts/static-vc-relax.pyx',
  'scripts/read_incomplete_ensemble.py'
],
install_dir: get_option('bindir')
)
