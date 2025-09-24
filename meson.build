project('python-sscha',
  ['c','fortran'],
  version : '1.4.2',
  license: 'GPL',
  meson_version: '>= 1.1.0', # <- set min version of meson.
  default_options : [
        'warning_level=1',
        'buildtype=release',
        'fortran_args=-O2',
        'fortran_args=-cpp'
    ])

  # --- System and Python Dependencies ---
  # Imports the Meson Python module to handle extensions and Python installation.
  # --- Opciones de compilación ---
  use_mkl = get_option('use_mkl')

  # --- System and Python Dependencies ---
  # Find the necessary installations
  py = import('python').find_installation(pure: false)
  py_dep = py.dependency()
  fc = meson.get_compiler('fortran')

  # Additional dependencies

  # --- BLAS/LAPACK Dependencies ---
  # Meson will attempt to automatically detect LAPACK and BLAS.
  # You don't need 'lapack_opt' as in numpy.distutils.
  # If detection fails, you can pass configuration options to Meson
  # (e.g., -Dblas_args='-L/path/to/blas -lblas').
  openblas_dep = dependency('openblas', required: false)
  if use_mkl
    mkl_dep = dependency('mkl', required: true)
    blas_dep = mkl_dep
    lapack_dep = mkl_dep
  else
    lapack_dep = dependency('lapack', required: true)
    blas_dep = dependency('blas', required: true)
  endif

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
  # Find the quadmath library, if available
  quadmath_dep = fc.find_library('quadmath', required: false)

  # --- NUMPY CONFIGURATION ---
  # Gets the path to the NumPy and f2py header directories using the Python command
  incdir_numpy = run_command(py,
      ['-c', 'import numpy; print(numpy.get_include())'],
      check : true
  ).stdout().strip()

  # f2py also requires the fortranobject.h header
  incdir_f2py = run_command(py,
      ['-c', 'import numpy.f2py; print(numpy.f2py.get_include())'],
      check : true
  ).stdout().strip()

  inc_np = include_directories(incdir_numpy, incdir_f2py)
  np_dep = declare_dependency(include_directories: inc_np)

  inc_f2py = include_directories(incdir_f2py)
  fortranobject_c = incdir_f2py / 'fortranobject.c'

  # --- END OF NUMPY CONFIGURATION ---
  if openblas_dep.found()
    message('openblas environment detected correctly.')
    list_dep = [py_dep, mpi_dep, quadmath_dep, openblas_dep, lapack_dep, blas_dep]
  else
    warning('No openblas found.')
    list_dep = [py_dep, mpi_dep, quadmath_dep, lapack_dep, blas_dep]
  endif

# --- Common Flags ---
common_fortran_flags = ['-cpp', '-fopenmp']
common_link_flags = ['-fopenmp'] # For OpenMP

# --- SCHAModules Module (Fortran) ---
# This compiles and links your Fortran module with the Python, LAPACK, BLAS, and MPI dependencies
# Note: The Fortran extension will generate a .so or .pyd file that Python can import.

# --- SCHAModules extension ---

schamodules_sources = [
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
  'SCHAModules/get_odd_straight_with_v4.f90'
  # ,'SCHAModules/module_polarization.f90'
]
## Generate the C wrapper with f2py using a custom target
f2py_SCHAModules_target = custom_target('SCHAModules-f2py-wrapper',
    input: schamodules_sources,
    output: ['SCHAModulesmodule.c','SCHAModules-f2pywrappers.f','SCHAModules-f2pywrappers2.f90'],
#    command: [py.full_path(), '-m', 'numpy.f2py', '--backend', 'meson',
#        '--dep', 'mpi', '--quiet', '-c', '@INPUT0@', '-m', 'SCHAModules'],
    command : [py, '-m', 'numpy.f2py', '@INPUT@', '-m', 'SCHAModules', '--lower'],
    install: false
)
py.extension_module('SCHAModules',
   [
    schamodules_sources,
     f2py_SCHAModules_target, fortranobject_c],
    include_directories: inc_np,
    dependencies: list_dep,
    install: true
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
# Define la ruta de la librería
#lib_dir = 'Modules'
#lib_name = 'Calculator.py'
#lib_path = python.find_sources(lib_dir, lib_name)

# Instala la librería en el sitio de paquetes de Python
#python.install_sources(
#    ['Modules/Calculator.py'],
#    pure: true,
#    subdir: 'sscha'
#)

# install_subdir('Modules', install_dir : join_paths(get_option('prefix'), 'share', 'python-sscha', 'Modules'))
# install_subdir('Modules', install_dir : 'share', strip_directory : false)   <-- Not the best option
# install_subdir('Modules', install_dir : py.get_install_dir()) <-- .get_install_dir() gives the python installation place for libraries

py.install_sources([
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
  'Modules/Utilities.py'
  ],
  subdir: 'sscha',
)

# --- Installing Scripts ---
# Meson is great for installing scripts and making them executable.
# Scripts will be installed in the `bin` directory of the Python environment (e.g., venv/bin/).

py.install_sources([
  'scripts/sscha.x',
  'scripts/cluster_check.x',
  'scripts/plot_frequencies.py',
  'scripts/sscha-plot-data.py',
  'scripts/static-vc-relax.pyx',
  'scripts/read_incomplete_ensemble.py'
])

# --- Copy package data (e.g. *.jl) ---
# This emulates package_data={"": ["*.jl"]}
# Create a 'sscha' subdirectory within the Python installation directory
# and copy the .jl files there.
install_data(
  ['Modules/fourier_gradient.jl'], # List the .jl files you need
  install_dir : py.get_install_dir() / 'sscha'
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


# Set the tests by pytest.
pytest_exe = find_program('pytest', required: false)

if pytest_exe.found()
  test('pytest', pytest_exe,
       args : ['-v', '-m', 'not release', 'tests'],
       workdir : meson.project_source_root())
else
  message('pytest not found; pytest tests are skipped.')
endif
