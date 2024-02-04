from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
import os
import sys
import subprocess

class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)

class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError('CMake must be installed to build the following extensions: ' +
                               ', '.join(e.name for e in self.extensions))

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        cmake_args = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
                      '-DPYTHON_EXECUTABLE=' + sys.executable]

        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]

        # Assuming CMakeLists.txt is in the same directory as setup.py
        os.chdir(ext.sourcedir)
        self.spawn(['cmake', '.'] + cmake_args)
        if not self.dry_run:
            self.spawn(['cmake', '--build', '.'] + build_args)
        os.chdir(self.get_ext_fullpath(ext.name))

setup(
    name='python-sscha',
    version='1.4.0',
    author='Lorenzo Monacelli',
    description='Python implementation of the sscha code',
    long_description='',
    ext_modules=[CMakeExtension('SCHAModules', sourcedir='')],
    cmdclass=dict(build_ext=CMakeBuild),
    packages=find_packages(),
    package_dir={'sscha': 'Modules'},
    include_package_data=True,
    package_data={"": ["*.jl"]},
    install_requires=['numpy', 'ase', 'scipy', 'cellconstructor', 'spglib', 'matplotlib'],
    scripts=['scripts/sscha', 'scripts/cluster_check.x', 'scripts/plot_frequencies.py',
             'scripts/sscha-plot-data.py', 'scripts/static-vc-relax.pyx', 'scripts/read_incomplete_ensemble.py'],
    license='GPLv3',
    zip_safe=False
)

