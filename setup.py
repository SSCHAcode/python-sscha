from distutils.core import setup
setup( name = "python-sscha",
       version = "0.1",
       description = "Python implementation of the sscha code",
       author = "Lorenzo Monacelli",
       url = "https://github.com/mesonepigreco/python-sscha",
       packages = ["python-sscha"],
       package_dir = {"python-sscha": "Modules"},
       install_requires = ["numpy", "ase", "scipy", "cellconstructor"],
       license = "GPLv3"
       )

def readme():
    with open("README.md") as f:
        return f.read()