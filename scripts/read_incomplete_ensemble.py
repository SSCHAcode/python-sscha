#!python
from  __future__ import print_function

try:
    import ase 
except:
    raise ImportError("Error, ase library needed to run this script.")

import sys, os
import numpy as np

import cellconstructor as CC 
import cellconstructor.Phonons

import sscha
import sscha.Ensemble 


# Create the units according to quantum-espresso
from ase.units import create_units
units = create_units("2006")#Rydberg, Bohr
Rydberg = units["Ry"]
Bohr = units["Bohr"]


INFO = """
This script reads an incomplete ensemble from a given working directory and saves it in the binary format

Usage: 
>>> read_incomplete_ensemble.py origin_dynmat PREFIX <save_dir> <population> [extension]

The script takes 4 or 5 positional arguments.
1) The original dynamical matrix. Only the prefix of the matrix is required.
The code will guess automatically how many irreducible q points must be read.
2) The calculation prefix. I.e. if we saved our espresso output files as 
data/ESP_X.pwo
where X is the number of the calculation, PREFIX must be 'data/ESP_' 
3) The directory that will store the ensemble.
4) The population index
5) [Optional] the extension for the calculation. (default = pwo, it is case sensitive)

For now it supports only quantum espresso calculation.
"""

if not len(sys.argv) in [5, 6]:
    print(INFO)
    print ("Error. Wrong number of arguments: %d", )
    sys.exit(1)

origin_dynmat = sys.argv[1]
PREFIX = sys.argv[2]
save_dir = sys.argv[3]
population = int(sys.argv[4])
extension = "pwo"
if len(sys.argv) == 6:
    extension = sys.argv[5]

# Read the dynamical matrix
nqirr = 0
while(os.path.exists(origin_dynmat + str(nqirr + 1))):
    nqirr += 1
if nqirr == 0:
    print("You must specify only the root of the dynamical matrix. The q star index will be added by this script.")
    raise IOError("Error while loading the dynamical matrix %s. Not found" % (origin_dynmat+ "1"))
dyn = CC.Phonons.Phonons(origin_dynmat, nqirr)

# Get the supercell
supercell = dyn.GetSupercell()

# Prepare the ensemble
ens = sscha.Ensemble.Ensemble(dyn, 0, supercell)

# Now get the output files
out_dir = os.path.dirname(PREFIX)
output_files = [out_dir + "/" + x for x in os.listdir(out_dir) if (PREFIX in out_dir + "/" + x) and ("."+extension in x)]
if len(output_files) == 0:
    print (INFO)
    raise IOError("Error, no good output file found. Please check your PREFIX argument")


N_configs = 0
energies = []
forces = []
stresses = []
xats = []
has_stress = True
for i, outfile in enumerate(output_files):
    try:
        ase_atm = ase.io.read(outfile)
    except: 
        continue
    struct = CC.Structure.Structure()
    struct.generate_from_ase_atoms(ase_atm) 

    # Read the energy and the force
    try:
        energy = ase_atm.get_total_energy() / Rydberg
        force = ase_atm.get_forces() / Rydberg
    except:
        continue

    # Try to read the stress
    if has_stress:
        try:
            stress = -ase_atm.get_stress(False)* Bohr**3 / Rydberg
        except:
            has_stress = False


    ens.structures.append(struct)
    forces.append(force)
    energies.append(energy)

    # Get the atom position
    xats.append(struct.coords.copy())

    if has_stress:
        stresses.append(stress) 
    N_configs += 1


ens.N = N_configs 
ens.energies = np.array(energies, dtype = np.float64)
ens.forces = np.array(forces, dtype = np.float64)
ens.xats = np.array(xats, dtype = np.float64)
ens.has_stress = has_stress
if has_stress:
    ens.stresses = np.array(stresses, dtype = np.float64)


# Save the ensemble
ens.save_bin(save_dir, population)