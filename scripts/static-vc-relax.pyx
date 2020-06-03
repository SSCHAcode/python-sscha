#!python
from __future__ import print_function

# -*- coding: utf-8 -*-

import argparse

import sys, os
import numpy as np
import difflib

try:
    import ase
    #from ase.calculators.espresso import Espresso
    from ase.optimize import BFGS
    from ase.units import GPa
    __USE_ASE__ = True
except:
    __USE_ASE__ = False
    

import cellconstructor as CC
import cellconstructor.Phonons

import sscha
import sscha.Optimizer, sscha.Calculator
import sscha.SchaMinimizer


INFO = """
STATIC VC RELAX
===============

This script is used to perform a static variable cell relaxation.
You can use the fixed volume options.
In this way it is easier to compare the Ab-initio vs sscha results.

The calculations are runned locally, even if a cluster is specified.
"""


PROG_NAME = "static-vc-relax"
VERSION = "0.1"

if not __USE_ASE__:
    sys.stderr.write("Error, ASE library not found in PYTHONPATH.\n")
    raise ImportError("Error, the ASE library is required to run this script.")


# Define a custom function to raise an error if the file does not exist
def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exists!" % arg)
    else:   
        return arg
    

# Prepare the parser for the command line
parser = argparse.ArgumentParser(prog = PROG_NAME,
                                 description = INFO,
                                 formatter_class = argparse.RawTextHelpFormatter)


# Add the arguments
parser.add_argument("-i", help="Path to the input file",
                    dest="filename", required = True, metavar = "FILE",
                    type = lambda x : is_valid_file(parser, x))

# Get the arguments from the command line
args = parser.parse_args()
filename = args.filename

# Define the good namelist variables
__STATIC_REL_HEADER__ = "static_vc_relax"
__STATIC_REL_MAX_KA__ = "max_ka"
__STATIC_REL_START_KA__ = "start_ka"
__STATIC_REL_MAX_FORCE__ = "force_thr"
__STATIC_REL_BULK_MODULUS__ = "bulk_modulus"
__STATIC_REL_PRESS_THR__ = "press_thr"
__STATIC_REL_TARGET_P__ = "press"
__STATIC_REL_V_FIXED__ = "fix_volume"
__STATIC_REL_SCF__ = "structure"
__STATIC_REL_SYM_THR__ = "symm_thr"
__STATIC_REL_USE_SYM__ = "use_symmetries"
__STATIC_REL_FILE__ = "save_file"

__STATIC_REL_RWK__ = [__STATIC_REL_SCF__, __STATIC_REL_PRESS_THR__, __STATIC_REL_MAX_FORCE__]
__STATIC_REL_KEYS__ = [__STATIC_REL_MAX_KA__, __STATIC_REL_START_KA__,
                       __STATIC_REL_MAX_FORCE__, __STATIC_REL_BULK_MODULUS__,
                       __STATIC_REL_PRESS_THR__, __STATIC_REL_TARGET_P__,
                       __STATIC_REL_V_FIXED__, __STATIC_REL_SCF__,
                       __STATIC_REL_SYM_THR__, __STATIC_REL_USE_SYM__,
                       __STATIC_REL_FILE__]


# Load the namespace and check for all the variables
data_relax = CC.Methods.read_namelist(filename)
calc = sscha.Calculator.prepare_calculator_from_namelist(data_relax)

# Load the keywords
if not __STATIC_REL_HEADER__ in data_relax.keys():
    raise ValueError("Error, the input file must contain the %s namespace." % __STATIC_REL_HEADER__)

c_info = data_relax[__STATIC_REL_HEADER__]
keys = c_info.keys()

# Check for syntax errors
for k in keys: 
    if not k in __STATIC_REL_KEYS__:
        sys.stderr.write( "Error with the key: "+ k+"\n")
        s = "Did you mean something like:" + str( difflib.get_close_matches(k, __STATIC_REL_KEYS__))
        sys.stderr.write(s + "\n")
        raise IOError("Error in "+__STATIC_REL_HEADER__+" namespace: key '" + k +"' not recognized.\n" + s)

for k in __STATIC_REL_RWK__:
    if not k in keys:
        raise IOError("Error, the key %s is requested in namespace %s." % (k,__STATIC_REL_HEADER__))


# Parse all the keys
start_from = c_info[__STATIC_REL_SCF__]
press_thr = c_info[__STATIC_REL_PRESS_THR__]
max_force = c_info[__STATIC_REL_MAX_FORCE__]
max_ka = 50
start_ka = 1
fix_volume = False
static_bulk_modulus = 100
target_p = 0
symm_thr = np.float64(1e-5)
use_symmetries = True
append_file = None

if __STATIC_REL_MAX_KA__ in keys:
    max_ka = int(c_info[__STATIC_REL_MAX_KA__])
if __STATIC_REL_START_KA__ in keys:
    start_ka = int(c_info[__STATIC_REL_START_KA__])
if __STATIC_REL_BULK_MODULUS__ in keys:
    static_bulk_modulus = c_info[__STATIC_REL_BULK_MODULUS__]
if __STATIC_REL_V_FIXED__ in keys:
    fix_volume = c_info[__STATIC_REL_V_FIXED__]
if __STATIC_REL_TARGET_P__ in keys:
    target_p = c_info[__STATIC_REL_TARGET_P__]
    if fix_volume:
        raise ValueError("Error, you cannot fix the volume and choose a target pressure.")
if __STATIC_REL_USE_SYM__ in keys:
    use_symmetries = bool(c_info[__STATIC_REL_USE_SYM__])
if __STATIC_REL_SYM_THR__ in keys:
    if not use_symmetries:
        raise ValueError("Error, symmetries disabled, key '%s' has no meaning." % __STATIC_REL_SYM_THR__)
    symm_thr = np.float64(c_info[__STATIC_REL_SYM_THR__])
if __STATIC_REL_FILE__ in keys:
    append_file = str(c_info[__STATIC_REL_FILE__])

struct = CC.Structure.Structure()
struct.read_scf(start_from)

# Get the symmetries
qe_sym = CC.symmetries.QE_Symmetry(struct, symm_thr)
qe_sym.SetupQPoint()
symmetries = qe_sym.GetSymmetries()

ase_struct = struct.get_ase_atoms()

# Setup the cell optimizer
cell_SD = sscha.Optimizer.UC_OPTIMIZER(struct.unit_cell)
static_bulk_modulus /= sscha.SchaMinimizer.__evA3_to_GPa__ 
cell_SD.alpha = 1 / (3 * static_bulk_modulus * np.linalg.det(struct.unit_cell))
cell_SD.min_alpha_factor = 0.1
cell_SD.reset_strain = True

# Perform the optimization at fixed volume
for i in range(start_ka, max_ka):
    # Setup the ASE calculation
    calc.set_label("espresso%05d" % i)
    ase_struct.set_calculator(calc)

    # Perform the fixed cell relaxation
    print ("Relaxation %d:" % i)
    relaxation = BFGS(ase_struct, trajectory = "struct-relax-%04d.traj" % i)
    relaxation.run(fmax = max_force)

    # Get the stress
    stress = - ase_struct.get_stress(False) 
    
    print ("Final stress [GPa]:")
    print (np.array_str(stress / GPa, precision=4, suppress_small=True))
    print() 
    print()


    # Perform the cell optimization
    new_struct = CC.Structure.Structure()
    new_struct.generate_from_ase_atoms(ase_struct)
    
    # Get the info from the structure
    energy = ase_struct.get_total_energy() / new_struct.N_atoms
    press = np.trace(stress) / (3 * GPa)
    tot_force = np.sqrt(np.sum(ase_struct.get_forces()**2)) / new_struct.N_atoms
    
    if not fix_volume:
        stress -= np.eye(3, dtype = np.float64) * target_p * GPa
    
    # Apply the symmetries on the structure
    if use_symmetries:
        print ("Imposing the simmetries on the final structure and stress.")
        new_struct.impose_symmetries(symmetries)
        
        # Impose the symmetries on the stress
        qe_sym.ApplySymmetryToMatrix(stress)
    
    cell = new_struct.unit_cell.copy()
    cell_SD.UpdateCell(cell, stress, fix_volume = fix_volume, verbose = False)
    new_struct.change_unit_cell(cell)
    ase_struct = new_struct.get_ase_atoms()
    
    # Get the volume per atom
    volume = np.abs(np.linalg.det(cell)) / new_struct.N_atoms

    # Save as SCF the new structure
    new_struct.save_scf("relaxed-%04d.scf" % i)

    # Check if it converged
    if fix_volume:
        conv =  np.eye(3, dtype = np.float64) * np.trace(stress/GPa)/3
    else:
        conv = np.zeros((3,3), dtype = np.float64) 
        
    if np.sqrt(np.sum( (conv - stress/GPa)**2)) < press_thr:
        break

new_struct.save_scf("final_relaxed.scf")
print ("Final structure saved in final_relaxed.scf")


if not append_file is None:
    # Add the header if it does not exists
    if not os.path.exists(append_file):
        myfile = open(append_file, "w")
        myfile.write("# Pressure [GPa], Volume per atom [A^3], energy per atom [eV], average quadratic force [eV/A]\n")
        myfile.close()             
                     
    with open(append_file, "a") as myfile:
        myfile.write("%16.8f %16.8f %16.8f %16.8f\n" % (press, volume, energy, tot_force))
        myfile.close()

    print ("Info on the last minimization saved in '%s'" % append_file)


print ("Done.")
    
