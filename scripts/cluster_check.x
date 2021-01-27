#!python

import sys, os
import numpy as np

import cellconstructor as CC
import cellconstructor.Structure
import cellconstructor.Phonons

import sscha
import sscha.Ensemble
import sscha.Cluster
import sscha.Calculator
import argparse

PROG_NAME = "check_cluster.x"
VERSION = "CHECK_CLUSTER, alpha version, not yet released."

info = """
CHECK CLUSTER SCRIPT
====================

This small utility will read an input file containing the cluster.
It will check for connectivity, and submit a small trial job and job array
to help you to configure your submission.

"""

# Define a custom function to raise an error if the file does not exist
def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exists!" % arg)
    else:   
        return arg


# Prepare the parser for the command line
parser = argparse.ArgumentParser(prog = PROG_NAME,
                                 description = info,
                                 version = VERSION, formatter_class = argparse.RawTextHelpFormatter)

# Add the arguments
parser.add_argument("-i", help="Path to the input file",
                    dest="filename", required = True, metavar = "FILE",
                    type = lambda x : is_valid_file(parser, x))
                    


# Get the arguments from the command line
args = parser.parse_args()

# Check if the mandatory arguments are defined
namelist = CC.Methods.read_namelist(args.filename)

# Setup the cluster
cluster = sscha.Cluster.Cluster()
cluster.setup_from_namelist(namelist)

# Check if there is a calculator
is_calc = False
if sscha.Calculator.__CALCULATOR_HEAD__ in namelist:
    calc = sscha.Calculator.prepare_calculator_from_namelist(namelist)
    is_calc = True

# Check the cluster connectivity
print "Checking communication..."
comm = cluster.CheckCommunication()
if not comm:
    print "Error, while trying to communicate with the cluster."
    print "Please, make sure your setup is correct."
    exit(1)
else:
    print "Communication established."

# If a calculator has been setup, perform a trial calculation
if is_calc:
    # Lets take a simple H system
    Hsys = CC.Structure.Structure(2)
    Hsys.atoms = ["H", "H"]
    Hsys.coords[0, :] = [0,0,0]
    Hsys.coords[1, :] = [1, 1, 1]
    Hsys.unit_cell[0, :] = [2,0,0]
    Hsys.unit_cell[1, :] = [0,2,0]
    Hsys.unit_cell[2, :] = [0,0,2]
    Hsys.has_unit_cell = True
    
    # Trial dyn
    dyn = CC.Phonons.Phonons(Hsys)
    
    # An ensemble of 10 equal elements
    ensemble = sscha.Ensemble.Ensemble(dyn, 0, (1,1,1))
    N_structs = 10
    ensemble.structures = [Hsys.copy() for x in xrange(N_structs)]
    ensemble.N = N_structs
    ensemble.energies = np.zeros(N_structs, dtype = np.float64)
    ensemble.forces = np.zeros( (N_structs, 2, 3), dtype = np.float64)
    
    # Submit the calculation as a single
    print "Submitting the calculation as a single..."
    cluster.compute_ensemble(ensemble, calc, False)
    
    print "Calculation ended."
    print "Printing energies (should be all equal):"
    print "\n".join(["%16.8f Ry" % x for x in ensemble.energies])
    print ""
    
    # Submit the parallel calculation
    print "Submitting the calculation as a batch..."
    print "(this is the way for big ensembles)"
    cluster.compute_ensemble_batch(ensemble, calc, False)
    
    print "Calculation ended."
    print "Printing energies (should be all equal):"
    print "\n".join(["%16.8f Ry" % x for x in ensemble.energies])
    print ""
    
print "Done."