#!python
# -*- coding: utf-8 -*-

from __future__ import print_function
from __future__ import division

"""
This is part of the program python-sscha
Copyright (C) 2018  Lorenzo Monacelli

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>. 
"""
import os
import time

import cellconstructor as CC
import cellconstructor.Methods

import sscha, sscha.SchaMinimizer
import sscha.Ensemble
import sscha.Relax
import sscha.Utilities
import argparse

PROG_NAME = "sscha"
VERSION = "1.1"
DESCRIPTION = r"""




* * * * * * * * * * * * * * * * * * * * * * * * * * *
*                                                   *
* STOCHASTIC SELF-CONSISTENT HARMONIC APPROXIMATION *
*                                                   *
* * * * * * * * * * * * * * * * * * * * * * * * * * *



Wellcome to the SSCHA code. 
This stand-alone script performs SSCHA minimization 
according to the options parsed from the input file.
Note, you can also run the SSCHA with the python-APIs. 
The python-APIs give access to much more customizable options,
not available with the stand-alone app.

Refer to the documentation for more details.


"""


# Define a custom function to raise an error if the file does not exist
def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exists!" % arg)
    else:   
        return arg

# Prepare the parser for the command line
parser = argparse.ArgumentParser(prog = PROG_NAME,
                                 description = DESCRIPTION,
                                 formatter_class = argparse.RawTextHelpFormatter)

# Add the arguments
parser.add_argument("-i", help="Path to the input file",
                    dest="filename", required = True, metavar = "FILE",
                    type = lambda x : is_valid_file(parser, x))
parser.add_argument("--plot-results", 
                    help="Plot the results",
                    dest="plot", action = "store_true")

parser.add_argument("--save-data", dest="save_data", metavar = "FILE",
                    help = "Save the minimization details")


# Get the arguments from the command line
args = parser.parse_args()

# Get the input filename and initialize the sscha minimizer.
minim = sscha.SchaMinimizer.SSCHA_Minimizer()

# Open the input file
fp = open(args.filename, "r")
line_list = fp.readlines()
fp.close()

# Check if the mandatory arguments are defined
namelist = CC.Methods.read_namelist(line_list)
if not sscha.SchaMinimizer.__SCHA_NAMELIST__ in namelist.keys():
    raise IOError("Error, %s namelist required." % sscha.SchaMinimizer.__SCHA_NAMELIST__)

# Check if the relax is hypothized
if sscha.Relax.__RELAX_NAMESPACE__ in namelist.keys():
    relax = sscha.Relax.SSCHA()
    print ("Found a %s namespace!" % sscha.Relax.__RELAX_NAMESPACE__)
    print ("Performing a whole relaxation.")
    
    relax.setup_from_namelist(namelist)
    
    # Save the final dynamical matrix
    print ("Saving the final dynamical matrix in %s" % (relax.minim.dyn_path + "_popuation%d_" % relax.minim.population))
    relax.minim.dyn.save_qe(relax.minim.dyn_path + "_population%d_" % minim.population)
    
    # Plot the results
    relax.minim.plot_results(args.save_data, args.plot)
    
    exit(0)


# Parse the input file
t1 = time.time()
print ("Loading everything...")
minim.setup_from_namelist(namelist)
t2 = time.time()
print ("Loaded in %.4f seconds." % (t2 - t1))

# Start the minimization
print ()
print ("Initialization...")
minim.init()
print ("System initialized.")
minim.print_info()
print ()

# If present, setup the custom functions
if sscha.Utilities.__UTILS_NAMESPACE__ in namelist:
    cfs = sscha.Utilities.get_custom_functions_from_namelist(namelist, minim.dyn)    
    minim.run(custom_function_pre = cfs[0], 
              custom_function_gradient = cfs[1],
              custom_function_post = cfs[2])
else:
    minim.run()
    
print ()
print ("Minimization ended.")
minim.finalize(2)

# Save the final dynamical matrix
print ("Saving the final dynamical matrix in %s" % ('final_sscha_dyn_pop%d' % minim.population))
print ('Saving the final centroids into {}'.format('final_sscha_structure_pop%d.scf' % minim.population))
print ()
minim.dyn.save_qe('final_sscha_dyn_pop%d' % minim.population)
minim.dyn.structure.save_scf('final_sscha_structure_pop%d.scf' % minim.population)

# Plot the results
minim.plot_results(args.save_data, args.plot)
