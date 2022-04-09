import numpy as np
import matplotlib.pyplot as plt
import sys, os


DESCRIPTION = '''
This program plots the results of a SSCHA minimization 
saved with the Utilities.IOInfo class.

It plots both frequencies and the free energy minimization.

Many files can be specified, in this case, they are plotted one after the other,
separated by vertical dashed lines.

Usage:
 
 >>> sscha-plot-data.py  file1  ...

The code looks for data file called file.dat and file.freqs, 
containing, respectively, the minimization data and the auxiliary frequencies.

These files are generated only by python-sscha >= 1.2.

'''


def main():
    print(DESCRIPTION)

    if len(sys.argv) < 2:
        ERROR_MSG = '''
Error, you need to specify at least one file.
Exit failure!
'''
        print(ERROR_MSG)
        raise ValueError(ERROR_MSG)
    
    freqs_files = [sys.argv[x] + '.freqs' for x in range(1, len(sys.argv))]
    minim_files = [sys.argv[x] + '.dat' for x in range(1, len(sys.argv))]

    # Check if all the files exist
    for f, m in zip(freqs_files, minim_files):
        assert os.path.exists(f), 'Error, file {} not found.'.format(f)
        assert os.path.exists(m), 'Error, file {} not found.'.format(m)

    # Load all the data
    freqs_data = np.concatenate([f for f in freqs_files])
    minim_data = np.concatenate([f for f in minim_files])
    
    
    
    
