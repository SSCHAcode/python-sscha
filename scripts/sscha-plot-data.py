import numpy as np
import matplotlib.pyplot as plt
import sys, os

import cellconstructor as CC, cellconstructor.Units


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
    
    plt.rcParams["font.family"] = "Liberation Serif"
    LBL_FS = 12
    DPI = 120
    TITLE_FS = 15

    freqs_files = [sys.argv[x] + '.freqs' for x in range(1, len(sys.argv))]
    minim_files = [sys.argv[x] + '.dat' for x in range(1, len(sys.argv))]

    # Check if all the files exist
    for f, m in zip(freqs_files, minim_files):
        assert os.path.exists(f), 'Error, file {} not found.'.format(f)
        assert os.path.exists(m), 'Error, file {} not found.'.format(m)

    # Load all the data
    freqs_data = np.concatenate([f for f in freqs_files])
    minim_data = np.concatenate([f for f in minim_files])
    
    fig_data, axarr = plt.subplots(nrows=2, ncols = 2, sharex = True, dpi = DPI)
    
    # Plot the data
    axarr[0,0].fill_between(minim_data[:,0], minim_data[:,1] - minim_data[:, 2]*.5 ,
                            minim_data[:, 1] + minim_data[:, 2] * .5, color = "aquamarine")
    axarr[0,0].plot(minim_data[:,0], minim_data[:,1], color = "k")
    axarr[0,0].set_ylabel("Free energy / unit cell [meV]", fontsize = LBL_FS)


    axarr[1,0].fill_between(minim_data[:,0], minim_data[:,3] - minim_data[:, 4]*.5 ,
                            minim_data[:, 3] + minim_data[:, 4] * .5, color = "aquamarine")
    axarr[1,0].plot(minim_data[:,0], minim_data[:,3], color = "k")
    axarr[1,0].set_ylabel("FC gradient", fontsize = LBL_FS)

    axarr[1,1].fill_between(minim_data[:,0], minim_data[:,5] - minim_data[:, 6]*.5 ,
                            minim_data[:, 5] + minim_data[:, 6] * .5, color = "aquamarine")
    axarr[1,1].plot(minim_data[:,0], minim_data[:,5], color = "k")
    axarr[1,1].set_ylabel("Structure gradient", fontsize = LBL_FS)
    axarr[1,1].set_xlabel("Good minimization steps", fontsize = LBL_FS)


    axarr[0,1].plot(minim_data[:,0], minim_data[:,7], color = "k")
    axarr[0,1].set_ylabel("Effective sample size", fontsize = LBL_FS)
    axarr[0,1].set_xlabel("Good minimization steps", fontsize = LBL_FS)
    fig_data.tight_layout()
    

    # Now plot the frequencies
    fig_freqs = plt.figure(dpi = DPI)
    ax = plt.gca()

    N_points, Nw = np.shape(freqs_data)

    for i in range(Nw):
        ax.plot(freqs_data[:, i] * CC.Units.RY_TO_CM)

    ax.set_xlabel("Good minimization steps", fontsize = LBL_FS)
    ax.set_ylabel("Frequency [cm-1]", fontsize = LBL_FS)
    ax.set_title("Frequcency evolution", fontsize = TITLE_FS)
    fig_freqs.tight_layout()

    plt.show()
    
if __name__ == "__main__":
    main()