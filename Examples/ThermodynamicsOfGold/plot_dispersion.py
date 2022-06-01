import cellconstructor as CC, cellconstructor.Phonons, cellconstructor.ForceTensor
import ase, ase.dft.kpoints

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors

import sys, os

NQIRR = 13
#CMAP = "Spectral_r"
PATH = "GXWXKGL"
N_POINTS = 1000

SPECIAL_POINTS = {"G": [0,0,0],
                  "X": [0, .5, .5],
                  "L": [.5, .5, .5],
                  "W": [.25, .75, .5],
                  "K": [3/8., 3/4., 3/8.]}

# Load the harmonic and sscha phonons
harmonic_dyn = CC.Phonons.Phonons('harmonic_dyn', NQIRR)
sscha_dyn = CC.Phonons.Phonons('sscha_T300_dyn', NQIRR)

# Get the band path
qpath, data = CC.Methods.get_bandpath(harmonic_dyn.structure.unit_cell,
                                      PATH,
                                      SPECIAL_POINTS,
                                      N_POINTS)
xaxis, xticks, xlabels = data # Info to plot correclty the x axis

# Get the phonon dispersion along the path
harmonic_dispersion = CC.ForceTensor.get_phonons_in_qpath(harmonic_dyn, qpath)
sscha_dispersion = CC.ForceTensor.get_phonons_in_qpath(sscha_dyn, qpath)

nmodes = harmonic_dyn.structure.N_atoms * 3

# Plot the two dispersions
plt.figure(dpi = 150)
ax = plt.gca()

for i in range(nmodes):
    lbl=None
    lblsscha = None
    if i == 0:
        lbl = 'Harmonic'
        lblsscha = 'SSCHA'
        
    ax.plot(xaxis, harmonic_dispersion[:,i], color = 'k', ls = 'dashed', label = lbl)
    ax.plot(xaxis, sscha_dispersion[:,i], color = 'r', label = lblsscha)

# Plot vertical lines for each high symmetry points
for x in xticks:
    ax.axvline(x, 0, 1, color = "k", lw = 0.4)
ax.axhline(0, 0, 1, color = 'k', ls = ':', lw = 0.4)

ax.legend()
    
# Set the x labels to the high symmetry points
ax.set_xticks(xticks)
ax.set_xticklabels(xlabels)

ax.set_xlabel("Q path")
ax.set_ylabel("Phonons [cm-1]")

plt.tight_layout()
plt.savefig("dispersion.png")
plt.show()

    

