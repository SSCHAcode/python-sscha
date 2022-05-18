import cellconstructor as CC, cellconstructor.Phonons
import numpy as np
import matplotlib.pyplot as plt

import scipy, scipy.optimize

import sys, os

"""
This simple scripts plot the thermal expansion.
It either directly plots the file produced at the end of thermal_expansion.py
or it processes the dynamical matrices generated throughout the minimization.

It also fits the V(T) cuve and estimates the volumetric thermal expansion coefficient

alpha = dV/dT / V

At 300 K
"""


# Load all the dynamical matrices and compute volume
DIRECTORY = "thermal_expansion"
FILE = os.path.join(DIRECTORY, "thermal_expansion.dat")


# Check if the file with the thermal expansion data exists
if not os.path.exists( FILE):
    # If the simulation is not ended, load the volume from the dynamical matrices
    
    # Get the dynamical matrices
    all_dyn_files = [x for x in os.listdir(DIRECTORY) if "sscha" in x and x.endswith("dyn1")]
    temperatures = [float(x.split("_")[-2].replace("T", "")) for x in all_dyn_files]

    # Now sort in order of temperature
    sortmask = np.argsort(temperatures)
    all_dyn_files = [all_dyn_files[x] for x in sortmask]
    temperatures = np.sort(temperatures)

    volumes = np.zeros_like(temperatures)

    for i, fname in enumerate(all_dyn_files):
        # Load the dynamical matrix
        # The full_name means that we specify the name including the tail 1
        dyn = CC.Phonons.Phonons(os.path.join(DIRECTORY, fname), full_name = True)
        volumes[i] = dyn.get_volumes()

else:
    # Load the data from the final data file
    temperatures, volumes = np.loadtxt(FILE, unpack = True)


# Prepare the figure and plot the V(T) from the sscha data
plt.figure(dpi = 150)
plt.scatter(temperatures, volumes, label = "SSCHA data")

# Fit the data with a quadratic curve
def parabola(x, a, b, c):
    return a + b*x + c*x**2
def diff_parab(x, a, b, c):
    return b + 2*c*x

popt, pcov = scipy.optimize.curve_fit(parabola, temperatures, volumes,
                                      p0 = [0,0,0])

# Evaluate the volume thermal expansion
vol_thermal_expansion = diff_parab(300, *popt) / parabola(300, *popt)
print("Vol thermal expansion: {} x 10^6  K^-1".format(vol_thermal_expansion * 1e6))
plt.text(0.6, 0.2, r"$\alpha_v = "+"{:.1f}".format(vol_thermal_expansion*1e6)+r"\times 10^6 $ K$^{-1}$",
         transform = plt.gca().transAxes)


# Plot the fit
_t_ = np.linspace(np.min(temperatures), np.max(temperatures), 1000)
plt.plot(_t_, parabola(_t_, *popt), label = "Fit")

# Adjust the plot adding labels, legend, and saving in eps
plt.xlabel("Temperature [K]")
plt.ylabel(r"Volume [$\AA^3$]")
plt.legend()
plt.tight_layout()
plt.savefig("thermal_expansion.eps")
plt.show()
