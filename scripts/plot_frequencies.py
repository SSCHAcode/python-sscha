#!python
from __future__ import print_function

from numpy import *
from matplotlib.pyplot import *
import sys

RY_TO_CM = 109691.40235

print()
print()
print("PLOTTING THE FREQIENCIES")
print("========================")
print()
print()



if len(sys.argv) < 2:
    print("Specify the frequency files on command line.")
    print("If more than one, the frequencies will be concatenated in the order they are specified.")
    exit(1)


print("Loading the data...")

freq_data = concatenate([loadtxt(sys.argv[x]) for x in range(1, len(sys.argv))])

N_points, Nw = shape(freq_data)

print("Plotting...")

figure(dpi = 200)
for i in range(Nw):
    plot(freq_data[:, i] * RY_TO_CM)

xlabel("Step")
ylabel("Frequency [cm-1]")
title("Frequcency evolution")
tight_layout()


print("Done.")

show()
