#!python

from numpy import *
from matplotlib.pyplot import *
import sys

RY_TO_CM = 109691.40235

if len(sys.argv) < 2:
    print ("Specify the frequency files.")
    exit(1)

freq_data = concatenate([loadtxt(sys.argv[x]) for x in range(1, len(sys.argv))])

N_points, Nw = shape(freq_data)

figure(dpi = 200)
for i in range(Nw):
    plot(freq_data[:, i] * RY_TO_CM)

xlabel("Step")
ylabel("Frequency [cm-1]")
title("Frequcency evolution")
tight_layout()
show()
