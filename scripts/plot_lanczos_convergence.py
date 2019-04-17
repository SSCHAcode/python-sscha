#!python
from __future__ import print_function
from __future__ import division

import sys, os

import numpy as np
import matplotlib.pyplot as plt

import cellconstructor as CC
import cellconstructor.Phonons

import sscha, sscha.DynamicalLanczos
import sscha.Ensemble

# Check how many files are here
files = [x for x in  os.listdir(".") if "LANCZOS_STEP" in x and ".npz" in x]

data = []
nat = 0
print ("Reading the lanczos files (N = %d)..." % len(files))
for f in files:
    l_tmp = sscha.DynamicalLanczos.Lanczos()
    l_tmp.load_status(f)
    data.append(l_tmp)
    nat = l_tmp.nat

# Order the data
l_steps = []
for d in data:
    l_steps.append(len(d.a_coeffs))

sort_mask = np.argsort(l_steps)
sorted_data = [data[x] for x in sort_mask]
sorted_steps = [l_steps[x] for x in sort_mask]

# Now get the static converge
W_START = 0
W_END = 10000 / CC.Phonons.RY_TO_CM
NW = 10000
SMEARING = 10 / CC.Phonons.RY_TO_CM

print ("Computing the static responce...")
freqs = np.zeros((len(data), 3*nat))
dynamical = np.zeros((len(data), NW))
dynamical_lenmann = np.zeros((len(data), NW))
w_array = np.linspace(W_START, W_END, NW)

for i, dat in enumerate(sorted_data):
    fc_odd = dat.get_static_odd_fc(False)
    fc_odd /= np.sqrt(np.outer(dat.m, dat.m))
    w, p = np.linalg.eigh(fc_odd)
    sign_mask = np.sign(w)
    w = sign_mask * np.sqrt(np.abs(w))
    w *= CC.Phonons.RY_TO_CM
    freqs[i, :] = w 

    dynamical_lenmann[i, :] = dat.get_spectral_function_from_Lenmann(w_array, SMEARING, False)
    dynamical[i, :] = -np.imag( dat.get_green_function_continued_fraction(w_array, use_terminator = True, smearing = 0))
    
    if i % 10 == 0:
        sys.stderr.write("\rProgress %% %d" % (i * 100 / len(sorted_data)))
        sys.stderr.flush()
sys.stderr.write("\n")

print ("Plotting the results...")

plt.figure(dpi = 160)
plt.title("Static Odd")
for i in range(3*nat):
    plt.plot(sorted_steps, freqs[:, i])

plt.xlabel("Lanczos step")
plt.ylabel("frequencies [cm-1]")
plt.tight_layout()

plt.figure(dpi = 160)
plt.title("Green function from Lenmann")
plt.imshow(np.log(dynamical_lenmann), aspect = "auto", origin = "upper", extent = [W_START*CC.Phonons.RY_TO_CM, W_END*CC.Phonons.RY_TO_CM, max(l_steps), min(l_steps)])
plt.xlabel("Frequency [cm-1]")
plt.ylabel("Steps")
plt.colorbar()
plt.tight_layout()


plt.figure(dpi = 160)
plt.title("Green function from continued fraction")
plt.imshow(np.log(dynamical), aspect = "auto", origin = "upper", extent = [W_START*CC.Phonons.RY_TO_CM, W_END*CC.Phonons.RY_TO_CM, max(l_steps), min(l_steps)])
plt.xlabel("Frequency [cm-1]")
plt.ylabel("Steps")
plt.colorbar()
plt.tight_layout()

plt.show()