#!python


import sys, os

import cellconstructor as CC
import cellconstructor.Phonons

import sscha, sscha.DynamicalLanczos
import sscha.Ensemble

# Create a fake dyn
fake_struct = CC.Structure.Structure(5)
dyn = CC.Phonons.Phonons(fake_struct)
ens = sscha.Ensemble.Ensemble(dyn, 0)

# Check how many files are here
files = [x for x in  os.listdir(".") if "LANCZOS_STEP" in x and ".npz" in x]

data = []
for f in files:
    l_tmp = sscha.DynamicalLanczos.Lanczos(ens)
    l_tmp.load_status(f)
    data.append(l_tmp)

# TODO: TO be completed