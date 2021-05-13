# -*- coding: utf-8 -*-
from __future__ import print_function
"""
This example uses the marconi cluster
to submit a force and energy calculation of a subset of 10 configurations
from the simple ensemble
"""

import ase
from ase.calculators.espresso import Espresso

import cellconstructor as CC
import cellconstructor.Phonons

import sscha
import sscha.Ensemble
import sscha.Cluster

# Load the dynamical matrix
dyn = CC.Phonons.Phonons("../ensemble_data_test/dyn")
ens = sscha.Ensemble.Ensemble(dyn, 0)

# Generate a random ensemble
ens.generate(10)

# Prepare the espresso calculator
# NOTE: these files should be located in the cluster $HOME/espresso/pseudo directory
pseudo = {"H": "H.pbe-rrkjus_psl.0.1.UPF",
                    "O": "O.pbe-n-rrkjus_psl.0.1.UPF"}

input_data = {"ecutwfc" : 45,
              "ecutrho" : 45*8,
              "conv_thr" : 1e-8,
              "occupations" : "fixed",
              "tstress" : True,
              "diskio" : "none",
              "tprnfor" : True}
KPTS = (3,3,2) # K points

calc = Espresso(input_data = input_data, pseudopotentials = pseudo, kpts = KPTS)

# Prepare the cluster
# marconi (setted from .ssh_config) => pippo@login.marconi.cineca.it
# No pwd, login with private key
cluster = sscha.Cluster.Cluster("marconi", partition_name="knl_usr_prod",
                                binary="$HOME/qe-6.2.1/bin/pw.x -npool $NPOOL -i PREFIX.pwi > PREFIX.pwo")

# Setup the working directory
cluster.workdir = "/marconi_work/IscrC_HydPD/tmp_calc"
cluster.n_nodes = 1
cluster.n_cpu = 32
cluster.account_name = "IscrB_COMRED"
cluster.n_pool = 4
cluster.load_modules = """module load autoload intel
module load autoload intelmpi
module load mkl
module load fftw"""


print("Sending the ensemble for the calculation.")
cluster.compute_ensemble(ens, calc)

ens.save_bin(".")
print("Ensemble saved.")
