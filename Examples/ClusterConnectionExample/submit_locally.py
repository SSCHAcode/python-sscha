# -*- coding: utf-8 -*-
from __future__ import print_function
"""
This example uses the marconi cluster
to submit a force and energy calculation of a subset of 10 configurations
from the simple ensemble
"""

import ase

import cellconstructor as CC
import cellconstructor.Phonons

from cellconstructor.calculators import Espresso

import sscha
import sscha.Ensemble
import sscha.Cluster

# Load the dynamical matrix
dyn = CC.Phonons.Phonons("dyn")
ens = sscha.Ensemble.Ensemble(dyn, 0)

# Generate a random ensemble
ens.generate(10)

# Prepare the espresso calculator
# NOTE: these files should be located in the cluster $HOME/espresso/pseudo directory
pseudo = {"H": "H.upf",
            "O": "O.upf"}

input_data = {
        "control" : {
            "disk_io" : "none",
              "tprnfor" : True,
              "tstress" : True,
            },
        "system" : {
            "ecutwfc" : 45,
            "input_dft" : "pbe",
              "ecutrho" : 45*8,
              "occupations" : "fixed"},
        "electrons" : {
              "conv_thr" : 1e-8
              }
        }
KPTS = (3,3,2) # K points

calc = Espresso(input_data = input_data, pseudopotentials = pseudo, kpts = KPTS)

# Prepare the cluster
# marconi (setted from .ssh_config) => pippo@login.marconi.cineca.it
# No pwd, login with private key
cluster = sscha.LocalCluster.LocalCluster("localhost", 
                                binary="pw.x -npool NPOOL -i PREFIX.pwi > PREFIX.pwo")

# Setup the working directory
cluster.workdir = "$SCRATCH/sscha_prova"
cluster.n_nodes = 1
cluster.n_cpu = 32
cluster.account_name = "IscrC_TDSTO"
cluster.n_pool = 4
cluster.partition_name = "g100_usr_prod"
cluster.load_modules = """
module load profile/chem-phys
module load autoload qe

export OMP_NUM_THREADS=1
"""

cluster.setup_workdir()

print("Sending the ensemble for the calculation.")
cluster.compute_ensemble(ens, calc)

ens.save_bin(".")
print("Ensemble saved.")
