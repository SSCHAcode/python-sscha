from __future__ import print_function

import ase
from ase.calculators.espresso import Espresso

import cellconstructor as CC
import cellconstructor.Structure

import sys, os

pseudo = {"H":"H.pbe-rrkjus_psl.1.0.0.UPF",
         "S" : "S.pbe-nl-rrkjus_psl.1.0.0.UPF"}

input_data = {"ecutwfc" : 35,
              "ecutrho" : 350,
              "occupations" : "smearing",
              "input_dft" : "blyp",
              "mixing_beta" : 0.2,
              "conv_thr" : 1e-9,
              "degauss" : 0.02,
              "smearing" : "mp",
              "pseudo_dir" : "."}


calc = Espresso(pseudopotentials =pseudo,
                input_data = input_data,
                kspacing = 0.06)

struct = CC.Structure.Structure()
struct.read_scf("H3S.scf")

ase_struct = struct.get_ase_atoms()
ase_struct.set_calculator(calc)

print("Computing the total energy...")
total_energy = ase_struct.get_total_energy()


phonons = """
&inputph
   ldisp = .true.
   nq1 = 2
   nq2 = 2
   nq3 = 2
&end
"""

with open("ph.in", "w") as f:
    f.writelines(phonons)
    
print("Computing the harmonic dynamical matrix")
os.system("/usr/bin/mpirun -np 4 ph.x -npool 4 -i ph.in > ph.out")
