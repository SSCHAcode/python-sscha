from __future__ import print_function

INFO = """
Get the symmetries in the polarization basis, and compare the two symmetrization procedures
followed by the C and Python code step by step.

Parameters:

 1) mode_a : int
 2) mode_b : int
 3) mode_c : int
 4) file : string (name of the file to be analized)

It will compare the symmetrization in the block of degeneracies identified by these modes.
"""

import numpy as np
import sys, os

import cellconstructor as CC
import cellconstructor.Phonons

import sscha, sscha.DynamicalLanczos
import sscha.Ensemble


def test_check_symmetries():
    if len(sys.argv) != 5:
        print(INFO)
        exit()

    mode_a = int(sys.argv[1])
    mode_b = int(sys.argv[2])
    mode_c = int(sys.argv[3])
    fname = sys.argv[4]

    if not os.path.exists(fname):
        raise IOError("Error, file '{}' not found.".format(fname))

    # Load the dynamical matrix and the ensemble
    dyn = CC.Phonons.Phonons("../dyn_gen_pop1_", 3)
    ens = sscha.Ensemble.Ensemble(dyn, 0 , dyn.GetSupercell())
    ens.load_bin("..", 1)

    lanczos = sscha.DynamicalLanczos.Lanczos(ens, mode = 2)
    lanczos.prepare_symmetrization()

    N_sym = lanczos.symmetries.shape[0]
    
    # read the file
    print("Reading the input file...")
    file_format = open(fname, "r")
    all_lines = file_format.readlines()
    file_format.close()

    print("Looking at the symmetries...")

    # Get the modes
    deg_space_a = []
    deg_space_b = []
    deg_space_c = []
    for line in all_lines:
        data = line.strip().split()

        # Break the cycle after readed the modes
        if  len(deg_space_a) * len(deg_space_b) * len(deg_space_c) == 0:
            if data[0] == "Mode":
                if int(data[1]) == mode_a:
                    list_deg = eval(line.strip()[line.strip().rfind("["):])
                    deg_space_a = list_deg
                if int(data[1]) == mode_b:
                    list_deg = eval(line.strip()[line.strip().rfind("["):])
                    deg_space_b = list_deg
                if int(data[1]) == mode_c:
                    list_deg = eval(line.strip()[line.strip().rfind("["):])
                    deg_space_c = list_deg
        else:
            break
                    

    deg_space_a.sort()
    deg_space_b.sort()
    deg_space_c.sort()

    print("Degenerate space:")
    print("a = ", deg_space_a)
    print("b = ", deg_space_b)
    print("c = ", deg_space_c)
    
    deg_dimension = len(deg_space_a) * len(deg_space_b) * len(deg_space_c)
    symmetries = np.zeros( (N_sym, deg_dimension, deg_dimension), dtype = np.double)
    new_a = 0
    new_b = 0
    new_c = 0
    for lline in all_lines:
        line = lline.strip().replace(",","")
        data = line.split()

        if data[0] == "Deg" and data[1] == "space":
            new_a = int(data[7])
            new_b = int(data[10])
            new_c = int(data[13])

   
            continue

        # Avoid reading further if we are not in the modes of interest
        if not new_a in deg_space_a:
            continue
        if not new_b in deg_space_b:
            continue
        if not new_c in deg_space_c:
            continue

        #print("Good looking! {} {} {}", new_a, new_b, new_c)

        if data[0] == "IN_VEC_OUT_DYN:":
            #print("Here!")
            if int(data[-1]) == 1:
                break # Avoid analyzing the transposed

            
            index_a = deg_space_a.index(new_a)
            index_b = deg_space_b.index(new_b)
            index_c = deg_space_c.index(new_c)

            sym_index = int(data[17])
            a, b, c = [int(x) for x in line[line.find("[")+1 : line.find("]")].split()]

            
            i1 = index_a * len(deg_space_b) * len(deg_space_c) + index_b * len(deg_space_c) + index_c
            
            index_a = deg_space_a.index(a)
            index_b = deg_space_b.index(b)
            index_c = deg_space_c.index(c)
            i2 = index_a * len(deg_space_b) * len(deg_space_c) + index_b * len(deg_space_c) + index_c


            #print("Adding element [{}, {}] = {}".format(i1, i2, float(data[3])))
            #print("{} {} {} | {} {} {}".format(new_a, new_b, new_c, a, b, c))
            symmetries[sym_index, i1, i2] = float(data[3])


    print("Creating the symmetries...")
    # Get the symmetries with the python way
    if not os.path.exists("data_sym"):
        os.makedirs("data_sym")
    
    for i in range(N_sym):
        sym_mat_py = np.einsum("ab,cd, ef ->acebdf", lanczos.symmetries[i, np.min(deg_space_a) : np.max(deg_space_a)+1, np.min(deg_space_a) : np.max(deg_space_a)+1],
                               lanczos.symmetries[i,  np.min(deg_space_b) : np.max(deg_space_b)+1,  np.min(deg_space_b) : np.max(deg_space_b)+1],
                               lanczos.symmetries[i,  np.min(deg_space_c) : np.max(deg_space_c)+1, np.min(deg_space_c) : np.max(deg_space_c)+1]).reshape((deg_dimension, deg_dimension))

        np.savetxt(os.path.join("data_sym", "py_symmetry_{:04d}.dat".format(i)), sym_mat_py)
        np.savetxt(os.path.join("data_sym", "c_symmetry_{:04d}.dat".format(i)), symmetries[i, :, :])


    print("Please, check the results!")

if __name__ == "__main__":
    test_check_symmetries()
