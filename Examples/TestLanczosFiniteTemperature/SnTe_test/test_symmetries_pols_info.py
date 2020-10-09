import sscha, sscha.DynamicalLanczos, sscha.Ensemble
import cellconstructor as CC, cellconstructor.Phonons

import numpy as np
import spglib

dyn = CC.Phonons.Phonons("SnTe_sscha", 3)
dyn.Symmetrize(use_spglib = True)

w, pols = dyn.DiagonalizeSupercell()

ss = dyn.structure.generate_supercell(dyn.GetSupercell())
syms = CC.symmetries.GetSymmetriesFromSPGLIB(spglib.get_symmetry(ss.get_ase_atoms()))                                             
syms_pols_cc = CC.symmetries.GetSymmetriesOnModes(syms, ss, pols)

# Test the symmetries here
npols = len(w)
n_syms = len(syms)

tocm = CC.Units.RY_TO_CM

print("Test the symmetries...")
for i in range(n_syms):
    stop = False
    for x in range(npols):
        for y in range(npols):
            
            # Check weak condition
            alpha_x = syms_pols_cc[i, :, x] * w 
            alpha_y = syms_pols_cc[i, :, y] * w 

            Talpha_x = syms_pols_cc[i, x, :] * w
            Talpha_y = syms_pols_cc[i, y, :] * w

            # The scalar product between them must be 0
            scalar_dot = alpha_x.dot(alpha_y)
            Tscalar_dot = Talpha_x.dot(Talpha_y)

            #print("Sym {} | Modes {} | {} => scalar 1 = {} | scalar 2 = {}".format(i, x, y, scalar_dot, Tscalar_dot))

            if  x != y:
                assert np.abs(scalar_dot) < 1e-4, "Scalar DOT: {}".format(scalar_dot)
                
            if np.abs(w[x] - w[y]) > 1e-10:
                if syms_pols_cc[i, x, y] > 1e-10:
                    print("Error on symmetry {}: w[{}] = {} | w[{}] = {} => S = {}".format(i, x, w[x] * tocm,
                                                                                           y, w[y] * tocm,
                                                                                           syms_pols_cc[i, x, y]))
                    stop = True
                    

    assert not stop, "Error, the symmetry {} violated the degeneracies.".format(i)
                

ens = sscha.Ensemble.Ensemble(dyn, 250, dyn.GetSupercell())
ens.load_bin("ensemble", 1)

lanczos = sscha.DynamicalLanczos.Lanczos(ens, mode = 2)
lanczos.ignore_v4 = True
lanczos.ignore_harmonic = True
lanczos.ignore_v3 = False
lanczos.prepare_symmetrization(no_sym = False)

Ns, Np, __ = lanczos.symmetries.shape
print("SYM:", Ns)

not_mask = np.ones((Np, Np), dtype = bool)
for j in range(Np):
    ndim = lanczos.N_degeneracy[j]


    for kx in range( ndim):
        x = lanczos.degenerate_space[j][kx]
        for ky in range(ndim):
            y = lanczos.degenerate_space[j][ky]
            not_mask[x, y] = False



# Print the symmetrization info
for i in range(lanczos.n_modes):
    print("Mode {} | freq = {:10.3} cm-1 | Degeneracy = ".format(i, lanczos.w[i] * CC.Units.RY_TO_CM), lanczos.degenerate_space[i])

    

for i in range(Ns):

    new_sym = lanczos.symmetries[i, not_mask]
    thr = np.max(np.abs(new_sym))
    print("Sym {} satisfy degeneracies by {}".format(i, thr))

    if thr > 1e-12:
        # Check which is the element that gives problem
        for x in range(Np):
            for y in range(Np):
                if not not_mask[x,y]:
                    continue

                if np.abs(lanczos.symmetries[i, x, y]) > 1e-12:
                    print("The two crossing are ({}, {}) with symmetry = {}".format(x,y, lanczos.symmetries[i, x, y]))
                    print("The frequencies are: {} => {} cm-1 | {} => {} cm-1".format(x, lanczos.w[x] * CC.Units.RY_TO_CM,
                                                                                      y, lanczos.w[y] * CC.Units.RY_TO_CM))
    assert thr < 1e-12, "Degenerate space violated by {}".format(thr)
    
