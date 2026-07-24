import cellconstructor as CC
import cellconstructor.Phonons
import cellconstructor.Structure
import cellconstructor.symmetries
import cellconstructor.Methods

import SCHAModules

#import FModules

import spglib
import numpy as np
import sys
import time

"""
This module contains the functions employed for the symmetry classification of atoms and
wave-vector orbits.
"""

def map_singlet(dyn, symprec=1e-5, verbose=False):
    """
    Classifies atomic singlets, and returns the atomic map for each symmetry operation.

    Parameters
    ----------
        - dyn: object
            Cellconstructor dynamical matrix.
        - symprec: float
            Tolerance parameter for spglib in symmetry detection.
            Defaults to 1e-5.
        - verbose: bool
            If True prints information during execution.
            Defaults to False.
    
    Returns
    -------
        - mapping: np.ndarray
            Atomic mapping for each symmetry operation. Dimension [Natom_sc,Nsym].
        - orbit1s: np.ndarray
            Symmetry mapping for each atomic mapping. Dimension [Natom_sc,Nsym].
        - rot_cart: np.ndarray
            Variable containing all the symmetry operations of the crystal. Dimension [Nsym,3,3]
        - map_uc: np.ndarray
            An array that works as a tool to know which index i we need to consider, taking
            index j as reference (unit-cell). Dimension [Natom_sc, Natom_sc].
        - map_tr: np.ndarray
            This array says which is the translation employed to map atom index i using index
            j as reference (unit-cell). Dimension [Natom_sc, Natom_sc].
        - T_list: np.ndarray
            Contains all the translations in unit-cell fractional units. Dimension [Nq,3].
        - T_list_frac: np.ndarray
            Contains all the translations in supercell fractional units. Dimension [Nq,3].
    """

    start_time = time.time()
    sg = spglib.get_spacegroup(dyn.structure.get_spglib_cell(), symprec)

    if verbose:
        print("Initial SG=", sg)
        print("===== STARTING SINGLET CLASSIFICATION =====")

    spg_syms = spglib.get_symmetry(dyn.structure.get_spglib_cell(), symprec)
    sym_uc = cellconstructor.symmetries.GetSymmetriesFromSPGLIB(spg_syms, regolarize=False)

    Nsym = len(sym_uc)

    # Obtain number of supercells. GetSupercell returns the modulation as 
    # 3 dimensional list, sc_size. 
    sc_size = dyn.GetSupercell() 
    Nsupercell = np.prod(sc_size)

    # Obtain point group symmetries in cryst coord respect sc: 
    sym_list = np.zeros((Nsym,3,4),dtype=np.float64)
    for isym in range(Nsym):
        sym_list[isym,:,:] = sym_uc[isym]
        for ll in range(3): # Transl respect sc
            sym_list[isym,ll,3] = sym_list[isym,ll,3]/sc_size[ll]

    # Create an object from the Phonons class of the SC
    dyn_sc = dyn.GenerateSupercellDyn(sc_size)

    # Get the symmetries of the supercell
    spg_syms_sc = spglib.get_symmetry(dyn_sc.structure.get_spglib_cell(), symprec)
    sym_list_sc = cellconstructor.symmetries.GetSymmetriesFromSPGLIB(spg_syms_sc, regolarize= False)
    
    # Obtain symmetries that are pure translations:
    translation = np.zeros((Nsupercell,3,4),dtype=np.float64)
    for i in range(Nsupercell):
        translation[i,:,:] = sym_list_sc[i*Nsym] # The first is always a translation
    
    T_list = np.empty([sc_size[0]*sc_size[1]*sc_size[2], 3], dtype=np.float64)
    T_list_frac = np.empty([sc_size[0]*sc_size[1]*sc_size[2], 3], dtype=np.float64)
    for Tx in range(sc_size[0]):
        for Ty in range(sc_size[1]):
            for Tz in range(sc_size[2]):
                index = Tx*sc_size[2]*sc_size[1]+Ty*sc_size[2]+Tz
                T_list[index] = np.array([Tx, Ty, Tz])
                T_list_frac[index] = np.array([Tx/sc_size[0],Ty/sc_size[1],Tz/sc_size[2]])
    
    map_uc, map_tr = _map_unitcell(dyn_sc, T_list_frac)

    #Obtain the rotations in cartesian coord:
    rot_cart = np.zeros((Nsym,3,3),dtype=np.float64)
    for isym in range(Nsym):
        rot_cart[isym] = cellconstructor.Methods.convert_matrix_cart_cryst2(sym_list[isym,:,:3], dyn.structure.unit_cell, cryst_to_cart = True)
        # Set elements smaller than the threshold to 0.0
        rot_cart[isym][np.abs(rot_cart[isym]) < 1e-12] = 0.0
    
    Natoms_sc = dyn_sc.structure.N_atoms
    
    singlet=np.zeros(1,dtype=np.intc)
    singlet_sym=np.zeros(1,dtype=np.intc)
    
    mapping = np.zeros((Natoms_sc, Nsym), dtype='<i4')
    orbit1s = np.zeros((Natoms_sc, Nsym), dtype='<i4')
    Nref_singlet=0
    for i in range(Natoms_sc):
        singlet = i
        Nref_singlet+=1
        Nequiv=0
        for isym in range(Nsym):
            sym_struct = dyn_sc.structure.copy()
            sym_struct.apply_symmetry(sym_list[isym],delete_original= True)
            # For each symmetry operation, we find which is equivalent to the i-th atom.
            irt = np.array(sym_struct.get_equivalent_atoms(dyn_sc.structure), dtype =np.intc)
            singlet_sym=irt[singlet]
            mapping[Nref_singlet-1, Nequiv] = singlet_sym
            orbit1s[Nref_singlet-1, Nequiv] = isym
            Nequiv+=1

    end_time = time.time()
    execution_time = end_time - start_time
    if verbose:
        print("execution_time in singlet recognition:", execution_time, " s")
        print("===============")
        print(" ")

    sys.stdout.flush()

    return(mapping, rot_cart, map_uc, map_tr, T_list, T_list_frac)


def recognize_doublet(dyn, orbit1a, map_uc, symprec=1e-5, verbose=False):
    """
    Classifies atomic doublets, and returns the atomic map for each symmetry operation,
    as well as the required information to reconstruct the 2nd order FCs from the symmetry independent
    doublets.

    Parameters
    ----------
        - dyn: object
            Cellconstructor dynamical matrix.
        - orbit1a: np.ndarray
            Atomic mapping for each symmetry operation. Dimension [Natom_sc,Nsym].
        - map_uc: An array that works as a tool to know which index i we need to consider, taking
            index j as reference (unit-cell). Dimension [Natom_sc, Natom_sc].
        - symprec: float
            Tolerance parameter for spglib in symmetry detection.
            Defaults to 1e-5.
        - verbose: bool
            If True prints information during execution.
            Defaults to False.

    Returns
    -------
        - orbit2a: np.ndarray
            Doublet mapping for each symmetry operation. Dimension [Nref2,tbd,2].
            In the second axis, valid inforation is up to the norbit[ref2]-th element.
        - orbit2s: np.ndarray
            Symmetry mapping for each doublet mapping. Dimension [Nref2,tbd,2].
            In the second axis, valid inforation is up to the norbit2[ref2]-th element.
        - norbit2: np.ndarray
            Number of equivalent doublets for each of the orbits. Dimension [Nref2].
        - indep_fc: np.ndarray
            Which are the independent cartesian elements from the 3x3 FC matrix. Dimension [Nref2,9].
            In the second axis, valid information is up to the n_indep_fc[ref2]-th element.
        - n_indep_fc: np.ndarray
            Number of independent cartesian elements for each reference doublet. Dimension [Nref2].
        - tensor: np.ndarray
            Tensor from which we will reconstruct the 3x3 FC matrix, knowing just the independet
            elements. In the case of triplet and quadruplets, we will return the kernel itself.
            Dimension [Nref2, tbd, 9, tbd]
    """

    start_time = time.time()

    if verbose:
        print("===== STARTING DOUBLET CLASSIFICATION =====")

    spg_syms = spglib.get_symmetry(dyn.structure.get_spglib_cell(), symprec)
    sym_uc = cellconstructor.symmetries.GetSymmetriesFromSPGLIB(spg_syms, regolarize=False)

    Nsym = len(sym_uc)

    # Obtain number of supercells. GetSupercell returns the modulation as 
    # 3 dimensional list, sc_size. 
    sc_size = dyn.GetSupercell() 
    Nsupercell = np.prod(sc_size)

    # Obtain point group symmetries in cryst coord respect sc: 
    sym_list = np.zeros((Nsym,3,4),dtype=np.float64)
    for isym in range(Nsym):
        sym_list[isym,:,:] = sym_uc[isym]
        for ll in range(3): # Transl respect sc
            sym_list[isym,ll,3] = sym_list[isym,ll,3]/sc_size[ll]
    
    # Create an object from the Phonons class of the SC
    dyn_sc = dyn.GenerateSupercellDyn(sc_size)
    
    # Get the symmetries of the supercell
    spg_syms_sc = spglib.get_symmetry(dyn_sc.structure.get_spglib_cell(), symprec)
    sym_list_sc = cellconstructor.symmetries.GetSymmetriesFromSPGLIB(spg_syms_sc, regolarize= False)
    
    # Obtain symmetries that are pure translations:
    translation = np.zeros((Nsupercell,3,4),dtype=np.float64)
    for i in range(Nsupercell):
        translation[i,:,:] = sym_list_sc[i*Nsym] # Lehena beti da translazioa
    
    
    #Obtain the rotations in cartesian coord:
    rot_cart = np.zeros((Nsym,3,3),dtype=np.float64)
    for isym in range(Nsym):
        rot_cart[isym] = cellconstructor.Methods.convert_matrix_cart_cryst2(sym_list[isym,:,:3], dyn.structure.unit_cell, cryst_to_cart = True)
        # Set elements smaller than the threshold to 0.0
        rot_cart[isym][np.abs(rot_cart[isym]) < 1e-12] = 0.0
    Natoms = dyn.structure.N_atoms
    Natoms_sc = dyn_sc.structure.N_atoms

    permutations=np.array([[0,1], [1,0]], dtype=np.intc)

    doublet=np.zeros(2,dtype=np.intc)
    doublet_perm=np.zeros(2,dtype=np.intc)
    doublet_sym=np.zeros(2,dtype=np.intc)
    
    # Initial size of array
    array_size = 16
    
    #Reference doublets:
    ref2 = np.zeros((array_size,2),dtype=np.intc)
    nref2 = 0
    # List with the number of equivalent doublets in each orbit
    norbit=np.zeros(array_size,dtype=np.intc)

    orbit2a = np.zeros((array_size, 2*Nsym, 2), dtype='<i4')
    orbit2s = np.zeros((array_size, 2*Nsym, 2), dtype='<i4')
    equilist = np.zeros((2*Nsym,2),dtype=np.intc)

    #List of all doublets:
    tot2=Natoms*Natoms_sc
    all2= np.zeros((tot2,2),dtype=np.intc)
    nall2=0

    #Create Rotation tensor: R_{αβ}^{α'β'}= R_{α}^{α'}*R_{β}^{β'}
    Rot = np.zeros((2,Nsym,9,9),dtype=np.float64) # Following Eq. (45), 2nd index refers to a single index combining α'β'
    #In order to take into account perm of indexes:
    cart_index = np.zeros(2,dtype=np.intc)
    
    t0 = time.time()
    for iperm in range(2):
        for isym in range(Nsym):
            for alphaprime in range(3):
                for betaprime in range(3):
                    indexprime = 3*alphaprime+betaprime
                    for alpha in range(3):
                        cart_index[0] = alpha
                        for beta in range(3):
                            cart_index[1] = beta
                            index = 3*alpha+beta
                            alphaperm=cart_index[permutations[iperm,0]]
                            betaperm =cart_index[permutations[iperm,1]]
                            Rot[iperm,isym,indexprime,index] = rot_cart[isym,alphaprime,alphaperm]*rot_cart[isym,betaprime,betaperm]

    R = np.zeros((array_size,2*Nsym,9,9),dtype=np.float64) #Rot matrix that goes from ref. doublet in an orbit to the rest
    
    #Gauss-Jordan elimination:
    constrain=np.zeros((2*Nsym*9,9),dtype=np.float64)
    #Overall results:
    kernel= np.zeros((array_size,9,9),dtype=np.float64)
    # List with the number of the independent FC elements for each ref. doub. 
    n_indep_fc = np.zeros((array_size),dtype=np.intc)
    # List with specific independent FC elements (possible values from 0 to 8) for each ref. doub.
    indep_fc = np.zeros((array_size,9),dtype=np.intc)
    
    #Matrix in which we apply Gauus-Jordan: M = R-I :
    M=Rot.copy()
    nontrivial=np.zeros((2,Nsym,9),dtype=np.intc)
    
    for iperm in range(2):
        for isym in range(Nsym):
            for indexprime in range(9):
                M[iperm,isym,indexprime,indexprime]-=1.0 #Possible values in the diagonal of M: 0,-1,-2
                for index in range(9):
                    if abs(M[iperm,isym,indexprime,index])>1e-12: # As a trivial constrain is a row of zeroes in R2
                        nontrivial[iperm,isym,indexprime]=1
                    else:
                        M[iperm,isym,indexprime,index]=0.0
    
    for ii in range(Natoms): # Unit cell, rest of the atoms (doublets) can be achieved by translations
        for jj in range(Natoms_sc):
            doublet[0]=ii
            doublet[1]=jj
            if _doublet_in_list(doublet, all2, nall2):
                continue
            for ll in range(2):
                ref2[nref2, ll] = doublet[ll]
            nref2+=1
            if (nref2 == array_size):
                ref2 = np.concatenate((ref2,ref2), axis = 0)
                norbit = np.concatenate((norbit, norbit), axis = 0)
                orbit2a = np.concatenate((orbit2a,orbit2a), axis = 0)
                orbit2s = np.concatenate((orbit2s,orbit2s), axis = 0)
                R = np.concatenate((R,R), axis = 0)
                kernel = np.concatenate((kernel,kernel), axis = 0)
                n_indep_fc = np.concatenate((n_indep_fc,n_indep_fc), axis = 0)
                indep_fc = np.concatenate((indep_fc,indep_fc), axis = 0)
                array_size<<=1 # Double the size
            Nequiv=0 # Number of equiv. doublets in the orbit of this ref. doublet
            equilist[:,:] = 0.0
            nconstrain = 0 # Amount of contrains for this ref. doublet
            constrain[:,:] = 0.0 #Restart
            for iperm in range(2):
                doublet_perm[0]=doublet[permutations[iperm,0]]
                doublet_perm[1]=doublet[permutations[iperm,1]]
                for isym in range(Nsym):
                    doublet_sym = [orbit1a[doublet_perm[0],isym], orbit1a[doublet_perm[1],isym]]
                    if (doublet_sym[0] >= Natoms): #First atom has to be always in the first unit cell
                         doublet_sym[1] = map_uc[doublet_sym[0],doublet_sym[1]]
                         doublet_sym[0] = map_uc[doublet_sym[0],doublet_sym[0]]
                    if Nequiv==0 or not(_doublets_are_equal(doublet,doublet_sym) or _doublet_in_list(doublet_sym, equilist, Nequiv)):
                        equilist[Nequiv] = doublet_sym
                        orbit2a[nref2-1, Nequiv] = doublet_sym
                        all2[nall2] = doublet_sym
                        orbit2s[nref2-1, Nequiv] = [iperm,isym]
                        R[nref2-1,Nequiv] = Rot[iperm,isym]
                        Nequiv+=1
                        nall2+=1                           
                    if _doublets_are_equal(doublet,doublet_sym):
                        for indexprime in range(9):
                            if nontrivial[iperm,isym,indexprime]: #1=true
                                for ll in range(9):
                                    constrain[nconstrain,ll] = M[iperm,isym,indexprime,ll]
                                nconstrain+=1
    
            norbit[nref2-1] = Nequiv
            # Gauss-Jordan Elimination for a doublet:
            constrain_reduced = constrain[:max(nconstrain,9),:9]
            kern, indep = _gauss_jordan(constrain_reduced)
            if verbose:
                print("Reference doublet: (", ii,",", jj,").")
                print("Orbit size:", Nequiv,". Number of constrains:",nconstrain,". Number of independent elements:", indep.shape[0])
                print(" ")
            # Save results of the elimination for this doublet:
            for iaux in range(9):
                for jaux in range(9):
                    kernel[nref2-1, iaux, jaux]  = kern[iaux,jaux]
            n_indep_fc[nref2-1] = indep.shape[0]
            for iaux in range(indep.shape[0]):
                indep_fc[nref2-1,iaux] = indep[iaux]
    
            if (nall2 == tot2): # All the doublets have been analyzed; exit nested loop
                break 
        else: # This is executed if the inner loop completes without breaking
            continue
        break
    
    tensor = np.zeros((nref2,max(norbit),9,max(n_indep_fc)),dtype=np.float64)
    for i in range(nref2):
        for j in range(norbit[i]):
            for indexprime in range(9): 
                for index in range(n_indep_fc[i]):
                    tensor[i,j,indexprime,index] = np.sum(R[i,j,indexprime,:] * kernel[i,:,index])

    end_time = time.time()
    execution_time = end_time - start_time

    if verbose:
        print(" ")
        print("Total doublets:", nall2)
        print("Number of Orbits:",nref2)
        print("Total number of independent elements:", sum(n_indep_fc[:nref2]))
        print(" ")

        print("Saving null-space and rotation matrices...") 
        print("execution_time in doublet recognition:", execution_time, " s")
        print("===== DOUBLET CLASSIFICATION and GJ ELIMINATION FINISHED ======")
        print(" ")

    sys.stdout.flush()

    return(orbit2a[:nref2], orbit2s[:nref2], norbit[:nref2], indep_fc[:nref2], n_indep_fc[:nref2], tensor[:nref2])

def recognize_triplet(dyn, orbit1a, map_uc, symprec=1e-5, verbose=False):
    """
    Classifies atomic triplets, and returns the atomic map for each symmetry operation,
    as well as the required information to reconstruct the 3rd order FCs from the symmetry independent
    triplets.

    Parameters
    ----------
        - dyn: object
            Cellconstructor dynamical matrix.
        - orbit1a: np.ndarray
            Atomic mapping for each symmetry operation. Dimension [Natom_sc,Nsym].
        - map_uc: An array that works as a tool to know which index i we need to consider, taking
            index j as reference (unit-cell). Dimension [Natom_sc, Natom_sc].
        - symprec: float
            Tolerance parameter for spglib in symmetry detection.
            Defaults to 1e-5.
        - verbose: bool
            If True prints information during execution.
            Defaults to False.

    Returns
    -------
        - orbit3a: np.ndarray
            Triplet mapping for each symmetry operation. Dimension [Nref3,tbd,3].
            In the second axis, valid inforation is up to the norbit[ref2]-th element.
        - orbit3s: np.ndarray
            Symmetry mapping for each doublet mapping. Dimension [Nref3,tbd,2].
            In the second axis, valid inforation is up to the norbit2[ref2]-th element.
        - norbit: np.ndarray
            Number of equivalent doublets for each of the orbits. Dimension [Nref3].
        - indep_fc: np.ndarray
            Which are the independent cartesian elements from the 3x3x3 FC matrix. Dimension [Nref3,27].
            In the second axis, valid information is up to the n_indep_fc[ref2]-th element.
        - n_indep_fc: np.ndarray
            Number of independent cartesian elements for each reference triplet. Dimension [Nref3].
        - kernel: np.ndarray
            Kernel from which we will reconstruct the 3x3x3 FC matrix, knowing just the independet
            elements. Dimension [Nref3, 27, 27]
        - Rot: np.ndarray
            Contains all the third order cartesian symmetry operations for all permutions.
            Dimension [6,Nsym,27,27]
        - mapping_triplet: np.ndarray
            Maps an arbitrary triplet to its reference triplet, permutaion, and symmtry operation after
            the classification. Dimension [Natom,Natom_sc,Natom_sc,3]
    """
    start_time = time.time()

    if verbose:
        print("===== STARTING TRIPLET CLASSIFICATION =====")

    # Get the symmetries of the unit cell
    spg_syms = spglib.get_symmetry(dyn.structure.get_spglib_cell(), symprec)
    sym_uc = cellconstructor.symmetries.GetSymmetriesFromSPGLIB(spg_syms, regolarize= False) #In crystalline coord respect uc (trans!)
    nsym = len(sym_uc)
    
    # Obtain the Sc size (list of 3 integers)
    sc_size = dyn.GetSupercell()
    nsupercell = np.prod(sc_size)
    
    # Obtain point group symmetries in cryst coord respect sc: 
    sym_list = np.zeros((nsym,3,4),dtype=np.float64)
    for isym in range(nsym):
        sym_list[isym,:,:] = sym_uc[isym]
        for ll in range(3): # Transl respect sc
            sym_list[isym,ll,3] = sym_list[isym,ll,3]/sc_size[ll]
    
    # Create an object from the Phonons class of the SC
    dyn_supercell = dyn.GenerateSupercellDyn(sc_size)
    
    # Get the symmetries of the supercell
    spg_syms_sc = spglib.get_symmetry(dyn_supercell.structure.get_spglib_cell(), symprec)
    sym_list_sc = cellconstructor.symmetries.GetSymmetriesFromSPGLIB(spg_syms_sc, regolarize= False)
    
    # Obtain symmetries that are pure translations:
    translation = np.zeros((nsupercell,3,4),dtype=np.float64)
    for i in range(nsupercell):
        translation[i,:,:] = sym_list_sc[i*nsym]
    
    
    #Obtain the rotations in cartesian coord:
    rot_cart = np.zeros((nsym,3,3),dtype=np.float64)
    for isym in range(nsym):
        rot_cart[isym] = cellconstructor.Methods.convert_matrix_cart_cryst2(sym_list[isym,:,:3], dyn.structure.unit_cell, cryst_to_cart = True)
        # Set elements smaller than the threshold to 0.0
        rot_cart[isym][np.abs(rot_cart[isym]) < 1e-12] = 0.0
    
    # Number of atoms
    nat = dyn.structure.N_atoms
    ntot = dyn_supercell.structure.N_atoms
    
    
    permutations=np.array([[0,1,2], [1,0,2], [2,1,0], [0,2,1], [1,2,0], [2,0,1]], dtype=np.intc)
    
    # Initial size of array
    array_size = 16
    
    #Reference triplets:
    ref3 = np.zeros((array_size,3),dtype=np.intc)
    nref3 = 0

    #List of all triplets:
    tot3=nat*ntot**2
    nall3=0
    
    #Create Rotation tensor: R_{αβγ}^{α'β'γ'}= R_{α}^{α'}*R_{β}^{β'}*R_{γ}^{γ'}
    Rot = np.zeros((6,nsym,27,27),dtype=np.float64) # Following Eq. (45), 2nd index refers to a single index combining α'β'γ'
    
    #In order to take into account perm of indexes:
    cart_index = np.zeros(3,dtype=np.intc)
    for iperm in range(6):
        for isym in range(nsym):
            for alphaprime in range(3):
                for betaprime in range(3):
                    for gammaprime in range(3):
                        indexprime = (3*alphaprime+betaprime)*3+gammaprime
                        for alpha in range(3):
                            cart_index[0] = alpha
                            for beta in range(3):
                                cart_index[1] = beta
                                for gamma in range(3):
                                    cart_index[2] = gamma
                                    index = (3*alpha+beta)*3+gamma
                                    alphaperm=cart_index[permutations[iperm,0]]
                                    betaperm =cart_index[permutations[iperm,1]]
                                    gammaperm=cart_index[permutations[iperm,2]]
                                    Rot[iperm,isym,indexprime,index] = rot_cart[isym,alphaprime,alphaperm]*rot_cart[isym,betaprime,betaperm]*rot_cart[isym,gammaprime,gammaperm]
    

    #Matrix in which we apply Gauus-Jordan: M = R-I :
    M=Rot.copy()
    nontrivial=np.zeros((6,nsym,27),dtype=np.intc)
    
    for iperm in range(6):
        for isym in range(nsym):
            for indexprime in range(27):
                M[iperm,isym,indexprime,indexprime]-=1.0 #Possible values in the diagonal of M: 0,-1,-2
                for index in range(27):
                    if abs(M[iperm,isym,indexprime,index])>1e-12: # As a trivial constrain is a row of zeroes in R2
                        nontrivial[iperm,isym,indexprime]=1
                    else:
                        M[iperm,isym,indexprime,index]=0.0

    nref3 = SCHAModules.module_hess.get_nref3(nat,ntot,tot3,nsym,orbit1a,map_uc)
    orbit3a, orbit3s, norbit, indep_fc, n_indep_fc, kernel, mapping_triplet = SCHAModules.module_hess.recognize_triplet(nat,ntot,tot3,nref3,nsym,orbit1a,map_uc,nontrivial,M,verbose)

    end_time = time.time()
    execution_time = end_time - start_time

    # Print detailed info
    if verbose:
        print(" ")
        print("Total triplets:", tot3)
        print("Number of Orbits:", nref3)
        print("Total number of independent elements:", sum(n_indep_fc[:nref3]))
        print(" ")


        print("execution_time in triplet recognition:", execution_time, " s")
        print("===== TRIPLET CLASSIFICATION and GJ ELIMINATION FINISHED ======")
        print(" ")

    sys.stdout.flush()
    
    return(orbit3a[:nref3], orbit3s[:nref3], norbit[:nref3], indep_fc[:nref3], n_indep_fc[:nref3], kernel[:nref3], Rot, mapping_triplet)

def recognize_quadruplet(dyn, orbit1a, map_uc, verbose=False, symprec=1e-5):
    """
    Classifies atomic quadruplets, and returns the atomic map for each symmetry operation,
    as well as the required information to reconstruct the 4th order FCs from the symmetry independent
    quadruplets.

    Parameters
    ----------
        - dyn: object
            Cellconstructor dynamical matrix.
        - orbit1a: np.ndarray
            Atomic mapping for each symmetry operation. Dimension [Natom_sc,Nsym].
        - map_uc: An array that works as a tool to know which index i we need to consider, taking
            index j as reference (unit-cell). Dimension [Natom_sc, Natom_sc].
        - symprec: float
            Tolerance parameter for spglib in symmetry detection.
            Defaults to 1e-5.
        - verbose: bool
            If True prints information during execution.
            Defaults to False.

    Returns
    -------
        - orbit4a: np.ndarray
            Triplet mapping for each symmetry operation. Dimension [Nref4,tbd,4].
            In the second axis, valid inforation is up to the norbit[ref2]-th element.
        - orbit4s: np.ndarray
            Symmetry mapping for each doublet mapping. Dimension [Nref4,tbd,2].
            In the second axis, valid inforation is up to the norbit4[ref4]-th element.
        - norbit: np.ndarray
            Number of equivalent doublets for each of the orbits. Dimension [Nref4].
        - indep_fc: np.ndarray
            Which are the independent cartesian elements from the 3x3x3x3 FC matrix. Dimension [Nref4,81].
            In the second axis, valid information is up to the n_indep_fc[ref4]-th element.
        - n_indep_fc: np.ndarray
            Number of independent cartesian elements for each reference triplet. Dimension [Nref4].
        - kernel: np.ndarray
            Kernel from which we will reconstruct the 3x3x3x3 FC matrix, knowing just the independet
            elements. Dimension [Nref4, 81, 81]
        - Rot: np.ndarray
            Contains all the third order cartesian symmetry operations for all permutions.
            Dimension [24,Nsym,81,81]
        - mapping_quadruplets: np.ndarray
            Maps an arbitrary triplet to its reference triplet, permutaion, and symmtry operation after
            the classification. Dimension [Natom,Natom_sc,Natom_sc,3]
    """

    start_time = time.time()

    if verbose:
        print("===== STARTING QUADRUPLET CLASSIFICATION =====")

    # Get the symmetries of the unit cell
    spg_syms = spglib.get_symmetry(dyn.structure.get_spglib_cell(), symprec)
    sym_uc = cellconstructor.symmetries.GetSymmetriesFromSPGLIB(spg_syms, regolarize= False) #In crystalline coord respect uc (trans!)
    nsym = len(sym_uc)
    
    # Obtain the Sc size (list of 3 integers)
    sc_size = dyn.GetSupercell()
    nsupercell = np.prod(sc_size)
    
    # Obtain point group symmetries in cryst coord respect sc: 
    sym_list = np.zeros((nsym,3,4),dtype=np.float64)
    for isym in range(nsym):
        sym_list[isym,:,:] = sym_uc[isym]
        for ll in range(3): # Transl respect sc
            sym_list[isym,ll,3] = sym_list[isym,ll,3]/sc_size[ll]
    
    # Create an object from the Phonons class of the SC
    dyn_supercell = dyn.GenerateSupercellDyn(sc_size)
    
    # Get the symmetries of the supercell
    spg_syms_sc = spglib.get_symmetry(dyn_supercell.structure.get_spglib_cell(), symprec)
    sym_list_sc = cellconstructor.symmetries.GetSymmetriesFromSPGLIB(spg_syms_sc, regolarize= False)
    
    # Obtain symmetries that are pure translations:
    translation = np.zeros((nsupercell,3,4),dtype=np.float64)
    for i in range(nsupercell):
        translation[i,:,:] = sym_list_sc[i*nsym]
    
    
    #Obtain the rotations in cartesian coord:
    rot_cart = np.zeros((nsym,3,3),dtype=np.float64)
    for isym in range(nsym):
        rot_cart[isym] = cellconstructor.Methods.convert_matrix_cart_cryst2(sym_list[isym,:,:3], dyn.structure.unit_cell, cryst_to_cart = True)
        # Set elements smaller than the threshold to 0.0
        rot_cart[isym][np.abs(rot_cart[isym]) < 1e-12] = 0.0
    
    # Number of atoms
    nat = dyn.structure.N_atoms
    ntot = dyn_supercell.structure.N_atoms

    tot4=nat*ntot**3
    
    #Create Rotation tensor: R_{αβγ}^{α'β'γ'}= R_{α}^{α'}*R_{β}^{β'}*R_{γ}^{γ'}
    Rot = np.zeros((24,nsym,81,81),dtype=np.float64) # Following Eq. (45), 2nd index refers to a single index combining α'β'γ'
    #In order to take into account perm of indexes:
    Rot = SCHAModules.module_hess.generate_rot4(rot_cart)
    #Matrix in which we appllly Gauus-Jordan: M = R-I :
    M=Rot.copy()
    nontrivial=np.zeros((24,nsym,81),dtype=np.intc)

    for iperm in range(24):
        for isym in range(nsym):
            for indexprime in range(81):
                M[iperm,isym,indexprime,indexprime]-=1.0 #Possible values in the diagonal of M: 0,-1,-2
                for index in range(81):
                    if abs(M[iperm,isym,indexprime,index])>1e-12: # As a trivial constrain is a row of zeroes in R2
                        nontrivial[iperm,isym,indexprime]=1
                    else:
                        M[iperm,isym,indexprime,index]=0.0
    
    nref4 = SCHAModules.module_hess.get_nref4(nat,ntot,tot4,nsym,orbit1a,map_uc)
    orbit4a, orbit4s, norbit, indep_fc, n_indep_fc, kernel, mapping_quadruplet = SCHAModules.module_hess.recognize_quadruplet(nat,ntot,tot4,nref4,nsym,orbit1a,map_uc,nontrivial,M,verbose)

    end_time = time.time()
    execution_time = end_time - start_time

    if verbose:
        print(" ")
        print("Total quadruplets:", tot4)
        print("Number of Orbits:",nref4)
        print("Total number of independent elements:", sum(n_indep_fc[:nref4]))
        print(" ")

        print("execution_time in quadruplet recognition:", execution_time, " s")
        print("===== QUADRUPLET CLASSIFICATION and GJ ELIMINATION FINISHED ======")
        print(" ")
    return(orbit4a[:nref4], orbit4s[:nref4], norbit[:nref4], indep_fc[:nref4], n_indep_fc[:nref4], kernel[:nref4], Rot, mapping_quadruplet)

# Small functions required in main subroutines.

def _gauss_jordan(a):
    """
    ShengBTE: Specialized version of Gaussian elimination.
    """

    EPS =1e-10
    row=a.shape[0]
    col=a.shape[1]

    dependent=np.empty(col,dtype=np.intc)
    independent=np.empty(col,dtype=np.intc)
    b=np.zeros((col,col),dtype=np.double)

    irow=0
    ndependent=0
    nindependent=0
    for k in range(min(row,col)):
        for i in range(row):
            if abs(a[i,k])<EPS:
                a[i,k]=0.
            # else:
            #     print("First if:",k+1,i+1)
        for i in range(irow+1,row):
            if abs(a[i,k])-abs(a[irow,k])>EPS:
                for j in range(k,col):
                    tmp=a[irow,j]
                    a[irow,j]=a[i,j]
                    a[i,j]=tmp
        if abs(a[irow,k])>EPS:
            dependent[ndependent]=k
            ndependent+=1
            for j in range(col-1,k,-1):
                a[irow,j]/=a[irow,k]
            a[irow,k]=1.
            for i in range(row):
                if i==irow:
                    continue
                for j in range(col-1,k,-1):
                    a[i,j]-=a[i,k]*a[irow,j]/a[irow,k]
                a[i,k]=0.
            if irow<row-1:
                irow+=1
        else:
            independent[nindependent]=k
            nindependent+=1
    for j in range(nindependent):
        for i in range(ndependent):
            b[dependent[i],j]=-a[i,independent[j]]
        b[independent[j],j]=1.
    return (b,independent[:nindependent])

def _map_unitcell(dyn_supercell, T_list):
    """
    Create a map of atomic indeces after the application of the Translation operation that bring atoms in the SC (ii) to the unit cell
    """
    ntot = dyn_supercell.structure.N_atoms
    nsupercell = T_list.shape[0]
    nat = ntot/nsupercell

    map_uc = np.zeros((ntot,ntot),dtype=np.intc)
    map_tr = np.zeros((ntot), dtype=np.intc)
    translation = np.empty([len(T_list),3,4], dtype=np.float64)
    for i in range(len(T_list)):
        translation[i,:,:3] = np.identity(3)
        translation[i,:,3] = T_list[i]

    for ii in range(ntot):
        uc_index=ii%nat
        #Search the corresponding translation
        for j in range(nsupercell):
            shifted_struct  = dyn_supercell.structure.copy()
            shifted_struct.apply_symmetry(translation[j],delete_original=True)
            irt_shifted = np.array(shifted_struct.get_equivalent_atoms(dyn_supercell.structure), dtype =np.intc)
            index_shifted = irt_shifted[ii]
            if (index_shifted == uc_index): # If after T, the first atom is in the 1uc:
                for jj in range(ntot):
                    map_uc[ii,jj] = irt_shifted[jj]
                    map_tr[ii] = j
                break
    return map_uc, map_tr

def _doublet_in_list(doublet, llist, nlist):
    """
    Return True if doublet is found in llist[:,:nlist]. 
    """
    return any(np.all(doublet == llist[:nlist, :], axis=1))
    # axis = 1 along rows
    # any() function returns True if any element of an iterable is True

def _doublets_are_equal(doublet1,doublet2):
    """
    Return True if two doublets are equal and False otherwise.
    """
    return np.all(doublet1 == doublet2)  # np.all() checks if all elements in a given boolean array are True