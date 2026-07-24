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

def map_singlet(q_list, q_list_frac, rcell, rot_cart):
    """
    This functions performs the mapping between the commensurate wave-vectors of the
    supercell. 
    
    Parameters
    ----------
        - q_list: np.ndarray
            List of wave-vectors in cartesian coordinates. Dimension [Nq,3].
        - q_list_frac: np.ndarray
            List of wave-vectors in fractional coordinates. Dimension [Nq,3].
        - rcell: np.ndarray
            Reciprocal unit-cell. Dimension [3,3].
        - rot_cart: np.ndarray
            Point symmetry in cartesian coordinates. Dimension [Nsym,3,3].
    
    Returns
    -------
        - mapping: np.ndarray
            The mapping between symmetry related wave-vectors for each point
            point symmetry of the crystal. Dimension [Nq,Nsym]
        - orbit1a: np.ndarray
            Wave-vector classification in orbits. The classification is not
            unique in this case, because needs to be used to construct the P
            matrix (symmetry related elements appear more than once).
            Dimension [Nq,tbd,1].
        - orbit1s: np.ndarray
            Which symmetry makes the classification in orbits/stars of orbit1a.
            Dimension [Nq,tbd,1].
        - norbit: np.ndarray
            Number of wave-vectors in orbit. Dimension [Nrefq1].
        - its_zb: np.ndarray
            Checks whether the points in q_list are on zone border or not. This
            is later used to impose time-reversal symmetry.
    """

    q_list_frac_fixed = np.empty(q_list_frac.shape)
    its_zb = np.empty(q_list_frac.shape[0], dtype=np.int8)
    for qi,q in enumerate(q_list):
        q_list_frac_fixed[qi] = np.round(CC.Methods.cart_to_cryst(rcell, _map_q_to_1st_bz(rcell, q)),6)
        if np.all((np.abs(np.round(2*q_list_frac_fixed[qi], 6))%1)<1e-6):
            its_zb[qi] = 0
        else:
            for alpha in range(3):
                if np.abs(q_list[qi,alpha])>1e-6:
                    if q_list[qi,alpha] > 1e-6:
                        its_zb[qi] = 1 # Non zone border. Positive class.
                        break
                    else:
                        its_zb[qi] = 2 # Non zone border. Negative class.
                        break

    mapping = np.zeros([len(q_list), rot_cart.shape[0]], dtype=np.int32)
    for qi, q in enumerate(q_list):
        for isym in range(rot_cart.shape[0]):
            q_sym_cart = rot_cart[isym] @ q
            q_sym_cart_1bz = _map_q_to_1st_bz(rcell, q_sym_cart)
            q_sym = np.round(CC.Methods.cart_to_cryst(rcell, q_sym_cart_1bz),6)
            match = np.all(np.abs(q_list_frac_fixed-q_sym)<1e-3, axis=1)
            qii = np.where(match)
            mapping[qi,isym] = qii[0][0]

    orbit1a = np.zeros([len(q_list),rot_cart.shape[0],1], dtype=np.int32)
    orbit1s = np.zeros([len(q_list),rot_cart.shape[0],1], dtype=np.int32)

    # Loop over all q-points, and knowing the mapping build the star for each q-point
    # All the classifications are considered (equivalents are not excluded in next iters)

    for qi in range(len(q_list)):
        for isym,qsym in enumerate(mapping[qi,:]):
            orbit1a[qi,isym,0] = qsym
            orbit1s[qi,isym,0] = isym
    return mapping, orbit1a, orbit1s, its_zb

def recognize_doublet(q_list, mapping, verbose=False):
    """
    This function classifies wave-vector doublets in orbits.

    Parameters
    ----------
        - q_list: np.ndarray
            List of wave-vectors in cartesian coordinates. Dimension [Nq,3].
        - mapping: np.ndarray
            The mapping between symmetry related wave-vectors for each point
            point symmetry of the crystal. Dimension [Nq,Nsym]
        - verbose: bool
            If True prints information during execution.
            Defaults to False.
    Returns
    -------
        - orbit2a: np.ndarray
            Wave-vector doublet classification in orbits. Dimension [Nrefq2,tbd,2].
        - orbit2s: np.ndarray
            Which symmetry makes the classification in orbits/stars of orbit2a.
            Dimension [Nrefq2,tbd,2].
        - norbit: np.ndarray
            Number of wave-vectors in orbit. Dimension [Nrefq2].
        - nref2: int
            Number of reference q-doublets.
    """

    start_time = time.time()

    if verbose:
        print("===== STARTING Q-DOUBLET CLASSIFICATION =====")

    doublet = np.empty([2], dtype=np.int32)
    doublet_sym = np.empty([2], dtype=np.int32)
    doublet_perm = np.empty([2], dtype=np.int32)

    permutations = [[0,1],[1,0]]

    orbit2t = np.zeros([len(q_list)**2,len(q_list)**2,2], dtype=np.int32)
    orbit2o = np.zeros([len(q_list)**2,len(q_list)**2,2], dtype=np.int32)
    equilist=np.zeros([len(q_list)**2,2], dtype=np.int32)
    norbit2 = np.zeros([len(q_list)**2], dtype=np.int32)
    nref2 = 0
    nall2 = 0
    for qmu in range(len(q_list)):
        for qnu in range(len(q_list)):
            if qmu==0 and qnu == 0:
                orbit2t[nref2,0,:] = [0,0]
                orbit2o[nref2,0,:] = [0,0]
                nref2 += 1
                nall2 += 1
                norbit2[nref2-1] += 1
            else:
                doublet[:] = [qmu, qnu]
                if not _doublet_in_list(doublet,equilist,nall2):
                    equilist[nall2] = doublet
                    orbit2t[nref2,0,:] = doublet
                    orbit2o[nref2,0,:] = [0,0]
                    nref2 += 1
                    nall2 += 1
                    norbit2[nref2-1] += 1
                else:
                    continue
                equiv=1
                for iperm in range(len(permutations)):
                    doublet_perm[0] = doublet[permutations[iperm][0]]
                    doublet_perm[1] = doublet[permutations[iperm][1]]
                    for isym in range(mapping.shape[1]):
                        doublet_sym = np.array([mapping[doublet_perm[0],isym], mapping[doublet_perm[1],isym]], dtype=np.int32)
                        if not _doublet_in_list(doublet_sym,equilist,nall2):
                            equilist[nall2] = doublet_sym
                            orbit2t[nref2-1,equiv,:]=doublet_sym
                            orbit2o[nref2-1,equiv,:]=[iperm,isym]
                            equiv += 1
                            nall2 += 1
                            norbit2[nref2-1] += 1
            if verbose:
                print(f"Reference q-doublet: ({qmu},{qnu})")
                print(f"Orbit size: {equiv}\n")
    
    end_time = time.time()
    execution_time = end_time - start_time
    
    if verbose:
        print(" ")
        print("Total q-doublets:", nall2)
        print("Number of Orbits:", nref2)
        print(" ")
        print("execution_time in q-doublet recognition:", execution_time, " s")
        print("===== Q-DOUBLET CLASSIFICATION FINISHED ======")
        print(" ")

    return orbit2t[:nref2], orbit2o[:nref2], norbit2[:nref2], nref2


def recognize_quadruplet(q_list, mapping, verbose=False):
    """
    This function classifies wave-vector quadruplets in orbits. It is writen
    if fortran. This function is just an interface.

    Parameters
    ----------
        - q_list: np.ndarray
            List of wave-vectors in cartesian coordinates. Dimension [Nq,3].
        - mapping: np.ndarray
            The mapping between symmetry related wave-vectors for each point
            point symmetry of the crystal. Dimension [Nq,Nsym]
        - verbose: bool
            If True prints information during execution.
            Defaults to False.
    Returns
    -------
        - orbit4a: np.ndarray
            Wave-vector doublet classification in orbits. Dimension [Nrefq4,tbd,4].
        - orbit4s: np.ndarray
            Which symmetry makes the classification in orbits/stars of orbit4a.
            Dimension [Nrefq4,tbd,2].
        - norbit: np.ndarray
            Number of wave-vectors in orbit. Dimension [Nrefq4].
        - nrefq4: int
            Number of reference q-quadruplets.
    """

    start_time = time.time()

    if verbose:
        print("===== STARTING Q-QUADRUPLET CLASSIFICATION =====")

    nrefq4, norbitq4 = SCHAModules.module_hess.get_q_nref4(mapping)
    orbit4a, orbit4s, norbit, nrefq4 = SCHAModules.module_hess.recognize_q_quadruplet(nrefq4, norbitq4, q_list, mapping,verbose)
    
    end_time = time.time()
    execution_time = end_time - start_time

    if verbose:
        print(" ")
        print("Total q-doublets:", len(q_list)**4)
        print("Number of Orbits:", nrefq4)
        print(" ")
        print("execution_time in q-doublet recognition:", execution_time, " s")
        print("===== Q-QUADRUPLET CLASSIFICATION FINISHED ======")
        print(" ")

    sys.stdout.flush()

    return orbit4a[:nrefq4], orbit4s[:nrefq4], norbit[:nrefq4], nrefq4

def construct_Pmn(mapping, orbit1t, orbit1o, pol_vecs, rot_cart):
    """
    This function construct the P matrix.
    
    Parameters
    ----------
        - mapping: np.ndarray
            The mapping between symmetry related wave-vectors for each point
            point symmetry of the crystal. Dimension [Nq,Nsym]
        - orbit1a: np.ndarray
            Wave-vector classification in orbits. The classification is not
            unique in this case(symmetry related elements appear more than once).
            Dimension [Nq,tbd,1].
        - orbit1s: np.ndarray
            Which symmetry makes the classification in orbits/stars of orbit1a.
            Dimension [Nq,tbd,1].
        - pol_vecs: np.ndarray
            Polarization vectors used to construct the P matrix.
            Dimension [Nq,Nmode,Nmode_sc]
        - rot_cart: np.ndarray
            Point symmetry in cartesian coordinates. Dimension [Nsym,3,3].
    Returns
    -------
        - Pmn: np.ndarray
            Pmn matrix. Dimension [Nq,Nmode,Nmode,Nsym]
    """
    # This should go in fortran if its computationally expensive
    nmodes = pol_vecs.shape[1]
    nat_sc = int(pol_vecs.shape[2]/3)
    nsym = rot_cart.shape[0]
    rot_pol_vec = np.empty([nat_sc,3], dtype=np.complex128)
    Pmn = np.empty([orbit1t.shape[0],nmodes,nmodes,nsym], dtype = np.complex128)
    for si, star in enumerate(orbit1t):
        for qsi in range(nsym):
            isym = orbit1o[si,qsi,0]
            for mu in range(nmodes):
                for nu in range(nmodes):
                    for a in range(nat_sc):
                        b = mapping[a,isym]
                        for alpha in range(3):
                            ref_pol_vec = np.reshape(pol_vecs[star[0],nu,:],[nat_sc,3])
                            rot_pol_vec[b, alpha] = np.matmul(rot_cart[isym,alpha,:],np.transpose(ref_pol_vec[a,:]))
                        Pmn[si,mu,nu,qsi] = np.matmul(np.conjugate(pol_vecs[star[qsi],mu,:]),np.transpose(np.reshape(rot_pol_vec, [nat_sc*3])))[0]

    return Pmn

def find_degeneracies(wq):
    """
    This function looks for vibrational modes with the same wave-vector modulation
    and degenerate energies.
    
    Parameters
    ----------
        - wq: np.ndarray
            Eigenvalue of the dynamical matrix. Dimension [Nq,Nmode]
    Returns
    -------
        - degs: np.ndarray
            If the i-th and j-th mode are degenerate for a given wave-vector q:
            degs[qi,i,j] == True. Dimension [Nq,Nmode,Nmode].
    """
    
    iq = wq.shape[0]
    nmodes = wq.shape[1]
    degs = np.zeros([iq,nmodes,nmodes], dtype=bool)
    for q in range(iq):
        for mu in range(nmodes):
            for nu in range(nmodes):
                if np.abs(wq[q,mu]-wq[q,nu])<1e-6:
                    degs[q,mu,nu]=True
    return degs

def _map_q_to_1st_bz(rcell, q_cart, atol=1e-5):
    """
    Maps q_cart to the 1st BZ.
    On boundary ties (e.g., equidistant points), it deterministically
    prefers the vector that is lexicographically smaller in Cartesian coordinates.
    """
    q_cart = np.array(q_cart).flatten()
    rcell = np.array(rcell)

    inv_rcell = np.linalg.inv(rcell)
    q_frac = q_cart @ inv_rcell
    q_frac_wrapped = q_frac - np.round(q_frac)

    shifts = np.array([-1, 0, 1])
    grid_frac = np.stack(np.meshgrid(shifts, shifts, shifts), -1).reshape(-1, 3)

    candidate_fracs = q_frac_wrapped + grid_frac
    candidate_carts = candidate_fracs @ rcell

    distances_sq = np.sum(candidate_carts**2, axis=1)
    min_dist_sq = np.min(distances_sq)

    # This tolerance must be wider than the 6-digit precision of the input
    is_close = np.abs(distances_sq - min_dist_sq) < atol
    tied_candidates = candidate_carts[is_close]

    idx = np.lexsort((tied_candidates[:, 2], tied_candidates[:, 1], tied_candidates[:, 0]))

    return tied_candidates[idx[0]]

def _doublet_in_list(doublet, llist, nlist):
    """
    Return True if doublet is found in llist[:,:nlist]. 
    """
    return any(np.all(doublet == llist[:nlist, :], axis=1))
    # axis = 1 along rows
    # any() function returns True if any element of an iterable is True