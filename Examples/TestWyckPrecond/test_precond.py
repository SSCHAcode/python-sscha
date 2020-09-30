from __future__ import print_function
from __future__ import division
"""
This script test whether the wyckoff preconditioner is actually the inverse of the
Dynamical matrix or not
"""

import numpy as np

import cellconstructor as CC
import cellconstructor.Phonons
import cellconstructor.symmetries

import sscha, sscha.SchaMinimizer
import sys,os, pytest


def check_sum_rule():
    
    # Change to local directory (used to test the file)
    total_path = os.path.dirname(os.path.abspath(__file__))
    os.chdir(total_path)
    
    # Load the dynamical matrix
    dyn = CC.Phonons.Phonons("dyn")
    #dyn.structure.masses["O"] = dyn.structure.masses["H"]
    dyn.Symmetrize()

    # Lets try to manually compute the preconditioning
    w,pols =  dyn.DyagDinQ(0)

    _m_ = dyn.structure.get_masses_array()
    m = np.tile(_m_, (3,1)).T.ravel()
    msqrt = np.sqrt(m)
    m_sqrt = np.tile(msqrt, (3,1)).T

    trans = CC.Methods.get_translations(pols, _m_)
    p_trans = pols[:, trans]
    pols = pols[:, ~trans]
    w = w[~trans]

    # Get the whole metric of trans pols
    standard_metric = np.einsum("am, an", p_trans, p_trans)
    divided_metric = np.einsum("am, an", p_trans / m_sqrt, p_trans / m_sqrt)
    multiply_metric = np.einsum("am,an", p_trans * m_sqrt, p_trans * m_sqrt)

    print("Standard:")
    print(standard_metric)
    print("Divided metric:")
    print(divided_metric)
    print("Multiply metric:")
    print(multiply_metric)

    # Check now the block of how translations interacts with all the others
    m_sqrt1 = np.tile(msqrt, (len(w) , 1)).T
    
    # Get the whole metric of trans pols
    standard_metric = np.einsum("am, an", p_trans, pols)
    divided_metric = np.einsum("am, an", p_trans / m_sqrt, pols / m_sqrt1)
    multiply_metric = np.einsum("am,an", p_trans * m_sqrt, pols * m_sqrt1)
    mixed_metric = np.einsum("am,an", p_trans * m_sqrt, pols / m_sqrt1)

    print("Standard:")
    print(standard_metric)
    print("Divided metric:")
    print(divided_metric)
    print("Multiply metric:")
    print(multiply_metric)
    print("Mixed metric:")
    print(multiply_metric)


def test_precond_wyck():
    # Change to local directory (used to test the file)
    total_path = os.path.dirname(os.path.abspath(__file__))
    os.chdir(total_path)
    
    # Load the dynamical matrix
    dyn = CC.Phonons.Phonons("dyn")
    #dyn.structure.masses["O"] = dyn.structure.masses["H"]
    dyn.Symmetrize()

    # Get the preconditioner
    prec = sscha.SchaMinimizer.GetStructPrecond(dyn)

    # Assert hermitianity
    assert np.max(np.abs(prec - prec.T)) < 1e-8, "Preconditioner not hermitian"

    # Multiply the two matrix
    identity = np.real(prec.dot(dyn.dynmats[0]))

    print("the dyn recostruction")

    # Lets try to manually compute the preconditioning
    w,pols =  dyn.DyagDinQ(0)

    _m_ = dyn.structure.get_masses_array()
    m = np.tile(_m_, (3,1)).T.ravel()
    msqrt = np.sqrt(m)
    
    # Discard translations
    print("SHAPE M:", _m_.shape)
    print("SHAPE Pols:", pols.shape)
    trans = CC.Methods.get_translations(pols, _m_)
    pols = pols[:, ~trans]
    w = w[~trans]

    Ddyn = np.einsum("am,bm, m", pols, np.conj(pols), w**2)
    inv_Ddyn = np.einsum("am, bm, m",  pols, np.conj(pols), 1/w**2)

    
    
    identity_1 = Ddyn.dot(inv_Ddyn)
    eigvals, _ = np.linalg.eigh(identity_1)
    print(eigvals)
    
    assert np.max(np.abs(eigvals[3:] - 1)) < 1e-8

    # Now lets get bach to the full space
    Phi = Ddyn * np.outer(msqrt, msqrt)
    Phi_inv = inv_Ddyn / np.outer(msqrt, msqrt)
    identity_2 = Phi.dot(Phi_inv)

    Phi2 = Phi.copy()
    CC.symmetries.CustomASR(Phi2)
    assert np.max(np.abs(Phi - Phi2)) < 1e-8, "Phi violated ASR"
    
    Phi_inv2 = Phi_inv.copy()
    CC.symmetries.CustomASR(Phi_inv2)
    assert np.max(np.abs(Phi_inv - Phi_inv2)) < 1e-8, "Phi inv violated ASR"

    print(Phi - dyn.dynmats[0])
    assert np.max(np.abs(Phi - dyn.dynmats[0])) < 1e-8

    identity_3 = identity_2.copy()
    CC.symmetries.CustomASR(identity_2)

    assert np.max(np.abs(identity_2 - identity_3)) < 1e-8, "Dot product violated ASR"
    
    eigvals, _ = np.linalg.eigh(identity_2)
    print(eigvals)
    
    assert np.max(np.abs(eigvals[3:] - 1)) < 1e-8


    
    

    print(identity)

    # Check whether the eigenvalues of the identity are close to 1
    eigvals, _ = np.linalg.eigh(identity)
    print(eigvals)
    eigvals = eigvals[ eigvals > 0.1 ]
    
    print (np.sqrt(np.sum( (eigvals- 1)**2)))

if __name__ == "__main__":
    check_sum_rule()
    #test_precond_wyck()
