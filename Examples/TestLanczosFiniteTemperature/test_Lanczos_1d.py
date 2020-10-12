from __future__ import print_function
import sscha
import sscha.DynamicalLanczos
import numpy as np
import cellconstructor as CC
import cellconstructor.Units
import itertools
import matplotlib.pyplot as plt
import sys, os

import sscha, sscha.DynamicalLanczos


# Define the system
W_S =  np.array([100, 150])
N_MODES = len(W_S)

d3 = np.zeros((N_MODES, N_MODES, N_MODES))
#d3[0,0,0] = 230000
d3[1,1,0] = d3[1,0,1] = d3[0,1,1] = 400000
d3[1,0,0] = d3[0,0,1] = d3[0,1,0] = 700000
#d3[2,0,0] = d3[0,0,2] = d3[0,2,0] = 100000

# Temperature
TEMP = 10000# 800000.

def print_matrix(x):
    n1,n2 = np.shape(x)

    for i in range(n1):
        for j in range(n2):
            sys.stdout.write(" %10.3e " % x[i,j])
        sys.stdout.write("\n")

def green_function(w_array, lanc, smearing = 1):
    """
    Compute the green function with the bubble diagram
    """

    # Occupation factors
    nt = 0 * lanc.w
    if lanc.T > 0:
        nt = 1 / (np.exp(lanc.w / (CC.Units.K_B * lanc.T)) - 1) 


    G = np.zeros(len(w_array), dtype = np.complex128)

    G_total = np.zeros((len(w_array), lanc.n_modes, lanc.n_modes), dtype = np.complex128)

    for a,h in itertools.product(range(lanc.n_modes), range(lanc.n_modes)):
        sigma = np.zeros(len(w_array), dtype = np.complex128)
        if a == h:
            sigma[:] = lanc.w[a]**2 - w_array**2 + 1j*w_array*smearing

        for b, c in itertools.product(range(lanc.n_modes), range(lanc.n_modes)):
            bubble = d3[a, b, c] * d3[b,c,h]
            bubble *=  -1. /(4*lanc.w[b] * lanc.w[c])
            Fmu = (lanc.w[b] + lanc.w[c]) * (1 + nt[b] + nt[c]) / ( (lanc.w[b] + lanc.w[c])**2 - w_array**2 + 1j*smearing * w_array)
            Fmu -= (nt[c] - nt[b])*(lanc.w[c] - lanc.w[b]) / ( (lanc.w[c] - lanc.w[b])**2 - w_array**2 + 1j*smearing*w_array)
            bubble *= Fmu

            #print("Bubble:",bubble[0])

            sigma[:] += bubble

        G_total[:,a,h] = sigma

    for k in range(len(w_array)):
        G_good = np.linalg.inv(G_total[k, :, :])
        
        G[k] = np.trace(G_good)

    return G
            


def test_lanczos_2modes(verbose = False):
    # Change to the local directory
    total_path = os.path.dirname(os.path.abspath(__file__))
    os.chdir(total_path)


    lanc = sscha.DynamicalLanczos.Lanczos()

    lanc.w = W_S 
    #lanc.w = np.array([100])
    lanc.n_modes = len(lanc.w)

    lanc.T = TEMP
    lanc.ignore_v4 = True
    lanc.initialized = True # Avoid symmetry initialization


    L = lanc.get_full_L_operator(debug_d3 = d3)
    Lt = lanc.get_full_L_operator_FT(debug_d3 = d3)

            
            

    # exclude_exchange = np.array([0,1,2,3,5])
    # Lt_reduces = np.zeros((len(exclude_exchange), len(exclude_exchange)), dtype = np.double)
    # for i,j in enumerate(exclude_exchange):
    #     Lt_reduces[i,:] = Lt[j, exclude_exchange]

    N, _ = np.shape(L)
    if verbose:
        print(" L = ")
        print_matrix(L)
        print()
        print(" Lt = ")
        print_matrix(Lt)
        print()
        print(" Linv = ")
        print_matrix(np.linalg.inv(L))
        print()
        print(" Lt inv = ")
        print_matrix(np.linalg.inv(Lt))


    w_array = np.linspace(0, 400, 500)
    eps = 10
    Gw = np.zeros(len(w_array), dtype = np.double)
    Gr = np.zeros(len(w_array), dtype = np.double)
    Gw_T0 = np.zeros(len(w_array), dtype = np.double)

    for i,w in enumerate(w_array):
        G_all = np.linalg.inv(Lt - np.eye(Lt.shape[0]) * (w**2 - 1j*w*eps))
        Gtmp = np.linalg.inv(L - np.eye(L.shape[0]) * (w**2 - 1j*w*eps))



        Gr[i] = np.real(G_all[0,0])
        for k in range(lanc.n_modes):
            Gw[i] += - np.imag(G_all[k,k]) #+ G_all[1,1]
            Gw_T0[i] += -np.imag(Gtmp[k,k])

    # Get the green function with the Lanczos algorithm
    spectral_lanczos_ft = np.zeros(w_array.shape, dtype= np.double)
    for k in range(lanc.n_modes):

        # Reset the Lanczos
        lanc.reset()

        # Prepare the perturbation along the phonon mode
        lanc.psi = np.zeros(lanc.L_linop.shape[0], dtype = np.double)
        lanc.psi[k] = 1

        # Execute the biconjugate Lanczos algorithm
        # We use 1000 steps, but it will converge much faster
        lanc.run_FT(1000, save_dir = None)

        arnoldi_matrix = lanc.build_lanczos_matrix_from_coeffs()

        if verbose:
            print("The arnoldi matrix:")
            print_matrix(arnoldi_matrix)
            print("Eigenvalues:")
            w2 = np.linalg.eigvals(arnoldi_matrix)
            print(np.sign(w2) * np.sqrt(np.abs(w2)))
            print("All eigenvalues:")
            w2 = np.linalg.eigvals(Lt)
            print(np.sort(np.sign(w2) * np.sqrt(np.abs(w2))))


        # Get the green function of the mode
        # with the continued fraction
        gf_kk = lanc.get_green_function_continued_fraction(w_array,
                                                           use_terminator = False,
                                                           smearing = eps / 2)
        spectral_lanczos_ft += -np.imag(gf_kk)

    G_bianco = -np.imag(green_function(w_array, lanc, eps))

    nan_mask = np.isnan(G_bianco)
    nan_mask = (nan_mask | np.isnan(Gw))
    
    thr = np.max(np.abs(G_bianco[~nan_mask] - Gw[~nan_mask]))
    maxval =  np.max(G_bianco[~nan_mask])
    if verbose:
        plt.plot(w_array, Gw, label = "Lanczos")
        plt.plot(w_array, Gw_T0, label = "Lanc T=0", ls = "dashed")
        plt.plot(w_array, spectral_lanczos_ft, label = "Lanczos continued frac.",
                 ls = "-.")
        plt.plot(w_array, G_bianco, ls  = "dotted", label = "Many body", color =  "k")

        plt.legend()
        plt.tight_layout()

        print("Good until: {} | max_value: {}".format(thr, maxval))
        plt.show()

    assert thr < 1e-8 * maxval



if __name__ == "__main__":
    test_lanczos_2modes(verbose = True)
