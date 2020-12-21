from __future__ import print_function
from __future__ import division

import sscha, sscha.DynamicalLanczos
import numpy as np

import scipy
import scipy.optimize


import matplotlib.pyplot as plt



def sscha_1d(func, x, m = 1, return_res = False):
    # Get the sscha value
    #_xreal_ = random.normal(0, 1, size = 5000)
    _xreal_ = np.linspace(-5, 5, 10000)
    def complete_sscha(phi):
        if phi < 0:
            return 1000000 + phi**2
        kin_energy = 1 / (8. * m * phi**2)
        _x_ = x + _xreal_ * phi
        dx_ = _x_[1]- _x_[0]
        
        #pot_energy = sum(func(_x_)) / 5000
        pot_energy = np.sum( func(_x_) * np.exp(-(_x_ - x)**2 / (2*phi**2))) *dx_ / np.sqrt(2*np.pi*phi**2)
        
        return kin_energy + pot_energy
    
    # Minimize
    res = scipy.optimize.minimize(complete_sscha, 1, method ="BFGS")
    #print res.x, complete_sscha(res.x[0]), 1 / (8. * m * res.x[0]**2)
    if return_res:
        return complete_sscha(res.x[0]), res.x[0]
    return  complete_sscha(res.x[0])

def sscha_rho(_x_, x0, phi0):
    return np.exp(-(_x_ - x0)**2 / (2*phi0**2))/ np.sqrt(2*np.pi*phi0**2)

def sscha_w(m, phi0):
    return 1 / (2*m * phi0**2)

def sscha_full(func, x0 = 1, phi0 = 1, m = 1):
    _xreal_ = np.linspace(-5, 5, 10000)
    def complete_sscha(params):
        x, phi = params
        if phi < 0:
            return 1000000 + phi**2
        kin_energy = 1 / (8. * m * phi**2)
        _x_ = x + _xreal_ * phi
        dx_ = _x_[1]- _x_[0]
        
        #pot_energy = sum(func(_x_)) / 5000
        rho = sscha_rho(_x_, x, phi)
        pot_energy = np.sum( func(_x_) * rho) * dx_
        
        return kin_energy + pot_energy
    
    start_par = [x0, phi0]
    res = scipy.optimize.minimize(complete_sscha, start_par, method ="BFGS")
    
    x0 = res.x[0]
    phi0= res.x[1]
    energy = complete_sscha([x0, phi0])
    return energy, x0, phi0

def sscha_force(x, x0, phi0, m = 1):
    return -(x - x0) * sscha_w(m, phi0)**2 * m


def bo_func(x):
    new_x = x
    return 3*new_x**4 - 3*new_x**2 + 0.5*new_x**3

def bo_force(x):
    new_x = x
    f = 12*new_x**3 - 6*new_x + 3 / 2. * new_x**2
    return -f

def bo_curvature(x):
    new_x = x
    return 36*new_x**2 - 6 + 3 * new_x

def bo_third_deriv(x):
    new_x = x
    return 72*new_x + 3


def get_d3_d4(x0, phi0):
    """
    Compute the d3 and d4 to compare between the exact Lanczos and the face lanczos application
    """
    
    d4 = 72
    d3 = 72 * x0 + 3

    return d3, d4


def get_numerical_d3_d4(x0, phi0):
    N_Random = 100000

    x_pos = np.random.normal(loc= x0, scale = phi0, size = N_Random) 
    u = x_pos - x0
    f = bo_force(x_pos) 
    f_scha = sscha_force(x_pos, x0, phi0) 

    w = sscha_w(1, phi0)

    Y = 2*w

    d3 = -Y**2 *  np.sum(u * u* (f - f_scha)) / N_Random
    d4 = -Y**3 * np.sum(u* u* u* (f - f_scha)) / N_Random


    return d3, d4



def get_L_matrix(x0, phi0):
    """
    Get the Lanczos L matrix
    """

    L = np.zeros((3,3), dtype = np.double)

    d3, d4 = get_d3_d4(x0, phi0)
    w = sscha_w(1, phi0)

    L[1,1] = - d4 / (2 * w) - 2*w**2
    L[1,2] = -4 * w**2 
    L[1,0] = 4*w*d3 
    L[0,1] = d3 / (8 * w**2)
    L[0,0] = -w**2

    return -L


def apply_L(x0, phi0, psi_vector, transpose = False):
    L = get_L_matrix(x0, phi0)
    if transpose:
        L = L.T
    return L.dot(psi_vector)


def get_pert_force_potential(x0, phi0, pert_x0, pert_Y):
    """
    Compute numerically the derivative of the average force and curvature of the BO
    by applying finite differences.
    """
    _x_ = np.linspace(-10, 10, 10000)
    dx = _x_[1] - _x_[0]

    PHI  = sscha_w(1, phi0) ** 2

    f_av = np.sum( (bo_force(_x_) - sscha_force(_x_, x0, phi0))  * sscha_rho(_x_, x0, phi0)) * dx
    c_av = np.sum( (bo_curvature(_x_) - PHI) * sscha_rho(_x_, x0, phi0)) * dx

    dphi = -phi0**3 * pert_Y / 2

    delta = 0.001

    # Check if the force both are zero?
    print("SSCHA GRADIENTS:")
    print("<f> = {}".format(f_av))
    print("<d^2v/dr^2 - Phi> = {}".format(c_av))

    new_x0 = x0 + pert_x0 * delta
    new_phi = phi0 + dphi * delta
    PHI_new =  sscha_w(1, new_phi)**2

    f_av_new = np.sum( (bo_force(_x_) - sscha_force(_x_, x0, phi0)) * sscha_rho(_x_, new_x0, new_phi)) * dx
    c_av_new = np.sum( (bo_curvature(_x_) - PHI) * sscha_rho(_x_, new_x0, new_phi)) * dx


    # Get the derivative 
    new_f = (f_av_new - f_av) / delta   
    new_curv = (c_av_new - c_av) / delta


    # Recompute everything with a monte carlo algorithm
    N = 50000
    random_x = np.random.normal(loc = x0, scale = phi0, size = N)
    weights = (sscha_rho(random_x, new_x0, new_phi) / sscha_rho(random_x, x0, phi0) - 1) / delta 
    f_pert_new  = np.sum(weights * (bo_force(random_x) - sscha_force(random_x, x0, phi0))) / N
    f_pert_new_noscha  = np.sum(weights * (bo_force(random_x))) / N
    c_pert_new  = np.sum(weights * (bo_curvature(random_x) - PHI))/ N 

    print("With monte carlo we get:")
    print("<f> pert = {}".format(f_pert_new))
    print("<f> pert = {} | noscha".format(f_pert_new_noscha))
    print("<d2v/dr2> pert = {}".format(c_pert_new))



    return new_f, new_curv


# This is the exact solution of the Lanczos algorithm (analytical in 1D) at T = 0
def DynamicalSSCHA(func,w_array, x0, phi0, _x_, full_spectral = False, m = 1, smearing = 1e-2):
    """
    Compute the dynamical sscha equation with the Lanczos algorithm
    
    Parameters
    ----------
        func :  function 
            the BO landscape
        w_array : ndarray
            The frequencies
        x0 : float
            The x0 of the equilibrium sscha
        phi0 : float
            The variance of the equilibrium sscha
        _x_ : ndarray
            The x array
        Full_spectral : bool
            Return all the poles without weight, or weight on a standard probe.
        m : float
            The mass
        smearing : float
            The smearing for the spectral function
            
    Returns
    -------
        spectral_function : ndarray (size = w.shape, dtype = complex128)
            The spectral function
    """
    
    # Get the SSCHA quantities
    w =  sscha_w(m, phi0)
    D_sscha = w**2
    
    B = 2* w**2
    Chi = 1 / (4. * w**3)
    
    sym_pal = np.sqrt(B * Chi)
    
    # Get the third and fourth order integrals
    _x_ = np.linspace(x0  -5 * phi0, x0 + 5*phi0, 10000)
    dx = _x_[1] - _x_[0]
    rho = sscha_rho(_x_, x0=x0, phi0 = phi0)
    
    
    v1 = np.gradient(func(_x_)) / dx
    v2 = np.gradient(v1) / dx
    v3 = np.gradient(v2) / dx
    v4 = np.gradient(v3) / dx
    
    #figure()
    #xlim((-1.5,1.5))
    #ylim((-2,2))
    #plot(_x_, v1)
    #plot(_x_, v2, label = "v2")
    #plot(_x_, v3, label = "v3")
    #plot(_x_, v4, label = "v4")
    #legend()
    #show()
    
    d1 = np.sum(rho * v1) * dx
    d2 = np.sum(rho * v2) * dx / m
    print("Check v1: {}".format(d1))
    print("Check v2:")
    print("D2 = {} | d2 = {}".format(D_sscha, d2))
    print("Test rho: {}".format(np.sum(rho) * dx))
    d3 = np.sum(rho * v3) * dx / m**2
    d4 = np.sum(rho * v4) * dx / m**3
    
    print("d3 = {}".format(d3))
    print("d4 = {}".format(d4))
    print("sym_pal = {} | test = {}".format(sym_pal, sqrt(1 / (2*w))))
    
    print("w2 = {} | w3 = {} | w4 = {}".format(d2, d3 * sym_pal, d4*sym_pal**2))
    #d4 = 0
    #d3 = 0
    
    # Prepare the Lanczos matrix
    Lanczos_mat = np.zeros((2,2))
    Lanczos_mat[0,0] = D_sscha
    Lanczos_mat[1,0] = - d3 * sym_pal 
    Lanczos_mat[0,1] = Lanczos_mat[1,0]
    Lanczos_mat[1,1] = 2 * B + d4 * sym_pal**2 
    
    
    # Dyagonalize the matrix
    eigvals, eigvects = np.linalg.eigh(Lanczos_mat)
    
    print("spectral poles: ", np.sqrt(eigvals))
    print("Eigvect 0:", eigvects[:,0])
    print("Eitvect 1:", eigvects[:,1])
    
    spectral = np.zeros(np.shape(w_array), dtype = np.complex128)
    if not full_spectral:
        for i in range(2):
            spectral += 1 / (w_array**2 - eigvals[i] - 1j*smearing)
    else:
        for i in range(2):
            mat_elem = eigvects[0,i]**2
            spectral += mat_elem / ((w_array -1j*smearing)**2 - eigvals[i])
    print("Freq real: {} | expected: {}".format(np.sqrt(-1/np.real(spectral[0])), 
                                                np.sqrt(eigvals[0])))
    print("Only first: {}".format(np.sqrt(eigvals[0]) / np.abs(eigvects[0,0])))
    return spectral
    
    


def test_lanczos_1d(plot = False):
    """
    This is the proper testing function.
    Here we are going to apply the SSCHA in the one dimensional BO energy landscape.
    Then, we compute the resutls of the Lanczos algorithm, comparing the exact one with
    the one implemented in the sscha module of python-sscha (that this subroutine tests).
    """

    # Get the SSCHA solution of the 1D problem
    energy, x_c, phi = sscha_full(bo_func)

    # Get the frequency
    w = sscha_w(m = 1, phi0 = phi)

    # Initialize the lanczos algorithm
    lanc = sscha.DynamicalLanczos.Lanczos()

    # Initliaize the lanczos from a 1D problem
    lanc.w = np.zeros(1, dtype = np.double)
    lanc.w[0] = w
    lanc.n_modes = 1
    lanc.pols = np.ones((1,1), dtype = np.double)
    lanc.m = np.array([np.double(1)])

    # get the N random configuration
    # for the Lanczos application
    N_RANDOM = 100000
    x_ensemble = np.double(np.random.normal(x_c, phi, size = (N_RANDOM, 1)))
    f = bo_force(x_ensemble) - sscha_force(x_ensemble, x_c, phi, m = 1)
    u_ensemble = x_ensemble - x_c
    lanc.rho = np.ones(N_RANDOM, dtype = np.double)
    lanc.N = N_RANDOM
    lanc.N_eff = N_RANDOM
    lanc.X = u_ensemble
    lanc.Y = f
    lanc.psi = np.zeros(3, dtype = np.double)

    # Set the temperature
    lanc.T = 0

    # Prepare the 1D phonon green function
    lanc.prepare_perturbation(np.ones(1, dtype = np.double), masses_exp = 0)

    # Fake the algorithm with a false initialization
    lanc.initialized = True

     # Prepare the L as a linear operator (Prepare the possibility to transpose the matrix)
    def L_transp(psi):
        return lanc.apply_full_L(psi, transpose= True, force_FT = True)
    def L_direct(psi):
        return lanc.apply_full_L(psi, transpose= False, force_FT = True)
    lanc.L_linop = scipy.sparse.linalg.LinearOperator(shape = (len(lanc.psi), len(lanc.psi)), matvec = L_direct, rmatvec = L_transp, dtype = np.double)


    # Compute the Lanczos algorithm (biconjugate)
    f1, c1 = get_pert_force_potential(x_c, phi, 0, -1)


    print("Standard perturbation:")
    print("<f> pert = {}".format(f1))
    print("<d2v/dr2> pert = {}".format(c1))
    print()

    # Compare the application on psi
    lanc.psi = np.zeros(3, dtype = np.double)
    lanc.psi[0] = 1 

    L_anal_psi = apply_L(x_c, phi, lanc.psi)
    lanc.apply_full_L(force_FT= True) 


    print()
    print("TEST L APPLICATION")
    print("------------------")
    print()
    print("Application exact:")
    print(L_anal_psi)
    print("Application Direct:")
    print(lanc.psi)
    print()
    print("Test d3 and d4:")
    d3, d4 = get_d3_d4(x_c, phi)
    d3_num, d4_num = get_numerical_d3_d4(x_c, phi)
    print("D3 anal: {} | numeric: {}".format(d3, d3_num))
    print("D4 anal: {} | numeric: {}".format(d4, d4_num))

    print()
    print()
    print("TRANSPOSED")
    # Compare the application on psi
    lanc.psi = np.zeros(3, dtype = np.double)
    lanc.psi[0] = 1 

    L_anal_psi = apply_L(x_c, phi, lanc.psi, transpose = True)
    lanc.apply_full_L(force_FT= True, transpose = True) 

    print("Application exact:")
    print(L_anal_psi)
    print("Application Direct:")
    print(lanc.psi)
    print()


    lanc.run_FT(10)



    # Get the green function
    w_array = np.linspace(0, 10, 1000)
    gf = lanc.get_green_function_continued_fraction(w_array, smearing = 0.1, use_terminator = False)

    if plot:
        plt.plot(w_array, -np.imag(gf))
        ax = plt.gca()
        ax.vlines( sscha_w(1, phi), 0, np.max(-np.imag(gf)), ls = "--", color = "k")
        plt.tight_layout()
        plt.show()


if __name__ == "__main__":
    test_lanczos_1d(plot = True)