import numpy as np 

I = np.eye(3)
def FromGradToStress(uc, grad):
    Omega = np.linalg.det(uc)
    # du_deps = np.zeros( (3,3,3,3))

    # for i in range(3):
    #     for j in range(3):
    #         for x in range(3):
    #             for y in range(3):
    #                 du_deps[i,j,x,y] = 0.5 * (uc[i,y] * I[x,j] + uc[i,x] * I[y,j])
    
    # dF_deps = np.einsum("ab, abcd", grad, du_deps)

    # return - dF_deps / Omega

    return 2 * np.transpose(uc).dot(grad)

def FromStressToGrad(uc, stress):


    volume = np.linalg.det(uc)
    uc_inv = np.linalg.inv(uc)
    grad_mat = - volume * np.transpose(uc_inv).dot(stress)
    return grad_mat


UC = np.random.uniform( size=(3,3))
GRAD_UC = np.random.uniform( size=(3,3))
STRESS = FromGradToStress(UC, GRAD_UC)


print "UC:"
print UC
print "GRAD:"
print GRAD_UC
print "STRESS:"
print STRESS
print "BACK AGAIN:"
print FromStressToGrad(UC, STRESS)
