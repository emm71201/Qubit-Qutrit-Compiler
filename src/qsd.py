import numpy
import scipy
from utils import *
#from cosine_sine_decomposition import qubit_cosine_sine_decomposition, qutrit_cosine_sine_decomposition


def moduluz(z):
    return numpy.sqrt(z.real**2 + z.imag**2)

def chop_arr(arr, tol=1e-15):
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            if moduluz(arr[i, j]) < tol:
                arr[i,j] = 0
    return arr


def qubit_qsd(u1, u2):
    """ return a quantum shannon decomposition of a block diagonal unitary
    u1 and u2 are the blocks
     """
    dsq, V = scipy.linalg.eig(u1 @ u2.conj().T)
    Dsq = numpy.diag(dsq)
    D = scipy.linalg.sqrtm(Dsq)
    W = D @ V.conj().T @ u2

    if numpy.allclose(V @ D @ W, u1) and numpy.allclose(V @ D.conj().T @ W, u2):
        return V, D, W

    print("The decomposition failed")
    print("Exiting ...")
    return

def qutrit_qsd(u1, u2, u3):
    """ return a quantum shannon decomposition of a block diagonal unitary u1, u2 and u3 """
    r = u1.shape[0]
    one = numpy.eye(r)
    V,D,W = qubit_qsd(u1, u2)
    Vp,Dp, Wp = qubit_qsd(V,u3)
    Vpp,Dpp,Wpp = qubit_qsd(W,numpy.eye(u1.shape[0]))

    M = direct_sum(Vp,Vp,Vp) @ direct_sum(Dp, Dp, Dp.conj().T) @ direct_sum(Wp,Wp,Wp) @ \
        direct_sum(D,D.conj().T,one) @ direct_sum(Vpp,Vpp,Vpp) @ direct_sum(Dpp,Dpp,Dpp.conj().T) @ \
        direct_sum(Wpp,Wpp,Wpp)

    if numpy.allclose(M, direct_sum(u1,u2,u3)):
        return (Vp, Dp, Wp), D, (Vpp, Dpp, Wpp)

    print("The decomposition failed")
    print("Exiting ...")
    return






# if __name__ == "__main__":
#     dim = 4
#     u1 = scipy.stats.unitary_group.rvs(dim)
#     u2 = scipy.stats.unitary_group.rvs(dim)
#     u3 = scipy.stats.unitary_group.rvs(dim)

    # matrix = np.genfromtxt("matrices/s108fft-108dim.csv", delimiter=",", dtype=numpy.complex_)
    # (u1, x1, x2), csu2, (y1, y2), cs, (v1, g1, g2), csv2, (h1, h2) = qutrit_cosine_sine_decomposition(matrix)
    #
    # (Vp, Dp, Wp), D, (Vpp, Dpp, Wpp) = qutrit_qsd(u1,x1,x2)
    #
    # print(numpy.diag(Dpp)/numpy.pi)

    # (v1,v2), cs, (w1, w2) = scipy.linalg.cossin(u1, p = 2, q = 2, separate=True)
    #
    # print(v1)





