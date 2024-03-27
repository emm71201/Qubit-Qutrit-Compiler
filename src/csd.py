import scipy
import numpy as np
import numpy
from utils import *

# def direct_sum(*matrices):
#     # Calculate the total size of the resulting matrix
#     total_rows = sum(matrix.shape[0] for matrix in matrices)
#     total_cols = sum(matrix.shape[1] for matrix in matrices)
#
#     # Initialize the result matrix with zeros
#     result = np.zeros((total_rows, total_cols))
#
#     # Place each matrix in its block-diagonal position
#     row_start, col_start = 0, 0
#     for matrix in matrices:
#         rows, cols = matrix.shape
#         result[row_start:row_start+rows, col_start:col_start+cols] = matrix
#         row_start += rows
#         col_start += cols
#
#     return result

def qubit_cosine_sine_decomposition(matrix, separate=True):
    assert matrix.shape[0] % 2 == 0
    r = matrix.shape[0] // 2

    (u1, u2), cs, (v1, v2) = scipy.linalg.cossin(matrix, p=r, q=r, separate=separate)

    return (u1, u2), cs, (v1, v2)

def qutrit_cosine_sine_decomposition(matrix, separate=True, verify=True):

    assert matrix.shape[0] % 3 == 0

    r = matrix.shape[0]//3 # r = 3^(n-1)
    # first decomposition
    (u1,u2), cs, (v1,v2) = scipy.linalg.cossin(matrix, p=r, q=r, separate=separate)

    # second decomposition
    # dimension of u2 and v2 = 2*r
    rp = matrix.shape[0]//3 # now, we set rp = 3^(n-1) again

    # decompose u2 and v2
    (x1, x2), csu2, (y1,y2) = scipy.linalg.cossin(u2, p=rp, q=rp, separate=separate)
    (g1, g2), csv2, (h1,h2) = scipy.linalg.cossin(v2, p=rp, q=rp, separate=separate)

    if not verify:
        return (u1, x1, x2), csu2, (y1, y2), cs, (v1, g1, g2), csv2, (h1, h2)

    if not verify_qutrit(u1, x1, x2, csu2, y1,y2, cs, v1,g1,g2, csv2, h1,h2, r, matrix):
        print("Qutrit CSD Decomposition Failed")
        return

    return (u1, x1, x2), csu2, (y1, y2), cs, (v1,g1,g2), csv2, (h1,h2)

def verify_qubit(u1, u2, cs, v1, v2, op):

    U = direct_sum(u1, u2)
    CS = direct_sum(np.diag(np.cos(cs)), np.diag(np.cos(cs))) + \
         off_diagonal_direct_sum(np.diag(np.sin(cs)), -np.diag(np.sin(cs)))
    V = direct_sum(v1, v2)

    return np.allclose(U @ CS @ V, op)

def verify_qutrit(u1, x1, x2, csu2, y1,y2, cs, v1,g1,g2, csv2, h1,h2, r, op):
    id_matrix = np.eye(r)
    X = direct_sum(u1, x1, x2)
    csp = direct_sum(np.diag(np.cos(csu2)), np.diag(np.cos(csu2))) + \
          off_diagonal_direct_sum(np.diag(np.sin(csu2)), -np.diag(np.sin(csu2)))
    CSP = direct_sum(id_matrix, csp)
    Y = direct_sum(id_matrix, y1, y2)

    CS = direct_sum(np.diag(np.cos(cs)), np.zeros(shape=(r, r)), np.diag(np.cos(cs))) + \
         direct_sum(np.zeros(shape=(r, r)), id_matrix, np.zeros(shape=(r, r))) + \
         off_diagonal_direct_sum(np.diag(np.sin(cs)), np.zeros(shape=(r, r)), -np.diag(np.sin(cs)))

    V = direct_sum(v1, g1, g2)

    cspp = direct_sum(np.diag(np.cos(csv2)), np.diag(np.cos(csv2))) + \
           off_diagonal_direct_sum(np.diag(np.sin(csv2)), -np.diag(np.sin(csv2)))
    CSPP = direct_sum(id_matrix, cspp)
    H = direct_sum(id_matrix, h1, h2)

    return np.allclose(X @ CSP @ Y @ CS @ V @ CSPP @ H, op)


def decompose(matrix, r, separate=True):
    return scipy.linalg.cossin(matrix, p=r, q=r, separate=separate)



# if __name__ == "__main__":
#
#     matrix = np.genfromtxt("matrices/s108fft-108dim.csv", delimiter=",", dtype=numpy.complex_)
    #print(matrix.shape)
    # print(matrix)
    #dim = 27
    #matrix = scipy.stats.unitary_group.rvs(dim)
    #u, cs, v = decompose(matrix, matrix.shape[0]//3, separate=True)

    # for phi in cs:
    #     print(phi)



    #sympy.pprint(sympy.Matrix(u))
    # matrix = u[1]
    # print(matrix.shape)
    # u, cs, v = decompose(matrix, matrix.shape[0] // 2)
    # sympy.pprint(sympy.Matrix(cs))

    #(u1, x1, x2), csu2, (y1, y2), cs, (v1,g1,g2), csv2, (h1,h2) = qutrit_cosine_sine_decomposition(matrix)


