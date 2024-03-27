import numpy
import scipy
from utils import *

rtol = 1e-7
atol = 1e-10
def arccos(value):
    if numpy.isclose(value, 1, rtol=rtol, atol=atol):
        return 0
    if numpy.isclose(value, -1):
        return numpy.pi
    if abs(value) > 1:
        raise ValueError("Invalid value for arccos(%f)" % value)
    return numpy.arccos(value)

def arcsin(value):
    if numpy.isclose(value, 1, rtol=rtol, atol=atol):
        return numpy.pi/2
    if numpy.isclose(value, -1):
        return -numpy.pi/2
    if abs(value) > 1:
        raise ValueError("Invalid value for arccos(%f)" % value)
    return numpy.arcsin(value)


def rotz(alpha):
    z = numpy.array([[1,0],[0,-1]])
    return scipy.linalg.expm(-1j * alpha * z/2)

def roty(alpha):
    y = numpy.array([[0,-1j], [1j, 0]])
    return scipy.linalg.expm(-1j * alpha * y/2)

def abs(z):
    return numpy.sqrt(z.real**2 + z.imag**2)
def single_qubit_gate_decompose(matrix):
    if matrix.ndim == 1:
        matrix = numpy.diag(matrix)
    assert matrix.ndim == 2

    det = numpy.linalg.det(matrix)
    phi = (1/2) * numpy.arctan2(det.imag, det.real)
    V = numpy.exp(-1j * phi) * matrix
    #V = V.conj().transpose()
    detv = numpy.linalg.det(V)

    if not numpy.isclose(modulus(detv),1, rtol=rtol, atol=atol):
        print(matrix)
        print(numpy.linalg.det(matrix))
        print(modulus(detv))

    assert numpy.isclose(modulus(detv),1, rtol=rtol, atol=atol)

    theta1 = 2 * arccos(abs(V[0,0])) if abs(V[0,0]) >= abs(V[0,1]) else 2 * arcsin(abs(V[0,1]))

    if numpy.cos(theta1/2) != 0:
        val = V[1,1]/numpy.cos(theta1/2)
        tmp_sum = 2 * numpy.arctan2(val.imag, val.real)
    else:
        tmp_sum = 0

    if numpy.sin(theta1/2) != 0:
        val = V[1,0]/numpy.sin(theta1/2)
        tmp_diff = 2 * numpy.arctan2(val.imag, val.real)
    else:
        tmp_diff = 0

    theta0, theta2 = (tmp_sum - tmp_diff)/2, (tmp_sum + tmp_diff)/2

    if not numpy.allclose(numpy.exp(1j * phi) * rotz(theta2) @ roty(theta1) @ rotz(theta0), matrix):
        print("Single Qubit Gate ZYZ decomposition incorrect")
        print("Exiting ...")

        return

    return phi, theta0, theta1, theta2



if __name__ == "__main__":
    dim = 2
    #matrix = scipy.stats.unitary_group.rvs(dim)
    #matrix = numpy.load("matrices/zyz_test.npy")
    matrix = scipy.stats.unitary_group.rvs(dim)
    print("Matrix")
    print(matrix)
    print("#"*30)
    phi, theta0, theta1, theta2 = single_qubit_gate_decompose(matrix)
    print(phi, theta0, theta1, theta2)

