import numpy as np
from qiskit import *
from qiskit.extensions import *
from qiskit.compiler import transpile
import random
from numpy import array, matrix
from scipy.io import mmread, mmwrite

BO_fourier = mmread("FTMat.mtx")
#np.set_printoptions(threshold=sys.maxsize,linewidth=400)

q = QuantumRegister(8, 'q')

print(BO_fourier.shape)
print(BO_fourier.round(4))
print("==================")

qc = QuantumCircuit(q)
gateU = UnitaryGate(BO_fourier)
qc.append(gateU, q )

print(qc)

bgates=['rx','ry','rz','cx','ccx','cccx']
#bgates=['u3','cx','ccx']

result = transpile(qc, basis_gates=bgates, optimization_level=3)
#print(result)

print(dict(result.count_ops()))

#print(result.qasm())
