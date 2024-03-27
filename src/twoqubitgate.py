import qiskit
from qiskit.compiler import transpile
from gates import *
from quantum_circuit import *
import scipy
import numpy as np

def decompose(unitary, level, register):

    # det = np.linalg.det(unitary)
    # phi = np.arctan2(det.imag, det.real)    
    
    assert register[-2:] == [2,2]

    final_gates = []

    qskt_qc = qiskit.QuantumCircuit(2)
    qskt_qc.unitary(unitary, [0,1])
    result = transpile(qskt_qc, basis_gates=['rz', 'ry', 'rx', 'cx'],optimization_level=3)
    
    #print(result)


    for index, qiskit_gate in enumerate(result):
        qgate, qubits = qiskit_gate[0], qiskit_gate[1]
        if qgate.name == "cx":
            control, target = level + qubits[0]._index, level + qubits[1]._index

            Xgate = SingleQubitGate(target, pauliX, "X")

            CXgate = ControlledGate(control, 1, Xgate)
            final_gates.append(CXgate)
        
        else:
            angle = qgate.params[0]
            #print(angle/np.pi)
            target = level + qubits[0]._index
            
            if qgate.name == "rx":
                rgate = QubitXRotation(target, angle)
            elif qgate.name == "ry":
                rgate = QubitYRotation(target, angle)
            elif qgate.name == "rz":
                rgate = QubitZRotation(target, angle)
            
            final_gates.append(rgate)
    
    return QuantumCircuit(register, final_gates)
            


    


if __name__ == '__main__':
    dim = 4
    matrix = scipy.stats.unitary_group.rvs(dim)
    #matrix = np.array([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]])
    register = [2,2]
    level = 0

    qc = decompose(matrix, level, register)


#     qc.draw(outputfile_name="two.tex")


    #print(matrix)

    # qiskit_qc = qiskit.QuantumCircuit(2)
    # qiskit_qc.unitary(matrix, [0,1])
    # result = transpile(qiskit_qc, basis_gates=['rz', 'ry', 'rx', 'cx'],\
    #  optimization_level=3)
    # print(result)
    # for index, g in enumerate(result):
    #     for item in g:
    #         # if isinstance(item, list):
    #         #     pass
    #         # else:
    #         #     print(item.name, item.params)
    #         if isinstance(item, qiskit.circuit.Instruction):
    #             print(item, index)

    #     print("")
    # for index, g in enumerate(result):
    #     print()
    #     qubits = g.qubits
    #     print(qubits)