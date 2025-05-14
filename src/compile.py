import numpy as np
import scipy as sp
import gates
from quantum_circuit import QuantumCircuit
import diagonal_unitary
import decompose
import single_qubit_gate_decompose as sqg
import twoqubitgate
from utils import *
import tqdm

pauliZ = gates.Z
pauliX = gates.pauliX
pauliY = gates.pauliY

def compile(matrix, register):

    print("\n\t\ta. Performing the QSD decomposition\n")
    root = decompose.decompose(matrix, register)
    print("\n\t\tb. Assembling the Circuit\n")
    qc = make_circuit(root, register)
    print("\n\t\tc. Counting Gates\n")
    count, controls_dict = qc.count_gates()

    for item, val in count.items():
        if val != 0:
            print(f"\t\t\t{item} : {val}")

    return

def make_circuit(root, register):
    qc = QuantumCircuit(register)

    queue = [root]
    #print("Start making circuit")
    pbar = tqdm.tqdm(total=len(queue), desc=f"Assembling Quantum Circuit Progress")
    while len(queue) > 0:
        last_len_queue = len(queue)

        node = queue.pop(0)

        if node.type == 0:

            qc.insert_gate(gates.QubitXRotation(node.level, -np.pi / 2))
            D = np.array([np.exp(-1j * item) for item in node.op] + [np.exp(1j * item) for item in node.op])
            dpart = diagonal_unitary.make_gates(D, register, node.level)

            for new_gate in dpart.gates:
                qc.insert_gate(new_gate)
            qc.insert_gate(gates.QubitXRotation(node.level, np.pi / 2))

        elif node.type == 1:
            D = np.array([item for item in np.diag(node.op)] + [np.conjugate(item) for item in np.diag(node.op)])
            dpart = diagonal_unitary.make_gates(D, register, node.level)
            for new_gate in dpart.gates:
                qc.insert_gate(new_gate)

        elif node.type == 2:

            qc.insert_gate(gates.QutritXRotation(node.level, -np.pi / 2, "12"))
            D = np.array([1 for item in node.op] + [np.exp(-1j * item) for item in node.op] \
                         + [np.exp(1j * item) for item in node.op])
            dpart = diagonal_unitary.make_gates(D, register, node.level)
            for new_gate in dpart.gates:
                qc.insert_gate(new_gate)
            qc.insert_gate(gates.QutritXRotation(node.level, np.pi / 2, "12"))

        elif node.type == 3:

            qc.insert_gate(gates.QutritXRotation(node.level, -np.pi / 2, "02"))
            D = np.array([1 for item in node.op] + [np.exp(-1j * item) for item in node.op] \
                         + [np.exp(1j * item) for item in node.op])
            dpart = diagonal_unitary.make_gates(D, register, node.level)
            for new_gate in dpart.gates:
                qc.insert_gate(new_gate)
            qc.insert_gate(gates.QutritXRotation(node.level, np.pi / 2, "02"))


        elif node.type == 4:
            D = np.array([item for item in np.diag(node.op)] + [item for item in np.diag(node.op)]\
                + [np.conjugate(item) for item in np.diag(node.op)])
            dpart = diagonal_unitary.make_gates(D, register, node.level)
            for new_gate in dpart.gates:
                qc.insert_gate(new_gate)

        elif node.type == 5:
            D = np.array([item for item in np.diag(node.op)] + [np.conjugate(item) for item in np.diag(node.op)] \
                         + [1 for item in np.diag(node.op)])
            dpart = diagonal_unitary.make_gates(D, register, node.level)
            for new_gate in dpart.gates:
                qc.insert_gate(new_gate)


        # if node is of dimension four, I should throw it to qiskit
        elif node.dim == 4 and node.type == -1:
            two_qbit_circ = twoqubitgate.decompose(node.op, node.level, register)
            for gg in two_qbit_circ.gates:
                qc.insert_gate(gg)

        # if node is of dimension 2, we need to convert to ZYZ rotation
        elif node.dim == 2 and node.type == -1:

            print("This gate should not happen")



            target = node.level
            phi, theta0, theta1, theta2 =  sqg.single_qubit_gate_decompose(node.op)
            if not np.isclose(theta0, 0):
                qc.insert_gate(gates.QubitZRotation(target, theta0))
            if not np.isclose(theta1, 0):
                qc.insert_gate(gates.QubitYRotation(target, theta1))
            if not np.isclose(theta2, 0):
                qc.insert_gate(gates.QubitZRotation(target, theta2))

        elif node.type is None or node.type == -1:
            pass
        else:
            print("This should not happen")
            print(node)
            return

        keys = [key for key in node.children.keys()]
        keys.sort()
        for key in keys[::-1]:
            queue.append(node.children[key])

        #print("Len queue = ", len(queue), "Length QC = ", len(qc))
        if len(queue) > last_len_queue:
            pbar.reset(0)
            pbar.total = len(queue)
            pbar.postfix = f"Length Quantum Circuit: {len(qc)}"
            pbar.update(0)
        else:
            pbar.postfix = f"Length Quantum Circuit: {len(qc)}"
            pbar.update(1)
        #if len(queue) > last_len_queue:


    pbar.close()

    return qc


#if __name__ == "__main__":
    # matrix = np.genfromtxt("matrices/s108fft-108dim.csv", delimiter=",", dtype=np.complex128)
    # matrix = np.array(matrix)
    # register = [3,3,3,2,2]

    # dim = 6
    # register = [3,2]
    # matrix = sp.stats.unitary_group.rvs(dim)

    # root = decompose.decompose(matrix, register)
    # qc = make_circuit(root, register)

    #qc.draw()
    #count = qc.count_gates()
    
    # for item, val in count.items():
    #     if val != 0:
    #         print(f"{item} : {val}")

    # dim = 4
    # matrix = sp.stats.unitary_group.rvs(dim)
    # from qiskit import *
    # circ = QuantumCircuit(2)
    # circ.unitary(matrix, [0,1])
    # #circ.decompose()
    # for gate in circ.decompose():
    #     print(gate)

    # line 37 to line 143

