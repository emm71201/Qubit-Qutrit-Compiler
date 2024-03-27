import numpy
import scipy
import itertools
from gates import *
from utils import modulus
from utils import is_diagonal_matrix
from utils import chop
from quantum_circuit import QuantumCircuit

def f01(word, j):
    """
    word string
    int j: qudit where a Rz^01 rotation acts on
    return int 1, -1, or 0
    """
    if int(word[j]) == 0:
        return 1
    elif int(word[j]) == 1:
        return -1
    elif int(word[j]) == 2:
        return 0

    return

def f02(word, j):
    if int(word[j]) == 0:
        return 1
    elif int(word[j]) == 1:
        return 0
    elif int(word[j]) == 2:
        return -1

    return

def cx01(word, cntrl,target, cntrl_state):

    word = [str(ch) for ch in word]
    cntrl_state = str(cntrl_state)


    if word[cntrl] == cntrl_state:
        if word[target] == "0":
            word[target] = "1"
        elif word[target] == "1":
            word[target] = "0"

    return "".join(word)

def cx02(word, cntrl,target, cntrl_state):

    word = [str(ch) for ch in word]
    cntrl_state = str(cntrl_state)


    if word[cntrl] == cntrl_state:
        if word[target] == "0":
            word[target] = "2"
        elif word[target] == "2":
            word[target] = "0"

    return "".join(word)

def mergef(word, j, controls, rot_type):
    if j is None:
        return 1

    for cntr_item  in controls:
        cntrl, target, cntrl_state = cntr_item
        if rot_type == "01":
            word = cx01(word, cntrl, target, cntrl_state)
        if rot_type == "02":
            word = cx02(word, cntrl, target, cntrl_state)


    if rot_type == "01":
        return f01(word, j)
    if rot_type == "02":
        return f02(word, j)

    return

def get_angles(reduced_register, D):
    dim = numpy.prod(numpy.array(reduced_register))
    M = numpy.empty(shape=(dim,dim))
    if D.ndim == 2:
        if not is_diagonal_matrix(D):
            return
        D = numpy.diag(D)

    phases = -1j * numpy.log(D)

    curr_rot = 0
    allcontrols = {}
    rot_types = {}
    for cntrl_pattern in itertools.product(*map(range, reduced_register)):
        target = None
        rot_type = None
        controls = []
        for index, ch in enumerate(cntrl_pattern):
            if ch == 1 or ch == 2:
                if target is None:
                    target = index
                    if ch == 1: rot_type = "01"
                    if ch == 2: rot_type = "02"
                elif index > target:
                    cand_controls = (index, target, str(ch))
                    if cand_controls not in controls:
                        controls.append(cand_controls)
                    else:
                        controls.remove(cand_controls)
        M[:, curr_rot] = numpy.array([mergef(word, target, controls, rot_type) \
                           for word in itertools.product(*map(range, reduced_register))])
        allcontrols[curr_rot] = controls
        rot_types[curr_rot] = (rot_type, target)
        curr_rot += 1

    angles = scipy.linalg.solve(M, phases)

    return map(chop, angles), allcontrols, rot_types

def make_gates(D, register, level, tol=1e-12):

    reduced_register = register[level:]
    angles, allcontrols, rot_types = get_angles(reduced_register, D)
    allgates = []

    curr_rot = 0
    for angle in angles:
        if numpy.abs(angle) > tol:
            rot_type, target = rot_types[curr_rot]

            if target is not None:

                # build the single qudit rotation gate
                target += level
                if register[target] == 2:
                    rotation_gate = QubitZRotation(target, angle)
                    xgate = SingleQubitGate(target, pauliX, "X")
                if register[target] == 3:
                    rotation_gate = QutritZRotation(target, angle, rot_type)
                    if rot_type == "01": xgate = SingleQutritGate(target, X01, "X^{01}")
                    if rot_type == "02": xgate = SingleQutritGate(target, X02, "X^{02}")

                controls = allcontrols[curr_rot]
                tmp_control_gates = []
                for cntrl_item in controls:
                    cntrl, _, cntrl_state = cntrl_item
                    cgate = ControlledGate(cntrl+level, int(cntrl_state), xgate)

                    allgates.append(cgate)
                    tmp_control_gates.append(cgate)

                allgates.append(rotation_gate)
                for cgate in tmp_control_gates[::-1]:
                    allgates.append(cgate)

        curr_rot += 1

    return QuantumCircuit(register, allgates)
    #return allgates




if __name__ == "__main__":
    word = "1012"
    cntrl_qudit = 3
    target_type = 3
    register = numpy.array([3,3, 3])
    level = 0
    dim = numpy.prod(register[level:])

    U = numpy.exp(-1j * numpy.random.rand(dim))

    allgates = make_gates(U, register, level, tol=1e-12)
    allgates.draw()
    # for gate in all_gates:
    #     print(gate)







