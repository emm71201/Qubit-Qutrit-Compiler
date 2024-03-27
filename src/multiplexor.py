import numpy
import itertools
from gates import *

def get_controls(register, level):

    def pattern(qudit):
        return itertools.chain(["-"], [str(ch) for ch in range(qudit)][1:])

    return itertools.product(*map(pattern, register[level+1:]))

def perform_contrlx01(word, cntrls):
    """ compute the operation of control X01 on the top qubit
    word -- > string eg. 011
    cntrls --> string eg. 11
    return new word 011 in this case
     """
    cntrls = "-" + "".join([ch for ch in cntrls])
    assert len(word) == len(cntrls)
    word = [ch for ch in word]
    i = len(word) - 1
    while i > 0:
        if cntrls[i] == "-":
            pass
        elif cntrls[i] == "1" and word[i] == "1":
            word[0] = str((int(word[0]) + 1) % 2)

        elif cntrls[i] == "2" and word[i] == "2":
            word[0] = str((int(word[0]) + 1) % 2)

        i -= 1

    return "".join(word)

def generating_matrix(cntrls, word):
    word = "".join([str(ch) for ch in word])
    val = int(perform_contrlx01(word, cntrls) != word)
    return (-1)**val

def get_thetaMatrix(cntrls_list, words, dim):
    thetaMatrix = numpy.zeros(shape=(dim, dim))
    words = [next(words) for i in range(dim)]

    assert len(cntrls_list) == len(words)

    for i in range(dim):
        word = words[i]
        for j in range(dim):
            cntrls = cntrls_list[j]
            thetaMatrix[i,j] = generating_matrix(cntrls, word)

    return thetaMatrix

# Decompose multiplexor of type D + D^dagger + 1
def decompose_multiplexor01(D, level, register):

    words = itertools.product(*map(range, register[level:]))
    cntrls_list = [item for item in get_controls(register, level)]


    thetaMatrix = get_thetaMatrix(cntrls_list, words, D.shape[0])
    thetas = numpy.linalg.solve(thetaMatrix, 2*D)

    mygates = []
    current = 0


    for cntrls in cntrls_list:
        if register[level] == 2:
            gate = SingleQubitGate(level, pauliX, "CX")
            rgate = QubitZRotation(level, thetas[current])
        elif register[level] == 3:
            gate = SingleQutritGate(level, X01, "CX01")
            rgate = QutritZRotation(level, thetas[current], "01")

        cxqueue = []

        for ii in range(len(cntrls)):
            if cntrls[ii] == "-":
                pass
            else:
                state = int(cntrls[ii])
                control = ii + level + 1


                cgate = ControlledGate(control=control, state=state, gate=gate)
                mygates.append(cgate)
                cxqueue.append(cgate)



        mygates.append(rgate)

        for cgate in cxqueue[::-1]:
            mygates.append(cgate)

        current += 1

    for mygate in mygates:
        print(mygate)

    return cntrls_list, thetas




if __name__ == "__main__":
    register = [3,3,3,2,2]
    register = [3, 2, 2]
    #register = [2,2,2]
    D = numpy.random.rand(4)
    print(D)
    # for index, item in enumerate(itertools.product(*map(range, register[:]))):
    #     print(item)

    #get_controls(register, 0)

    cntrls, thetas = decompose_multiplexor01(D, 0, register)


