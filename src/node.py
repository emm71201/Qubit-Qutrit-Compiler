import numpy
import csd
import qsd
from utils import *

class Node:

    def __init__(self, op, level, type):
        """ unitary matrix class
        op: 2D numpy array
        level: int
        type: int
         """
        self.op = op
        self.level = level
        self.type = type
        self.dim = len(op)
        self.children = {}
        self.parent = None

    def get_dim(self):
        if self.type  == 0:
            return 2 * len(self.op)

    def __str__(self):
        return f"Tree Node: Level = {self.level} Dimension = {self.dim}"

    def is_identity(self):
        return numpy.allclose(self.op, numpy.eye(self.dim))

    def get_unitary(self, top_qudit="qubit"):
        if self.type is None:
            return self.op

        if self.type == -1:
            if top_qudit == "qubit":
                return direct_sum(self.op, self.op)
            if top_qudit == "qutrit":
                return direct_sum(self.op, self.op, self.op)

        if self.type == 0:
            C,S = numpy.diag(numpy.cos(self.op)), numpy.diag(numpy.sin(self.op))
            return direct_sum(C,C) + off_diagonal_direct_sum(S, -S)
        if self.type == 1:
            return direct_sum(self.op, self.op.conj().T)

        if self.type == 2:
            one = numpy.eye(len(self.op))
            tmp = Node(self.op, self.level, 0)
            return direct_sum(one, tmp.get_unitary())

        if self.type == 3:
            r = len(self.op)
            one = numpy.eye(r)
            zero = numpy.zeros(shape=(r,r))
            C, S = numpy.diag(numpy.cos(self.op)), numpy.diag(numpy.sin(self.op))
            return direct_sum(C, one, C) + off_diagonal_direct_sum(S, zero, -S)

        if self.type == 4:
            return direct_sum(self.op,self.op,self.op.conj())
        if self.type == 5:
            one = numpy.eye(len(self.op))
            return direct_sum(self.op,self.op.conj(), one)

    def multiply(self, *others, top_qudit):
        res = self.op
        for other in others:
            res = res @ other.get_unitary(top_qudit=top_qudit)
        return res


    def verify_equality(self, *others, top_qudit):
        res = numpy.eye(self.dim)
        for other in others:
            #print(other.type)
            try:
                res = res @ other.get_unitary(top_qudit=top_qudit)
            except:
                print(other.type)

        return numpy.allclose(self.op, res)

    def decompose(self, type="qubit", verify=True):
        """ decompose a unitary matrix
        type: qubit decomposition or qutrit decomposition
        """
        # qubit type
        if type == "qubit":
            (u1, u2), cs, (v1, v2) = csd.qubit_cosine_sine_decomposition(self.op, separate=True)
            V1, D1, W1 = qsd.qubit_qsd(u1,u2)
            V2, D2, W2 = qsd.qubit_qsd(v1, v2)

            D1 = Node(D1, self.level, 1)
            cs = Node(cs, self.level, 0)
            D2 = Node(D2, self.level, 1)

            V1 = Node(V1, self.level+1, -1)
            W1 = Node(W1, self.level+1, -1)
            V2 = Node(V2, self.level+1, -1)
            W2 = Node(W2, self.level+1, -1)

            if verify:
                if self.verify_equality(V1,D1,W1,cs,V2,D2,W2, top_qudit="qubit"):
                    for indx, new_node in enumerate([V1,D1,W1,cs,V2,D2,W2]):
                        self.children[indx] = new_node
                        new_node.parent = self
                    return True
                else:
                    print("Decomposition Failed")
                    return False

        # qutrit type
        if type == "qutrit":
            one = numpy.eye(self.dim//3)
            (u1, x1, x2), csu2, (y1, y2), cs, (v1, g1, g2), csv2, (h1, h2) = csd.qutrit_cosine_sine_decomposition(self.op, separate=True)

            (V1p, D1p, W1p), D1, (V1pp, D1pp, W1pp) = qsd.qutrit_qsd(u1,x1,x2)
            (V2p, D2p, W2p), D2, (V2pp, D2pp, W2pp) = qsd.qutrit_qsd(one,y1,y2)
            (V3p, D3p, W3p), D3, (V3pp, D3pp, W3pp) = qsd.qutrit_qsd(v1, g1, g2)
            (V4p, D4p, W4p), D4, (V4pp, D4pp, W4pp) = qsd.qutrit_qsd(one, h1, h2)

            cs = Node(cs, self.level, 3)
            csu2, csv2 = Node(csu2, self.level, 2), Node(csv2, self.level, 2)
            D1p, D1, D1pp = Node(D1p, self.level, 4), Node(D1, self.level, 5), Node(D1pp, self.level, 4)
            D2p, D2, D2pp = Node(D2p, self.level, 4), Node(D2, self.level, 5), Node(D2pp, self.level, 4)
            D3p, D3, D3pp = Node(D3p, self.level, 4), Node(D3, self.level, 5), Node(D3pp, self.level, 4)
            D4p, D4, D4pp = Node(D4p, self.level, 4), Node(D4, self.level, 5), Node(D4pp, self.level, 4)

            V1p, W1p, V1pp, W1pp = Node(V1p, self.level + 1, -1), Node(W1p, self.level + 1, -1), \
                Node(V1pp, self.level + 1, -1), Node(W1pp, self.level + 1, -1)
            V2p, W2p, V2pp, W2pp = Node(V2p, self.level + 1, -1), Node(W2p, self.level + 1, -1), \
                Node(V2pp, self.level + 1, -1), Node(W2pp, self.level + 1, -1)
            V3p, W3p, V3pp, W3pp = Node(V3p, self.level + 1, -1), Node(W3p, self.level + 1, -1), \
                Node(V3pp, self.level + 1, -1), Node(W3pp, self.level + 1, -1)
            V4p, W4p, V4pp, W4pp = Node(V4p, self.level + 1, -1), Node(W4p, self.level + 1, -1), \
                Node(V4pp, self.level + 1, -1), Node(W4pp, self.level + 1, -1)

            if verify:
                if self.verify_equality(V1p,D1p,W1p,D1,V1pp,D1pp,W1pp, csu2, V2p,D2p,W2p,D2,V2pp,D2pp,W2pp, cs, V3p,D3p,W3p,D3,V3pp,D3pp,W3pp,csv2, V4p,D4p,W4p,D4,V4pp,D4pp,W4pp, top_qudit="qutrit"):
                    for indx, new_node in enumerate([V1p,D1p,W1p,D1,V1pp,D1pp,W1pp, csu2, V2p,D2p,W2p,D2,V2pp,D2pp,W2pp, cs, V3p,D3p,W3p,D3,V3pp,D3pp,W3pp,\
                                                     csv2, V4p,D4p,W4p,D4,V4pp,D4pp,W4pp]):
                        self.children[indx] = new_node
                        new_node.parent = self
                else:
                    print("Decomposition Failed")

        return













