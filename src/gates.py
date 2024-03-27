import numpy
import scipy

# qubit setup
Z = numpy.array([[1, 0], [0, -1]])
pauliX = numpy.array([[0, 1], [1, 0]])
pauliY = numpy.array([[0, -1j], [1j, 0]])
one = numpy.eye(2)
# qutrit setup
Z01 = numpy.array([[1,0,0], [0,-1,0], [0,0,0]])
Z02 = numpy.array([[1,0,0], [0,0,0], [0,0,-1]])
X01 = numpy.array([[0,1,0], [1,0,0], [0,0,1]])
X02 = numpy.array([[0,0,1], [0,1,0], [1,0,0]])
X12 = numpy.array([[1,0,0], [0,0,1], [0,1,0]])

Y01 = numpy.array([[0,-1j,0], [1j,0,0], [0,0,1]])
Y02 = numpy.array([[0,0,-1j], [0,1,0], [1j,0,0]])
Y12 = numpy.array([[1,0,0], [0,0,-1j], [0,1j,0]])

ONE = numpy.eye(3)

class SingleQubitGate:
    def __init__(self, target=None, op=None, name=None):
        self.target = target
        self.op = op
        if name is None:
            self.name = 'single qubit gate'
        else:
            self.name = name

    def set_target(self, target):
        self.target = target
    def set_op(self, op):
        assert op.shape[0] == 2 and op.shape[1] == 2
        self.op = op

    def __str__(self):
        return f"Single Qubit Gate: Name = {self.name} Target = {self.target}"

    def to_latex(self):
        return r"\gate{" + f"{self.name}" + "}"

class SingleQutritGate:
    def __init__(self, target=None, op=None, name=None):
        self.target = target
        self.op = op
        if name is None:
            self.name = 'single qubit gate'
        else:
            self.name = name

    def set_target(self, target):
        self.target = target
    def set_op(self, op):
        assert op.shape[0] == 3 and op.shape[1] == 3
        self.op = op
    def __str__(self):
        return f"Single Qutrit Gate: Target = {self.target} name: {self.name}"
    def to_latex(self):
        return r"\gate{" + f"{self.name}" + "}"

class QubitZRotation(SingleQubitGate):

    def __init__(self, target, angle):
        op = scipy.linalg.expm(-1j * angle * Z /2)
        super().__init__(target, op)
        self.angle = angle

    def set_op(self, angle):
        self.op = scipy.linalg.expm(-1j * angle * Z /2)
        self.angle = angle

    def __str__(self):
        return f"Rz: Target = {self.target} Angle = {self.angle}"

    def to_latex(self):
        return  r"\gate{R_z" + f"({self.angle:.4f})" + "}"



class QubitYRotation(SingleQubitGate):

    def __init__(self, target, angle):
        op = scipy.linalg.expm(-1j * angle * pauliY /2)
        super().__init__(target, op)
        self.angle = angle

    def set_op(self, angle):
        self.op = scipy.linalg.expm(-1j * angle * pauliY /2)
        self.angle = angle

    def __str__(self):
        return f"Ry: Target = {self.target} Angle = {self.angle}"

    def to_latex(self):
        return  r"\gate{R_y" + f"({self.angle:.4f})" + "}"


class QubitXRotation(SingleQubitGate):

    def __init__(self, target, angle):
        op = scipy.linalg.expm(-1j * angle * pauliX /2)
        super().__init__(target, op)
        self.angle = angle

    def set_op(self, angle):
        self.op = scipy.linalg.expm(-1j * angle * pauliX /2)
        self.angle = angle

    def __str__(self):
        return f"Ry: Target = {self.target} Angle = {self.angle}"

    def to_latex(self):
        return  r"\gate{R_x" + f"({self.angle:.4f})" + "}"


class QutritZRotation(SingleQutritGate):

    def __init__(self, target, angle, type):
        if type == "01":
            sigma = Z01
        elif type == "02":
            sigma = Z02
        else:
            print(f"Invalid Qutrit Rotation Type {type}")
            print("Exiting ...")
            return
        op = scipy.linalg.expm(-1j * angle * sigma/2)
        super().__init__(target, op)
        self.angle = angle
        self.type = type

    def set_op(self, angle):
        if self.type == "01":
            sigma = Z01
        elif self.type == "02":
            sigma = Z02
        else:
            print(f"Invalid Qutrit Rotation Type {type} ")
            print("Exiting ...")
            return
        self.op = scipy.linalg.expm(-1j * angle * sigma/2)
        self.angle = angle

    def set_target(self, target):
        self.target = target

    def __str__(self):
        return f"Qutrit {self.type} Z Rotation: Target = {self.target} Angle = {self.angle}"

    def to_latex(self):
        return r"\gate{R" + "^{" + f"{self.type}" + "}" + "_z" + f"({self.angle:.4f})" + "}"


class QutritYRotation(SingleQutritGate):

    def __init__(self, target, angle, type):
        if type == "01":
            sigma = Y01
        elif type == "02":
            sigma = Y02
        elif type == "12":
            sigma = Y12
        else:
            print(f"Invalid Qutrit Rotation Type {type}")
            print("Exiting ...")
            return
        op = scipy.linalg.expm(-1j * angle * sigma/2)
        super().__init__(target, op)
        self.angle = angle
        self.type = type

    def set_op(self, angle):
        if type == "01":
            sigma = Y01
        elif type == "02":
            sigma = Y02
        elif type == "12":
            sigma = Y12
        else:
            print(f"Invalid Qutrit Rotation Type {type}")
            print("Exiting ...")
            return
        self.op = scipy.linalg.expm(-1j * angle * sigma/2)
        self.angle = angle

    def set_target(self, target):
        self.target = target

    def __str__(self):
        return f"Qutrit {self.type} Y Rotation: Target = {self.target} Angle = {self.angle}"

    def to_latex(self):
        return r"\gate{R" + "^{" + f"{self.type}" + "}" + "_y" + f"({self.angle:.4f})" + "}"

class QutritXRotation(SingleQutritGate):

    def __init__(self, target, angle, type):
        if type == "01":
            sigma = X01
        elif type == "02":
            sigma = X02
        elif type == "12":
            sigma = X12
        else:
            print(f"Invalid Qutrit Rotation Type {type}")
            print("Exiting ...")
            return
        op = scipy.linalg.expm(-1j * angle * sigma/2)
        super().__init__(target, op)
        self.angle = angle
        self.type = type

    def set_op(self, angle):
        if type == "01":
            sigma = X01
        elif type == "02":
            sigma = X02
        elif type == "12":
            sigma = X12
        else:
            print(f"Invalid Qutrit Rotation Type {type}")
            print("Exiting ...")
            return
        self.op = scipy.linalg.expm(-1j * angle * sigma/2)
        self.angle = angle

    def set_target(self, target):
        self.target = target

    def __str__(self):
        return f"Qutrit {self.type} X Rotation: Target = {self.target} Angle = {self.angle}"

    def to_latex(self):
        return r"\gate{R" + "^{" + f"{self.type}" + "}" + "_x" + f"({self.angle:.4f})" + "}"


class ControlledGate:

    def __init__(self, control, state, gate):
        self.control = control
        self.state = state
        self.gate = gate

    def __str__(self):
        return f"Controlled Gate: Control = {self.control} State = {self.state} " + self.gate.__str__()

if __name__ == "__main__":
    angle = numpy.pi
    type = "02"
    gate = QutritZRotation(0, angle, type)
    print(gate.to_latex())
    print(gate.op)



