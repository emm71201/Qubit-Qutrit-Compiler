from gates import *
import sys
class QuantumCircuit:

    def __init__(self, register, gates = [], names = None):
        """ register: list example [2,2,3,3] """
        self.nqubits = 0
        self.nqutrits = 0
        self.register = register
        for item in register:
            if item == 2:
                self.nqubits += 1
            if item == 3:
                self.nqutrits += 1
        self.gates = gates
        if names is None:
            names = []
            tmp1, tmp2 = 0, 0
            for item in register:
                if item == 2:
                    names.append(f"p_{tmp1}")
                    tmp1 += 1
                else:
                    names.append(f"q_{tmp2}")
                    tmp2 += 1
        self.names = names


    def insert_gate(self, gate, pos = -1):
        self.gates.insert(pos,gate)

    def __len__(self):
        return len(self.gates)

    def count_gates(self):
        gate_dict = {"R_x^{01}": 0, "R_y^{01}":0, "R_z^{01}":0,\
                     "R_x^{02}":0, "R_y^{02}":0, "R_z^{02}":0,\
                     "R_x^{12}":0, "R_y^{12}":0,\
                     "R_x":0, "R_y":0, "R_z":0}

        controls = []
        for gate in self.gates:
            if isinstance(gate,QubitXRotation):
                gate_dict["R_x"] += 1
            elif isinstance(gate,QubitYRotation):
                gate_dict["R_y"] += 1
            elif isinstance(gate,QubitZRotation):
                gate_dict["R_z"] += 1

            elif isinstance(gate, QutritXRotation):
                gate_dict["R_x^{" + gate.type + "}"] += 1
            elif isinstance(gate, QutritYRotation):
                gate_dict["R_y^{" + gate.type + "}"] += 1
            elif isinstance(gate, QutritZRotation):
                gate_dict["R_z^{" + gate.type + "}"] += 1

            elif isinstance(gate, ControlledGate):
                controls.append(gate)
                state = gate.state
                name = gate.gate.name
                key = f"C_{state}_{name}"
                gate_dict[key] = gate_dict.get(key, 0) + 1
                #gate_dict["Controlled X"] += 1

        controls_dict = {}
        for cgate in controls:
            c = self.register[cgate.control]
            t = self.register[cgate.gate.target]
            key = f"{c}{t}"
            controls_dict[key] = controls_dict.get(key, 0) + 1

        #print(controls_dict)

        return gate_dict, controls_dict




    def draw(self, outputfile_name = None, maxgate = None):
        if maxgate is None:
            maxgate = len(self)


        if outputfile_name is None:
            ofile = open("qc_output.tex", "w")
        else:
            ofile = open(outputfile_name, "w")

        with open("latex_template.tex", "r") as template:
            preamble = template.readlines()

        lines = {}
        for pr in preamble:
            ofile.write(pr)

        gate_count = 0

        start_scale_box = r"\scalebox{0.7}{" + "\n"
        end_scale_box = r"}" + "\n"
        start_quantikz= r"\begin{quantikz}[row sep={0.7cm,between origins}]" + "\n"
        end_quantikz = r"\end{quantikz}"
        for gate in self.gates:
            if gate_count == 0:
                ofile.write(start_scale_box)
                ofile.write(start_quantikz)
                for i in range(self.nqubits + self.nqutrits):
                    lines[i] = r"\lstick{" + "$" + self.names[i] + "$ " + "}"

            if isinstance(gate, SingleQubitGate) or isinstance(gate, SingleQutritGate):
                target = gate.target
                for i in range(self.nqubits + self.nqutrits):
                    if i == target:
                        lines[i] += " & " + gate.to_latex()
                    else:
                        lines[i]  += r" & \qw "

            if isinstance(gate, ControlledGate):
                target = gate.gate.target
                control = gate.control
                state = gate.state


                state_sign = "[wire style = {" + f"\"{state}\"" + "}]"
                lines[control] += r" & \ctrl" + state_sign + "{" + f"{target - control}" + "}"
                lines[target] += r"& " + gate.gate.to_latex()
                for iii in range(self.nqubits + self.nqutrits):
                    if iii != control and iii != target:
                        lines[iii] += r" & \qw "


            if gate_count == maxgate - 1:
                for i in range(self.nqubits + self.nqutrits):
                    lines[i] += r" & \\" + "\n"
                    ofile.write(lines[i])

                ofile.write(end_quantikz)
                ofile.write(end_scale_box)
            gate_count = (gate_count + 1) % maxgate

        ofile.write(r"\end{document}")

        ofile.close()




if __name__ == "__main__":
    register = [3,2]
    qc = QuantumCircuit(register)
    xgate1 = SingleQubitGate(0, pauliX, "X")
    xgate2 = SingleQutritGate(1, X01, "X^{01}")
    rzgate = QubitZRotation(1,3.14)
    cxgate1 = ControlledGate(0, 2, xgate2)
    qc.insert_gate(xgate1)
    qc.insert_gate(xgate2)
    qc.insert_gate(cxgate1)
    qc.insert_gate(rzgate)
    qc.draw()

