import sys
import os
import numpy as np
sys.path.append("src/")
import compile

if __name__ == "__main__":
    filename = None
    register = None
    for i in range(len(sys.argv)):
        if sys.argv[i] == "-f":
            filename = sys.argv[i+1]
            if not os.path.isfile(filename):
                raise Exception(f"File '{filename}' Does Not Exist\nExiting...")

        if sys.argv[i] == "-r":
            register = [int(ch) for ch in sys.argv[i+1].split(",")]

    if register[-1] == 3:
        raise Exception("The general single qutrit operator is not yet supported.\nAt present, the last regsiter has to be qubit.")



    print("*" * 100, "\n")
    print(" "*40, "Hybrid Qubit - Qutrit Compiler ...\n")
    print("*"*100)
    print("\n\t1. Reading and Validating data")
    print("\n\t\ta. Matrix File: ", filename)
    print("\n\t\tb. Registers: ", register)

    try:
        matrix = np.genfromtxt("matrices/s108fft-108dim.csv", delimiter=",", dtype=np.complex_)
        if np.prod(register) != matrix.shape[0]:
            raise Exception(f" Register dimension {np.prod(register)} does not match matrix dimension {matrix.shape[0]} ")
    except:
        print("Failed to read the matrix file.")

    print("\n\t2. Starting Compilation")

    compile.compile(matrix, register)

    print("\n\t3. Process Completed. Exiting.\n")
    print("*" * 100, "\n")

