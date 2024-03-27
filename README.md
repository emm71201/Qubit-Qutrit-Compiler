# Hybrid Quantum Circuit Compiler: Qubit & Qutrit

##  Quick Tutorial
1. Write unitary matrix to a CSV File
2. Define the register as a Comma separated sequence of 3 (qutrit) and 2 (qubit). 
3. Run the command `python3 main.py -f <Path to Matrix File> -r register`

Example: `python3 main.py -f matrices/s108fft-108dim.csv -r 3,3,3,2,2`

4. Possible Errors:
    * The Matrix File does not exist
    * Numpy Failed to Read the Matrix File
    * The dimension of the Matrix File does not match the product of the register
    * The Matrix provided is not Unitary
    * Packages not installed

