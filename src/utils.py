import numpy as np
import scipy as sp

def modulus(z):
    return np.sqrt(z.real**2 + z.imag**2)

def chop(z, tol=1e-12):
    """
    Returns 0 if the magnitude of the complex number z is less than tol,
    otherwise returns z.

    Parameters:
    - z: A complex number.
    - tol: Tolerance level.

    Returns:
    - A complex number or 0.
    """
    # Calculate the magnitude of the complex number
    magnitude = modulus(z)

    # Compare the magnitude to the tolerance
    if magnitude < tol:
        return 0

    if np.abs(z.imag) < tol:
        return z.real
    if np.abs(z.real) < tol:
        return 1j * z.imag
    return z
def direct_sum(*matrices):
    # Calculate the total size of the resulting matrix
    total_rows = sum(matrix.shape[0] for matrix in matrices)
    total_cols = sum(matrix.shape[1] for matrix in matrices)

    # Initialize the result matrix with zeros
    result = np.zeros((total_rows, total_cols), dtype=np.complex128)

    # Place each matrix in its block-diagonal position
    row_start, col_start = 0, 0
    for matrix in matrices:
        rows, cols = matrix.shape
        result[row_start:row_start+rows, col_start:col_start+cols] = matrix
        row_start += rows
        col_start += cols

    return result

def off_diagonal_direct_sum(*matrices):
    # Determine the total size of the resulting matrix
    total_rows = sum(matrix.shape[0] for matrix in matrices)
    total_cols = sum(matrix.shape[1] for matrix in matrices)

    # Initialize the result matrix with zeros
    result = np.zeros((total_rows, total_cols), dtype=np.complex128)

    # Place each matrix in its off-diagonal position
    row_end, col_start = total_rows, 0
    for matrix in matrices:
        rows, cols = matrix.shape
        row_start = row_end - rows
        result[row_start:row_end, col_start:col_start + cols] = matrix
        row_end = row_start
        col_start += cols

    return result

def is_diagonal_matrix(matrix, tol=1e-12):
    # Check if the matrix is square
    if matrix.shape[0] != matrix.shape[1]:
        return False

    # Create a boolean mask for the diagonal elements
    diagonal_mask = np.eye(matrix.shape[0], dtype=bool)

    # Invert the mask to get non-diagonal elements
    non_diagonal_mask = ~diagonal_mask

    # Check if the magnitude of all non-diagonal elements is close to zero within the tolerance
    # np.abs is used to get the magnitude of complex numbers
    if np.all(np.abs(matrix[non_diagonal_mask]) < tol):
        return True

    return False

def rot_matrix(angle, generator):
    return sp.linalg.expm(-1j * angle * generator)