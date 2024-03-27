import numpy as np

def is_diagonal_matrix_complex(matrix, tol=1e-12):
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

if __name__ == "__main__":

    # Example usage with complex numbers
    matrix = np.array([[1+1j, 0, 0], [0, 2+2j, 0], [0, 0, 3+3j]])
    print(is_diagonal_matrix_complex(matrix))  # This should print True

    matrix = np.array([[1+1j, 1e-13, 0], [0, 2+2j, 0], [1e-13j, 0, 3+3j]])
    print(is_diagonal_matrix_complex(matrix))  # This should print True

    matrix = np.array([[1+1j, 1e-10, 0], [0, 2+2j, 0], [0, 0, 3+3j]])
    print(is_diagonal_matrix_complex(matrix))  # This should print False, depending on tol
