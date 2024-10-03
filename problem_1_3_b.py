#HW1, problem 1.3.b
#Writer: Aryana Haghjoo
#Oct 1, 2024

import numpy as np
#First, I define three functions to multiply two matrices: matrix_multiply_loop, matrix_multiply_list, and matrix_multiply_numpy.
# Define the matrix multiplication routines
def matrix_multiply_loop(A, B):
    n_rows_A = len(A)
    n_cols_A = len(A[0])
    n_cols_B = len(B[0])
    
    result = [[0 for _ in range(n_cols_B)] for _ in range(n_rows_A)]
    
    for i in range(n_rows_A):
        for j in range(n_cols_B):
            for k in range(n_cols_A):
                result[i][j] += A[i][k] * B[k][j]
    
    return result

def matrix_multiply_list(A, B):
    n_cols_A = len(A[0])
    return [[sum(A[i][k] * B[k][j] for k in range(n_cols_A)) for j in range(len(B[0]))] for i in range(len(A))]

def matrix_multiply_numpy(A, B):
    return np.dot(A, B)

