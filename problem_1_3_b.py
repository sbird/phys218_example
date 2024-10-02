#HW1, problem 1.3.b
#Writer: Aryana Haghjoo
#Oct 1, 2024

import numpy as np
import time
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
'''
# Create a large random matrix and its inverse
size = 500  # Change size as needed
A = np.random.rand(size, size)
A_inv = np.linalg.inv(A)

# Convert A to a nested list for testing
A_nested = A.tolist()
A_inv_nested = A_inv.tolist()

# Test nested loops
start_time = time.time()
result_nested = matrix_multiply_loop(A_nested, A_inv_nested)
elapsed_nested = time.time() - start_time
print("Nested loops elapsed time:", elapsed_nested)

# Test list comprehension
start_time = time.time()
result_list_comp = matrix_multiply_list(A_nested, A_inv_nested)
elapsed_list_comp = time.time() - start_time
print("List comprehension elapsed time:", elapsed_list_comp)

# Test NumPy
start_time = time.time()
result_numpy = matrix_multiply_numpy(A, A_inv)
elapsed_numpy = time.time() - start_time
print("NumPy elapsed time:", elapsed_numpy)

# Check if results are close to the identity matrix
#identity = np.eye(size)

#print("Nested loops result close to identity:", np.allclose(result_nested, identity.tolist()))
#print("List comprehension result close to identity:", np.allclose(result_list_comp, identity.tolist()))
#print("NumPy result close to identity:", np.allclose(result_numpy, identity))
'''
