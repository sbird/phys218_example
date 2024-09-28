import numpy as np
import time

# Define the matrix multiplication routines
def matrix_multiply_nested(A, B):
    rows_A = len(A)
    cols_A = len(A[0])
    cols_B = len(B[0])
    
    result = [[0 for _ in range(cols_B)] for _ in range(rows_A)]
    
    for i in range(rows_A):
        for j in range(cols_B):
            for k in range(cols_A):
                result[i][j] += A[i][k] * B[k][j]
    
    return result

def matrix_multiply_list_comp(A, B):
    cols_A = len(A[0])
    return [[sum(A[i][k] * B[k][j] for k in range(cols_A)) for j in range(len(B[0]))] for i in range(len(A))]

def matrix_multiply_numpy(A, B):
    return np.dot(A, B)

# Create a large random matrix and its inverse
size = 500  # Change size as needed
A = np.random.rand(size, size)
A_inv = np.linalg.inv(A)

# Convert A to a nested list for testing
A_nested = A.tolist()
A_inv_nested = A_inv.tolist()

# Test nested loops
start_time = time.time()
result_nested = matrix_multiply_nested(A_nested, A_inv_nested)
elapsed_nested = time.time() - start_time
print("Nested loops elapsed time:", elapsed_nested)

# Test list comprehension
start_time = time.time()
result_list_comp = matrix_multiply_list_comp(A_nested, A_inv_nested)
elapsed_list_comp = time.time() - start_time
print("List comprehension elapsed time:", elapsed_list_comp)

# Test NumPy
start_time = time.time()
result_numpy = matrix_multiply_numpy(A, A_inv)
elapsed_numpy = time.time() - start_time
print("NumPy elapsed time:", elapsed_numpy)

# Check if results are close to the identity matrix
identity = np.eye(size)

print("Nested loops result close to identity:", np.allclose(result_nested, identity.tolist()))
print("List comprehension result close to identity:", np.allclose(result_list_comp, identity.tolist()))
print("NumPy result close to identity:", np.allclose(result_numpy, identity))


#Nested loops elapsed time: 15.469426393508911
#List comprehension elapsed time: 12.796960592269897
#NumPy elapsed time: 0.0021791458129882812
#Nested loops result close to identity: True
#List comprehension result close to identity: True
#NumPy result close to identity: True
