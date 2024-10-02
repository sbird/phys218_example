import numpy as np
import pint as pt
import time

##### 1a.
ureg = pt.UnitRegistry()
G = ureg.constant('6.67430e-11 meter**3 / kilogram / second**2')  
c = ureg.constant('299792458 meter / second')  
M_sun = 1.989e30 * ureg.kilogram  
R_s = (2 * G * M_sun) / c**2
print({R_s})

###### 1b.
## nested for loop 

def matrix_multiply_loops(A, B):
    result = [[0 for _ in range(len(B[0]))] for _ in range(len(A))]
    for i in range(len(A)):
        for j in range(len(B[0])):
            for k in range(len(B)):
                result[i][j] += A[i][k] * B[k][j]
    return result

###### 1c.
# list comprehension

def matrix_multiply_comprehension(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(len(A[0]))) for j in range(len(B[0]))] for i in range(len(A))]
def matrix_multiply_numpy(A, B):
    return np.dot(A, B)

def create_matrix_and_inverse(size):
    A = np.random.rand(size, size)
    A_inv = np.linalg.inv(A)
    return A, A_inv

# Function to check if a matrix is approximately an identity matrix
def is_identity(matrix):
    identity = np.eye(len(matrix))
    return np.allclose(matrix, identity)

def time_function(func, A, A_inv):
    start_time = time.time()
    result = func(A, A_inv)
    end_time = time.time()
    return result, end_time - start_time

size = 300  
A, A_inv = create_matrix_and_inverse(size)

A_list = A.tolist()
A_inv_list = A_inv.tolist()


result_loops, time_loops = time_function(matrix_multiply_loops, A_list, A_inv_list)
print(f"Loops method time: {time_loops:.4f} seconds")
print(f"Is identity (loops): {is_identity(np.array(result_loops))}")

result_comprehension, time_comprehension = time_function(matrix_multiply_comprehension, A_list, A_inv_list)
print(f"Comprehension method time: {time_comprehension:.4f} seconds")
print(f"Is identity (comprehension): {is_identity(np.array(result_comprehension))}")

result_numpy, time_numpy = time_function(matrix_multiply_numpy, A, A_inv)
print(f"NumPy method time: {time_numpy:.4f} seconds")
print(f"Is identity (NumPy): {is_identity(result_numpy)}")

