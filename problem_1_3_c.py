#HW1, problem 1.3.c
#Writer: Aryana Haghjoo
#Oct 1, 2024

import numpy as np
import time

#importing the functions from problem 1.3.b
from problem_1_3_b import matrix_multiply_loop, matrix_multiply_list, matrix_multiply_numpy

#testing the functions with a big random matrix
# We also calculate the inverse of this matrix
size = 500
A = np.random.rand(size, size)
A_inv = np.linalg.inv(A)

# since fundtions "matrix_multiply_loop" and "matrix_multiply_list" work with nested lists, we convert A and A_inv to nested lists
A_list = A.tolist()
A_inv_list = A_inv.tolist()

#testing "matrix_multiply_loop"
time_initial = time.time()
result_loop = matrix_multiply_loop(A_list, A_inv_list)
delta_time = time.time() - time_initial
print(f"Nested loops Calculation time is {delta_time:.2f} seconds.")

#testing "matrix_multiply_list"
time_initial = time.time()
result_list = matrix_multiply_list(A_list, A_inv_list)
delta_time = time.time() - time_initial
print(f"List Comprehension Calculation time is {delta_time:.2f} seconds.")

#testing "matrix_multiply_numpy"
time_initial = time.time()
result_numpy = matrix_multiply_numpy(A_list, A_inv_list)
delta_time = time.time() - time_initial
print(f"Numpy Calculation time is {delta_time:.2f} seconds.")