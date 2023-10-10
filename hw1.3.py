# First Question:
import pint

# Define the units and constants
ureg = pint.UnitRegistry()
G = 6.67430e-11 * ureg.meter**3 / (ureg.kilogram * ureg.second**2)  # Gravitational constant
M = 1.989e30 * ureg.kilogram  # Mass of the Sun
c = 3e8 * ureg.meter / ureg.second  # Speed of light

# Calculate the Schwarzschild radius
Schwarzschild_radius = 2 * G * M / c**2

# Print the result in meters
print(f"The Schwarzschild radius of the Sun is approximately {Schwarzschild_radius.to(ureg.meter)}")

# Second Question:
# for method:
def matrix_multiply_nested_loop(matrix1, matrix2):
    # Get the dimensions of the input matrices
    rows1, cols1 = len(matrix1), len(matrix1[0])
    rows2, cols2 = len(matrix2), len(matrix2[0])

    # Check if the matrices can be multiplied
    if cols1 != rows2:
        raise ValueError("Number of columns in matrix1 must be equal to the number of rows in matrix2")

    # Initialize the result matrix with zeros
    result = [[0 for _ in range(cols2)] for _ in range(rows1)]

    # Perform matrix multiplication using nested loops
    for i in range(rows1):
        for j in range(cols2):
            for k in range(cols1):
                result[i][j] += matrix1[i][k] * matrix2[k][j]

    return result
# list method:
def matrix_multiply_list_comprehension(matrix1, matrix2):
    # Get the dimensions of the input matrices
    rows1, cols1 = len(matrix1), len(matrix1[0])
    rows2, cols2 = len(matrix2), len(matrix2[0])

    # Check if the matrices can be multiplied
    if cols1 != rows2:
        raise ValueError("Number of columns in matrix1 must be equal to the number of rows in matrix2")

    # Perform matrix multiplication using list comprehension
    result = [[sum(matrix1[i][k] * matrix2[k][j] for k in range(cols1)) for j in range(cols2)] for i in range(rows1)]

    return result
# numpy method:
import numpy as np

def matrix_multiply_numpy(matrix1, matrix2):
    # Convert input matrices to NumPy arrays
    array1 = np.array(matrix1)
    array2 = np.array(matrix2)

    # Perform matrix multiplication using NumPy
    result = np.dot(array1, array2)

    # Convert the result back to a nested list (if needed)
    result = result.tolist()

    return result

##Third Question:
import numpy as np
import timeit

# Define the size of the square matrix
matrix_size = 100  # You can adjust the size as needed

# Create a random square matrix
matrix = np.random.rand(matrix_size, matrix_size)

# Calculate the inverse of the matrix
matrix_inverse = np.linalg.inv(matrix)

# Define a function to test each multiplication routine
def test_matrix_multiplication(matrix, matrix_inverse, multiplication_function):
    result = multiplication_function(matrix, matrix_inverse)
    # Check if the result is close to the identity matrix (within a small tolerance)
    return np.allclose(result, np.identity(matrix_size))

# Measure the execution time for each routine
nested_for_time = timeit.timeit(lambda: test_matrix_multiplication(matrix, matrix_inverse, matrix_multiply_nested_loop), number=10)
list_comp_time = timeit.timeit(lambda: test_matrix_multiplication(matrix, matrix_inverse, matrix_multiply_list_comprehension), number=10)
numpy_time = timeit.timeit(lambda: test_matrix_multiplication(matrix, matrix_inverse, matrix_multiply_numpy), number=10)

# Print the execution times
print("Execution time for matrix_multiply_nested_loop: {:.6f} seconds".format(nested_for_time))
print("Execution time for matrix_multiply_list_comprehension: {:.6f} seconds".format(list_comp_time))
print("Execution time for matrix_multiply_numpy: {:.6f} seconds".format(numpy_time))

