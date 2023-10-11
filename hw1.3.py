# First Question:
import pint

# Define the units required
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

def multiply_matrices (matrix_a,matrix_b):
    #Check if the matrices can be multiplied.
    if len(matrix_a[0]) != len(matrix_b):
        return "Can not multiply matrices. Invalid dimensions"
    #Initialize the result matrix with zeros.
    result = []
    for i in range(len(matrix_a)):
        row = []
        for j in range(len(matrix_b[0])):
            row.append(0)
        result.append(row)
    #Iterate through the rows of the first matrix:
    for i in range(len(matrix_a)):
        #Iterate through the columns of the second matrix:
        for j in range(len(matrix_b[0])):
            #Iterate through the elements of each row in the second matrix and do the multiplication:
            for k in range(len(matrix_b)):
                result[i][j] += matrix_a[i][k] * matrix_b[k][j]
    return result
# Example usage:
matrix_a = [[1, 2, 3], [4, 5, 6]]
matrix_b = [[7, 8], [9, 10], [11, 12]]
result= multiply_matrices (matrix_a,matrix_b)
print(result)


# list method:
def multiply_matrices_listcomprehension(matrix_a,matrix_b):
    #Check if the matrices can be multiplied.
    if len(matrix_a[0]) != len(matrix_b):
        return "Can not multiply matrices. Invalid dimensions"
    results=[[sum(matrix_a[i][k] * matrix_b[k][j] for k in range(len(matrix_b)))for j in range(len(matrix_b[0]))] for i in range(len(matrix_a))]
    return results
#Example usage:
matrix_a = [[1, 2, 3], [4, 5, 6]]
matrix_b = [[7, 8], [9, 10], [11, 12]]
results= multiply_matrices_listcomprehension (matrix_a,matrix_b)
print(results)

# numpy method:
import numpy as np

def multiply_numpy(matrix_a, matrix_b):
    # Convert input matrices to NumPy arrays
    array1 = np.array(matrix_a)
    array2 = np.array(matrix_b)
                    

    # Perform matrix multiplication using NumPy
    result = np.dot(array1, array2)

    # Convert the result back to a nested list (if needed)
    result = result.tolist()

    return result

#Example usage:
matrix_a = [[1, 2, 3], [4, 5, 6]]
matrix_b = [[7, 8], [9, 10], [11, 12]]
result= multiply_numpy (matrix_a,matrix_b)
print(result)


#Third Question:

import timeit

# Define the size of the square matrix
matrix_size = 100  # You can adjust the size as needed

# Create a random square matrix
matrix = np.random.rand(matrix_size, matrix_size)

# Calculate the inverse of the matrix
matrix_inverse = np.linalg.inv(matrix)

# Define a function to test a matrix multiplication routine
def test_matrix_multiplication(matrix, matrix_inverse, multiplication_function):
    result = multiplication_function(matrix, matrix_inverse)
    # Check if the result is close to the identity matrix (within a small tolerance)
    return np.allclose(result, np.identity(matrix_size))

# Perform the tests
nested_for_result = test_matrix_multiplication(matrix, matrix_inverse, multiply_matrices)
list_comp_result = test_matrix_multiplication(matrix, matrix_inverse, multiply_matrices_listcomprehension)
numpy_result = test_matrix_multiplication(matrix, matrix_inverse, multiply_numpy)

# Check if all tests passed and print a message
if nested_for_result and list_comp_result and numpy_result:
    print("All matrix multiplication routines passed the test.")
else:
    print("One or more matrix multiplication routines did not pass the test.")

# Measure the execution time for each routine
nested_for_time = timeit.timeit(lambda: test_matrix_multiplication(matrix, matrix_inverse, multiply_matrices), number=10)
list_comp_time = timeit.timeit(lambda: test_matrix_multiplication(matrix, matrix_inverse, multiply_matrices_listcomprehension), number=10)
numpy_time = timeit.timeit(lambda: test_matrix_multiplication(matrix, matrix_inverse, multiply_numpy), number=10)

# Print the execution times
print("Execution time for matrix_multiply_nested_loop: {:.6f} seconds".format(nested_for_time))
print("Execution time for matrix_multiply_list_comprehension: {:.6f} seconds".format(list_comp_time))
print("Execution time for matrix_multiply_numpy: {:.6f} seconds".format(numpy_time))