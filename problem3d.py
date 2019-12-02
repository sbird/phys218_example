""" Simple test for Problem 3c: multiply a large matrix with its inverse to obtain the identity. """

# Results: all of the matrix multiplication implementations in Problem 3c work
# BUT they have wildly varying degrees of efficiency:
# For loops take a very long time, as expected
# List comprehensions take about the same amount of time as for loops, since
#    they are basically just for loops in disguise
# The built-in numpy function takes almost no time at all

import numpy as np
import problem3c

A = np.random.randint(10., size=(100, 100))

A_INV = np.linalg.inv(A)

def test_for():
    """ test the for loop multiplication with a big matrix """
    print("Testing for loop multiplication")
    I = problem3c.for_mult(A, A_INV)
    for i in range(len(I[0])):
        assert 0.9999 <= I[i][i] <= 1.0001, "Result is not the identity matrix"
        for j in range(len(I[0])):
            if j != i:
                assert abs(I[i][j]) <= 10**-12, "Result is not the identity matrix"

def test_comp():
    """ test the list comprehension multiplication from problem 3c """
    print("Testing list comprehension multiplication")
    I = problem3c.comp_mult(A, A_INV)
    for i in range(len(I[0])):
        assert 0.9999 <= I[i][i] <= 1.0001, "Result is not the identity matrix"
        for j in range(len(I[0])):
            if j != i:
                assert abs(I[i][j]) <= 10**-12, "Result is not the identity matrix"

def test_np():
    """ test the built-in numpy matrix multiplication """
    print("Testing built-in numpy.matmult")
    I = problem3c.np_mult(A, A_INV)
    for i in range(len(I[0])):
        assert 0.9999 <= I[i][i] <= 1.0001, "Result is not the identity matrix"
        for j in range(len(I[0])):
            if j != i:
                assert abs(I[i][j]) <= 10**-12, "Result is not the identity matrix"

if __name__ == "__main__":
    test_for()
    test_comp()
    test_np()
