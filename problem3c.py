""" Write a simple routine to multiply two matrices together. Write one routine using nested for loops,
one using a list comprehension, and one using the built-in numpy matrix multiplication routines."""

import numpy as np

def for_mult(amatrix, bmatrix):
    """ amatrix and bmatrix are matrices (numpy arrays) that we will multiply using for loops """
    # first, make sure that we can actually multiply the matrices together
    a_rows = len(amatrix[:, 0])
    a_columns = len(amatrix[0, :])
    b_rows = len(bmatrix[:, 0])
    b_columns = len(bmatrix[0, :])
    assert a_columns == b_rows, "matrix dimensions do not match"
    cmatrix = np.zeros((a_rows, b_columns))
    # now we can multiply them by looping through the rows and columns:
    for i in range(a_rows):
        for j in range(b_columns):
            for k in range(b_rows):
                cmatrix[i, j] += amatrix[i, k] * bmatrix[k, j]
    return cmatrix

def comp_mult(amatrix, bmatrix):
    """ amatrix and bmatrix are matrices (numpy arrays) that we will multiply using list comprehensions """
    a_columns = range(len(amatrix[0, :]))
    b_rows = range(len(bmatrix[:, 0]))
    b_columns = range(len(bmatrix[0, :]))
    cmatrix = [[sum(amatrix[i, k] * bmatrix[k, j] for k in b_rows) for j in b_columns] for i in a_columns]
    return cmatrix

def np_mult(amatrix, bmatrix):
    """ multiply two matrices using built-in numpy function """
    cmatrix = np.matmul(amatrix, bmatrix)
    return cmatrix
