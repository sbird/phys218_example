import numpy as np
from matrices import *


def test_multiply_matrices():
    A = np.random.rand(50,50) + np.eye(50) # A 50x50 matrix containing random integers
    B = np.linalg.inv(A) # Matrix inverse

    assert np.allclose(multiply_matrices(A,B), np.identity(len(A))), "Matrix multiplication fails."

def test_multiply_matrices_list():
    A = np.random.rand(50,50) + np.eye(50) # A 50x50 matrix containing random integers
    B = np.linalg.inv(A) # Matrix inverse

    assert np.allclose(multiply_matrices_list(A,B), np.identity(len(A))), "Matrix multiplication fails."

def test_multiply_matrices_np():
    A = np.random.rand(50,50) + np.eye(50) # A 50x50 matrix containing random integers
    B = np.linalg.inv(A) # Matrix inverse

    assert np.allclose(multiply_matrices_np(A,B), np.identity(len(A))), "Matrix multiplication fails."