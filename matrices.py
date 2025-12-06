import time
import numpy as np


def multiply_matrices(A: list[list[float]], B: list[list[float]]) -> list[list[float]]:
    """
    Multiply two matrices A and B using nested for loops.
    
    :param A: Matrix of size l x m
    :type A: list[list[float]]
    :param B: Matrix of size m x n
    :type B: list[list[float]]
    :return: Matrix C of size l x n
    :rtype: list[list[float]]
    """

    assert len(A[0]) == len(B), "The columns of A do not match the rows of B."

    C = []
    for row in range(len(A)):
        C.append([0] * len(B[0]))

    for row in range(len(A)):
        for column in range(len(B[0])):
            for k in range(len(B)):
                C[row][column] += A[row][k] * B[k][column]

    return C

def multiply_matrices_list(A: list[list[float]], B: list[list[float]]) -> list[list[float]]:
    """
    Multiply two matrices A and B using nested list comprehension.
    
    :param A: Matrix of size l x m
    :type A: list[list[float]]
    :param B: Matrix of size m x n
    :type B: list[list[float]]
    :return: Matrix C of size l x n
    :rtype: list[list[float]]
    """

    assert len(A[0]) == len(B), "The columns of A do not match the rows of B."

    C = [[sum(a*b for a, b in zip(A_row, B_column)) for B_column in zip(*B)] for A_row in A]

    return C

def multiply_matrices_np(A: list[list[float]], B: list[list[float]]) -> list[list[float]]:
    """
    Multiply matrices using numpy.
    
    :param A: Matrix of size l x m
    :type A: list[list[float]]
    :param B: Matrix of size m x n
    :type B: list[list[float]]
    :return: Matrix C of size l x n
    :rtype: list[list[float]]
    """

    assert np.shape(A)[1] == np.shape(B)[0], "The columns of A do not match the rows of B."

    C = np.dot(A, B)

    return C

# compute time
A = np.random.rand(50,50) + np.eye(50) # A 50x50 matrix containing random integers
B = np.linalg.inv(A) # Matrix inverse

start_for = time.time()
print(np.allclose(multiply_matrices(A, B), np.identity(len(A))))
end_for = time.time()
print("For loop time:", end_for - start_for, "seconds")

start_list = time.time()
print(np.allclose(multiply_matrices_list(A, B), np.identity(len(A))))
end_list = time.time()
print("List time:", end_list - start_list, "seconds")

start_np = time.time()
print(np.allclose(multiply_matrices_np(A, B), np.identity(len(A))))
end_np = time.time()
print("np time:", end_np - start_np, "seconds")