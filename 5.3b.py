import numpy as np
"""
Matrix Multiplication using different methods
"""

"""Defining the matrices"""
A=[[1,2,3],[2,3,4],[1,4,2]]
B=[[1,4,2],[1,1,2],[2,3,1]]
print("Matrix A")
for i in range(3):
    for j in range(3):
        print(A[i][j], end = " ")
    print()
print("Matrix B")
for i in range(3):
    for j in range(3):
        print(B[i][j], end = " ")
    print()

def Nested_For(A,B):
    """Using nested for"""
    result=[[0,0,0],[0,0,0],[0,0,0]]
    for i in range(len(A)):
        for j in range(len(B[0])):
            for k in range(len(B)):
                result[i][j] += A[i][k] * B[k][j]
    for r in result:
        print(r)


def List(A,B):
    """Using list comprehension"""
    result = [[sum(a * b for a, b in zip(A_row, B_col))  for B_col in zip(*B)] for A_row in A]
    for r in result:
        print(r)

def Numpy(A,B):
    """Using Numpy"""
    print(np.dot(A,B))

print("Matrix multiplication using Nested for ")
Nested_For(A,B)
print("Matrix multiplication using List Comprehension ")
List(A,B)
print("Matrix Multiplication using Numpy ")
Numpy(A,B)
