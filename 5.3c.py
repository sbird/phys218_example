import numpy as np
import time
import math
"""Defining the matrices"""
C=np.array([[1,2,3,1,3],[2,3,4,5,5],[4,1,6,2,7],[1,2,6,1,9],[3,2,1,5,6]])
print("The first matrix is:")
print(C)
D=np.linalg.inv(C)
print("Its Inverse matrix is:")
print(D)

def Nested_For(A,B):
    """Using nested for"""
    start= time.time()
    result=np.zeros([5,5])
    for i in range(len(A)):
        for j in range(len(B[0])):
            for k in range(len(B)):
                result[i][j] += A[i][k] * B[k][j]
    result[0][0]=math.ceil(result[0][0])
    for r in result:
        print(r.astype(int))
    end=time.time()
    print("Time taken in this routine:"+str(end-start)+" s")


def List(A,B):
    """Using list comprehension"""
    start= time.time()
    result = [[sum(a * b for a, b in zip(A_row, B_col))  for B_col in zip(*B)] for A_row in A]
    result[0][0]=math.ceil(result[0][0])
    result=np.array(result).astype(int)
    for r in result:
        print(r)
    end=time.time()
    print("Time taken in this routine:"+str(end-start) + " s")

#Using numpy
def Numpy(A,B):
    """Using Numpy"""
    start= time.time()
    res=np.dot(A,B)
    res[0][0]=math.ceil(res[0][0])
    print(res.astype(int))
    end=time.time()
    print("Time taken in this routine:"+str(end-start) + " s")

print("Matrix multiplication using Nested for ")
Nested_For(C,D)
print("Matrix multiplication using List Comprehension ")
List(C,D)
print("Matrix Multiplication using Numpy ")
Numpy(C,D)
