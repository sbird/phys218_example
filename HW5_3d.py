import numpy as np
from numpy.linalg import inv
import time


N= 100
A = np.random.randint(0, 9, size=(N,N))
B = np.linalg.inv(A)
#print("A = ",A)
#print("B = ",B)
#for loops
s = (N, N)
D = np.zeros(s)
start = time.process_time()
for i in range(0, N):
    for j in range(0, N):
        for k in range(0, N):
            D[i][j] += np.multiply(A[i][k], B[k][j])
print(time.process_time() - start)
print(D)
#list comprehension
start = time.process_time()
E = [[sum(a*b for a,b in zip(A_row,B_col)) for B_col in zip(*B)] for A_row in A]
print(time.process_time() - start)
for r in E:
    print(r)
#numpy
start = time.process_time()
C = np.matmul(A, B)
print(time.process_time() - start)
print(C)
#as we can see by the results the numpy routine is the fastest way:)