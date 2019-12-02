import numpy as np

N= 5
A = np.random.randint(0, 9, size=(N,N))
B = np.random.randint(0, 9, size=(N,N))
#for loops
s = (N, N)
D = np.zeros(s)
for i in range(0, N):
    for j in range(0, N):
        for k in range(0, N):
            D[i][j] += np.multiply(A[i][k], B[k][j])
print(D)
#list comprehension
E = [[sum(a*b for a,b in zip(A_row,B_col)) for B_col in zip(*B)] for A_row in A]
for r in E:
    print(r)
#numpy
C = np.matmul(A, B)
print(C)