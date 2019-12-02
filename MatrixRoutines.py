import numpy as np

def matrix_mult_loops(A, B):
    
    C = np.zeros([A.shape[1], B.shape[0]])

    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            for k in range(A.shape[0]):
                 C[i][j] += A[i][k]*B[k][j]
    
    return C 

def matrix_mult_list(A, B):
    
    C = [[sum(a*b for a,b in zip(A_row,B_col)) for B_col in zip(*B)] for A_row in A]

    return C

def matrix_mult_numpy(A,B):

    C = np.matmul(A,B)

    return C
