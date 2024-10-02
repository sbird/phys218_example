"Write a simple routine to multiply two matrices together in python"
import numpy as np
import time

"routine using nested for loops"
def matrix_mul_1(Matrix_A,Matrix_B):
    rowA=len(Matrix_A);colA=len(Matrix_A[0])
    rowB=len(Matrix_B);colB=len(Matrix_B[0])

    Matrix_C=[]
    for i in range(rowA):
        row_C=[]
        for j in range(colB):
            sum=0
            for m in range(colA):
                sum+=Matrix_A[i][m]*Matrix_B[m][j]
            row_C.append(sum)
        Matrix_C.append(row_C)

    return np.array(Matrix_C)

"routine using a list comprehension"
def matrix_mul_2(Matrix_A,Matrix_B):
   return np.array([[np.sum(np.array([Matrix_A[j][m]*Matrix_B[m][i] for m in range(len(Matrix_A[0]))])) for i in range(len(Matrix_B[0]))] for j in range(len(Matrix_A))])

"routine using the built-in numpy matrix multiplication"
def matrix_mul_3(Matrix_A,Matrix_B):
    return np.dot(Matrix_A,Matrix_B)


np.set_printoptions(threshold=10)
"test"
"generate dimension=dim diagonal matrix so that it is inversible."
dim=100
mA=np.diag([i+1 for i in range(dim)])
mB=np.linalg.inv(mA)
print("Matrix A is")
print(mA)
print("Matrix B is")
print(mB)

"test1"
print("test1 for routine using nested for loops")
print("Multiplication of A and B is")
stime=time.time()
print(matrix_mul_1(mA,mB))
etime=time.time()
print("time="+"%.5f"%(etime-stime)+" s")
print('\n')

"test2"
print("test2 using a list comprehension")
print("Multiplication of A and B is")
stime=time.time()
print(matrix_mul_2(mA,mB))
etime=time.time()
print("time="+"%.5f"%(etime-stime)+" s")
print('\n')

"test3"
print("test3 using the built-in numpy matrix multiplication")
print("Multiplication of A and B is")
stime=time.time()
print(matrix_mul_3(mA,mB))
etime=time.time()
print("time="+"%.5f"%(etime-stime)+" s")