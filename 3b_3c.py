'''
Problems 3b and 3c of Homework set 1 for PHYS 206
Following is the code to obtain the the matrix multiplication and
test it against the numpy routines
'''

import numpy as np
import time

#PART 1: NESTED FOR LOOP

def get_matrix_product_nested(A, B):
	'''
	Takes two matrices and gives its product using nested for loop
	'''
	assert len(A.shape) == 2 and len(B.shape) == 2, 'Please give 2D arrays'
	assert A.shape[1] == B.shape[0], 'Please give matrices which can be multiplied'
	


	C = np.zeros((A.shape[0], B.shape[1]))

	for j in range(B.shape[1]): #loop over the columns of the final matrix
		for i in range(A.shape[0]): #loop over the rows of the final matrix
			aik = A[i]
			bkj = B[:, j]
			cij = 0
			for k in range(A.shape[1]): #Summing all elements in the row - column multiplication for  each elements 
				cij += aik[k] * bkj[k]
			C[i][j] = cij
	return C


#PART 2: USING LIST 

def get_matrix_product_lc(A, B):
	'''
	Takes two matrices and gives a product using list comprehension technique
	'''
	C = [[sum([A[i][k]*B[:, j][k] for k in range(A.shape[1])]) for i in range(A.shape[0])] for j in range(B.shape[1])]

	return C




# Generating large random arrays:
size = 100 #This is the sixe of the array that you need

A = np.random.rand(100, 100)
while 1 == 1:
	if np.linalg.det(A) != 0:
		break
	else:
		A = np.random.rand(10, 10)

Ainv = np.linalg.inv(A)




# test of the time that it takes for the code to compile for different methods

t1 = time.time()

C1 = get_matrix_product_nested(A, Ainv)

t2 = time.time()

C2 = get_matrix_product_lc(A, Ainv)

t3 = time.time()

C3 = np.dot(A, Ainv)

t4 = time.time()

print(f'Time taken for nested loop is {round(t2-t1, 4)} seconds')
print(f'Time taken for list comprehension is {round(t3-t2, 4)} seconds')
print(f'Time taken for numpy is {round(t4-t3, 4)} seconds')




