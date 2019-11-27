#!/usr/bin/env python
"""Program to multiply 2  matrices together. (nxm)*(mxl)"""

import numpy as np 
import timeit

speed_test = True

def check_shape(A,B):
	"""Checking if it's possible to multiply these matrices together"""
	shape_A = np.shape(A)
	shape_B = np.shape(B)
	if shape_A[1] != shape_B[0]:
		raise ValueError("shapes {0} and {1} not aligned".format(shape_A,shape_B))

def matrix_mult_nested(A,B):
	"""Matrix multiplication using nested for loops"""
	check_shape(A,B)
	rows = np.shape(A)[0]
	cols = np.shape(B)[1]
	C = np.zeros((rows,cols)) #resultant matrix after multiplication
	for i in range(rows): 
		for j in range(cols): #loops give [i,j] location in matrix
			for k in range(cols): #each [i,j] location is a dot-product of a row of X and a column of Y
				C[i][j] += A[i][k]*B[k][j]
	return C

def matrix_mult_list_comp(A,B):
	"""Matrix multiplication using list comprehension"""
	check_shape(A,B)
	rows = np.shape(A)[0]
	cols = np.shape(B)[1]
	C = [[np.sum([A[i][k]*B[k][j] for k in range(cols)]) for j in range(cols)] for i in range(rows)]
	return C

def matrix_mult_numpy(A,B):
	"""Matrix multiplication using numpy dot function"""
	C = np.dot(A,B)
	return C

def time_elasped(function):
	start = timer()
	function
	end = timer()
	return end-start

#Example matrices
X = np.array([[1,2],[3,4],[5,6]])
Y = np.array([[5,6],[7,8]])

print("Example Matrices")
print("X * Y = \n{0} * \n{1}".format(X,Y))
print("Nested for loop:\n", matrix_mult_nested(X,Y))
print("List comprehension:\n", matrix_mult_list_comp(X,Y))
print("Numpy dot function:\n", matrix_mult_numpy(X,Y))


if speed_test: #speed test is optional
	#Large Matrix
	np.random.seed(1)
	large_X =  np.random.rand(20,20)
	large_X_inv = np.linalg.inv(large_X)
	print("\nSpeed test")

	#I only made these function because they don't have any arguments so I can use timeit on them easily
	def func1():
		return matrix_mult_nested(large_X,large_X_inv)
	def func2():
		return matrix_mult_list_comp(X,Y)
	def func3():
		return matrix_mult_numpy(X,Y)

	t_nested = np.min(timeit.repeat("func1()", setup="from __main__ import func1",repeat=5,number= 100))/100
	t_list_comp = np.min(timeit.repeat("func2()", setup="from __main__ import func2",repeat=5,number= 100))/100
	t_numpy = np.min(timeit.repeat("func3()", setup="from __main__ import func3",repeat=5,number= 100))/100
	print("Nested for loop:")
	print("{0:.3e} seconds".format(t_nested))
	print("List comprehension:")
	print("{0:.3e} seconds".format(t_list_comp))
	print("Numpy dot function:")
	print("{0:.3e} seconds".format(t_numpy))




#Code graveyard
	# print(timeit.repeat("func2()", setup="from __main__ import func2"))
	# print(timeit.repeat("func3()", setup="from __main__ import func3")
	# t_nested = time_elasped(matrix_mult_nested(large_X,large_X_inv))
	# t_list_comp = time_elasped(matrix_mult_list_comp(large_X,large_X_inv))
	# t_numpy = time_elasped(matrix_mult_numpy(large_X,large_X_inv))
	# print("Nested for loop:\n", t_nested*1e7)
	# #np.round(matrix_mult_nested(large_X,large_X_inv)))
	# print("List comprehension:\n", t_list_comp*1e7)
	# #np.round(matrix_mult_list_comp(large_X,large_X_inv)))
	# print("Numpy dot function:\n", t_numpy*1e7)
	# #np.round(matrix_mult_numpy(large_X,large_X_inv)))







