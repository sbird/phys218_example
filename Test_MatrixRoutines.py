import numpy as np
import MatrixRoutines
import time

#This test gives ~ identity matrix up to small rounding erros

#A = np.array([ [0.21141, 0.0088557, 0.8135, 0.688021, 0.590286],
#               [0.865644, 0.842035, 0.335114, 0.0222691, 0.963474],
#               [0.757321, 0.649845, 0.147549, 0.283823, 0.813359],
#               [0.150639, 0.917907, 0.241624, 0.154534, 0.923752],
#               [0.659617, 0.14401, 0.498398, 0.336458, 0.392511]])

#I = np.array([ [-0.657785, 0.0617966, 0.535632, -0.613242, 1.17083],
#               [-3.26064, -3.45166, -1.16934, 4.17959, 5.96285],
#               [-0.421356, 0.362159, -2.88574, 1.29481, 2.67726],
#               [-0.89975, -3.801, 1.97666, 1.68467, 2.62241],
#               [3.60801, 3.96089, 1.49873, -3.59111, -7.25504]])

A = np.random.randint(9, size=(10,10))
I = np.linalg.inv(A)

print("Matrix initialized, Starting multiplication")

start_time = time.time()
C1 = MatrixRoutines.matrix_mult_loops(A, I)
print("Loops takes: --- {0} seconds ---".format(time.time() - start_time))

start_time = time.time()
C2 = MatrixRoutines.matrix_mult_list(A, I)
print("Lists takes: --- {0} seconds ---".format(time.time() - start_time))

start_time = time.time()
C3 = MatrixRoutines.matrix_mult_numpy(A, I)
print("Numpy takes: --- {0} seconds ---".format(time.time() - start_time))

print("Loops: \n", C1)
print("List Comprehention: \n", C2)
print("Numpy: \n", C3)
