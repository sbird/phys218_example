# This is Q3_a

import numpy as np
import pint as pt

ureg = pt.UnitRegistry()
mass_of_sun = 1.989 * ureg.kilogram
G_constant = 6.67430 * ureg.newton * ureg.meter ** 2 / ureg.kilogram ** 2
c = ureg.speed_of_light
schwarzschild_r = (2 * G_constant * mass_of_sun) / (c ** 2)
print(schwarzschild_r.to(ureg.meter))
#######################################################################################

# This is Q3_b

# loops:

def matrix_multiply_loop(m1, m2):
    ro_1 = len(m1)
    ro_2 = len(m2)
    col_1 = len(m1[0])
    col_2 = len(m2[0])

    if col_1 != ro_2:
        print ("Error!")
    result = [[0] * col_2 for _ in range(ro_1)]
    for i in range(ro_1):
        for j in range(col_2):
            for k in range(col_1):
                result[i][j] += m1[i][k] * m2[k][j]
    return result


# list comprehension:

def matrix_multiply_comprehension(m1, m2):
    ro_1 = len(m1)
    ro_2 = len(m2)
    col_1 = len(m1[0])
    col_2 = len(m2[0])

    if col_1 != ro_2:
        print("Error!")

    result = [[sum(m1[i][k] * m2[k][j] for k in range(col_1)) for j in range(col_2)] for i in range(ro_1)]
    return result


# NumPy:

def matrix_multiply_numpy(m1, m2):
    return np.matmul(m1, m2)
#######################################################################################


# This is Q3_c

# Test

def test(routine):
    matrix = np.random.rand(1000, 1000)
    inverse = np.linalg.inv(matrix)
    time = timeit.timeit(lambda: routine(matrix, inverse), number=1)
    return time

print(matrix_multiply_numpy, "seconds")
print(matrix_multiply_comprehension, "seconds")
print(matrix_multiply_loop, "seconds")


