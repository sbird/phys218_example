import numpy as np
from matrix_multipliers import loop_multiplier, list_comp_muliplier, numpy_multiplier

m = np.random.rand(10,10)
n = np.linalg.inv(m)

print(loop_multiplier(m,n))

print(list_comp_muliplier(m,n))

print(numpy_multiplier(m,n))

