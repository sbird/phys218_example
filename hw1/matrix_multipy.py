import time
import numpy as np
# test
# a = np.random.rand(75, 100)
# b = np.random.rand(100, 80)

a = [[1,2,3], [1,2,5], [2,4,6]]
b = [[1,2,3,4], [3,2,4,5], [9,0,7,6]]


def matr_mult_nested(a,b):
    import numpy as np
    
    assert np.shape(a)[1]==np.shape(b)[0]
    c = np.zeros((np.shape(a)[0],np.shape(b)[1]))
    for i in range(np.shape(c)[0]):
        for j in range(np.shape(c)[1]):
            for k in range(np.shape(a)[1]):
                c[i][j] += a[i][k] * b[k][j]
    return c

matr_mult_comp = []



tn1 = time.time()
nested = matr_mult_nested(a, b)
tn2 = time.time()
tn = tn2 - tn1
print("Nested in ", tn, "seconds.")

tnp1 = time.time()
nump = np.matmul(a,b)
tnp2 = time.time()
tnp = tnp2 - tnp1
print("Numpy in ", tnp, "seconds.")

nested == nump