import time

import numpy as np

np.random.seed(2046)



def loop_matrix_multiply(m1, m2):
    result = np.zeros((m1.shape[0], m2.shape[1]))
    for i in range(m1.shape[0]):
        for j in range(m2.shape[1]):
            for k in range(m1.shape[1]):
                result[i, j] += m1[i, k] * m2[k, j]
    return result

def test_loop(m1, m2):
    start = time.time()
    assert np.allclose(loop_matrix_multiply(m1, m2), np.eye(m1.shape[0]))
    end = time.time()

    print(f"Loop Method: {end - start:.6f} seconds")

def list_matrix_multiply(m1, m2):
    return np.array([[sum(m1[i, :] * m2[:, j]) for j in range(m2.shape[1])] for i in range(m1.shape[0])])

def test_list(m1, m2):
    start = time.time()
    assert np.allclose(list_matrix_multiply(m1, m2), np.eye(m1.shape[0]))
    end = time.time()

    print(f"List Method: {end - start:.6f} seconds")

def numpy_matrix_multiply(m1, m2):
    return np.dot(m1, m2)

def test_numpy(m1, m2):
    start = time.time()
    assert np.allclose(numpy_matrix_multiply(m1, m2), np.eye(m1.shape[0]))
    end = time.time()

    print(f"Numpy Method: {end - start:.6f} seconds")

if __name__ == "__main__":
    m1 = np.random.rand(500, 500)
    m2 = np.linalg.inv(m1)
    
    test_loop(m1, m2)
    test_list(m1, m2)
    test_numpy(m1, m2)

