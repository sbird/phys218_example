from pint import UnitRegistry
import numpy as np
ureg = UnitRegistry()

#PART B

def loop_multiplier(X,Y):
    """X and Y are matrices being multiplied,
        this is a nested for loop technique"""
    if len(X[0]) != len(Y):
        print("Matrices are wrong dimensions to be multiplied.")
    else:
        Z = np.zeros((len(X),len(Y[0])))
        for i in range(len(X)):
            for j in range(len(Y[0])):
                for k in range(len(X[0])):
                    Z[i][j] += X[i][k] * Y[k][j]
                #rounding for the test
                Z[i][j] = round(Z[i][j])
        return(Z)


def list_comp_muliplier(X,Y):
    """X and Y are matrices being multiplied,
        this is a list comprehension technique"""
    if len(X[0]) != len(Y):
        print("Matrices are wrong dimensions to be multiplied.")
    else:
        Z = [[round(sum(i*j for i,j in zip(row,col)),4) for col in zip(*Y)] for row in X]
        return(Z)


def numpy_multiplier(X,Y):
    """X and Y are matrices being multiplied,
        this used numpy's built-in matrix multiplication"""
    if len(X[0]) != len(Y):
        print("Matrices are wrong dimensions to be multiplied.")
    else:
        Z = np.matmul(X,Y)
        #rounding for the test 
        for i in range(len(Z)):
            for j in range(len(Z[0])):
                Z[i][j] = round(Z[i][j])
        return(Z)