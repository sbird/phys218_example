#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 17:05:06 2020

@author: wenjunchang
"""
from pint import UnitRegistry
ureg = UnitRegistry()
ureg.define("Msolar = 1.98855*10**30 * kilogram")

#The unit of inputing mass is Msolar
def Calcu_Schwarzchild_radius(mass):
    mass = mass * ureg.Msolar
    r = ureg.gravitational_constant * mass * 2 / ureg.c**2
    return r.to_base_units()
    
import numpy as np
def matrix_multi_loop(matrixA, matrixB):
    solv = np.zeros((len(matrixA), len(matrixB[0])), dtype=float)
    for j in range(len(matrixA)):
        for i in range(len(matrixB[0])):
            for k in range(len(matrixA[0])):
                solv[j][i] += matrixA[j][k] * matrixB[k][i]                
    return solv

def matrix_multi_comp(matrixA, matrixB):
    solv = [[sum(a*b for a,b in zip(X_row,Y_col)) for Y_col in zip(*matrixB)] for X_row in matrixA]
    return np.array(solv)

def matrix_multi_np(matrixA, matrixB):
    return np.dot(matrixA, matrixB)

###Check    
import datetime

def check_matrix_loop(matrixA):
    A_T = np.linalg.inv(matrixA)
    starttime = datetime.datetime.now()
    solv = matrix_multi_loop(matrixA, A_T)
    endtime = datetime.datetime.now()
    return solv, (endtime - starttime)

def check_matrix_comp(matrixA):
    A_T = np.linalg.inv(matrixA)
    starttime = datetime.datetime.now()
    solv = matrix_multi_comp(matrixA, A_T)
    endtime = datetime.datetime.now()
    return solv, (endtime - starttime)

def check_matrix_np(matrixA):
    A_T = np.linalg.inv(matrixA)
    starttime = datetime.datetime.now()
    solv = matrix_multi_np(matrixA, A_T)
    endtime = datetime.datetime.now()
    return solv, (endtime - starttime)

#check
if __name__=='__main__':    
    
    mass = 1 #unit of Msolar
    radius = Calcu_Schwarzchild_radius(mass)
    print('The Schwarzchild radius of the sun is: ', radius)
    
    len_matrix = 20
    matrix_A = np.random.rand(len_matrix,len_matrix)
    solv_loop, time_loop = check_matrix_loop(matrix_A)
    solv_comp, time_comp = check_matrix_comp(matrix_A)
    solv_np, time_np = check_matrix_np(matrix_A)
    print('\ntime loop: ',time_loop)
    print('time comprehension: ',time_comp)
    print('time np multiplication: ',time_np)