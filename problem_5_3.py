import pint
def sr_of_earth():
    ureg = pint.UnitRegistry(system='mks')#set default system so that the unit of radius is m
    G = 6.67 * 10**(-11) * ureg.N * ureg.meter**2 / ureg.kg**2#gravitaitonal constant
    m = 2 * 10**30 * ureg.kg#solar mass
    c = 3 * 10**8 * ureg.meter / ureg.s#speed of light
    r = 2 * G * m / c**2#Schwarzschild radius
    print(r.to_base_units())
sr_of_earth()

import halo_mass_function as hmf
import numpy as np
def halonumber(n):
	halo = hmf.HaloMassFunction(redshift=0)#set redshift = 0
	h = 0.67#hubble constant
	mass = np.array([12 + i * 8 / n for i in range(0, n+1)]) / h#upper limit of halo mass is 20 Msolar/h
	dn = halo.dndm(mass)
	number = (sum(2 * dn) - dn[0] -dn[-1]) * 8 / (2 * n * h)#use sum to calculate integral.
	print("Hubble constant is %.2f.\nThe number of halos per (Mpc/h)^3 are %.0f"%(h, number))
halonumber(10)

def Nestedfor(A,B):
	m = A.shape[0]
	n = A.shape[1]
	p = B.shape[1]
	C = np.zeros((m,p))
	for i in range(m):
		for j in range(p): 
			for k in range(n):
				C[i][j] += A[i][k] * B[k][j]
	return C

def ListComprehension(A,B):
	m = A.shape[0]
	n = A.shape[1]
	p = B.shape[1]
	C = np.array([[sum([A[i][k] * B[k][j] for k in range(n)]) for j in range(p)] for i in range(m)])
	return C

def builtin(A,B):
	return np.dot(A,B)
import time
def test(m):
	A = np.random.random((m,m))
	B = np.linalg.inv(A)
	start = time.time()
	Nestedfor(A,B)
	end = time.time()
	print("Running time using nested for:%fs"%(end - start))
	start = time.time()
	ListComprehension(A,B)
	end = time.time()
	print("Running time using list comprehension:%fs"%(end - start))
	start = time.time()
	builtin(A,B)
	end = time.time()
	print("Running time using built-in numpy function:%fs"%(end - start))

test(50)