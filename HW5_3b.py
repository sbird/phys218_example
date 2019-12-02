import numpy as np
import halo_mass_function as HMF
import scipy as sp

#M_sun = 2 * pow(10, 30)
m = np.logspace(12, 20)
halo = HMF.HaloMassFunction(redshift=0)
h = halo.overden.hubble0
M = m/h

dndm_0 = halo.dndm_z(M, zz=0)
N = sp.integrate.trapz(dndm_0, M)

print("N = ",N)