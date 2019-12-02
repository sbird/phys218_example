""" Using halo_mass_function.py, write a function to compute the total number of halos at z=0 with a mass above 10^12 Msolar
according to one of the listed halo mass function formulae.
I formatted this as a script to easily run it from the command line, but the function works fine on its own.
The result: the number of halos with mass > 10^12 Msolar at z=0 is effectively 0. (as close as python can get) """

# The plan: get the mass function dn/dM and integrate for masses above 10^12 Msolar

import numpy as np
import scipy
import halo_mass_function as hm

def num_halos_above(z, msolar):
    halo = hm.HaloMassFunction(redshift=z)
    # halo.dndm requires mass in units of M_sun / h, so we need to divide the mass by hubble:
    h = halo.overden.hubble0
    mass = msolar / h
    dndM = lambda m: halo.dndm(m)
    nhalos = scipy.integrate.quad(dndM, mass, np.inf)
    # scipy integrating returns results for the upper and lower bounds so let's properly subtract them:
    return nhalos[0]-nhalos[1]

if __name__ == "__main__":
    z = 0
    msolar = 10**12
    nhalos = num_halos_above(z, msolar)
    print(nhalos)
