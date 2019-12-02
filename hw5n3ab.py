import numpy as np
import pint

ureg = pint.UnitRegistry()

# (a)
ureg.define('Solar_Mass = 2e30 * kilogram = Msolar')

M = 1 * ureg.Msolar
G = 1 * ureg.newtonian_constant_of_gravitation
c = 1* ureg.speed_of_light

rsch = G.to_base_units() * M.to_base_units() / c.to_base_units()**2 / 2
rsch.to_base_units()
print("5.3 a) rsch = ",rsch)

# (b)

import halo_mass_function as hmf
import scipy as sp

halo_pop = hmf.HaloMassFunction(redshift=0)

hubble = halo_pop.overden.hubble0

mass = np.logspace(12, 18, 50)
mass = np.divide(mass, hubble) # convert to Solar mass / h

dn_dm = halo_pop.dndm(mass)
N = sp.integrate.trapz(dn_dm,mass)

print("5.3 b) There are {0} halos above 10**12 solar mass per Mpc/h**3".format(N))
