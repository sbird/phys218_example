from pint import UnitRegistry
import numpy as np
ureg = UnitRegistry()


#defining variables 
G = ureg.newtonian_constant_of_gravitation
M_sun = 1.989e30 * ureg.kilogram
c = ureg.speed_of_light

#define a function
def schwarzs_rad(m):
    r = 2*G*m/c**2
    return r.to_base_units()

print("The Schwarzschild_radius of the sun is", schwarzs_rad(M_sun))
