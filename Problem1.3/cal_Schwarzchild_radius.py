"Problem1.3(a) Write a python script, using pint, which finds the Schwarzchild radius of the Sun, in m."
import numpy as np
import pint
from astropy import constants as const

units=pint.UnitRegistry()
def Cal_R_Sch(m_in_solarms):
    const_G=const.G.value*units.meter**3/(units.kilogram*units.second**2)
    const_solarmass=const.M_sun.value*units.kilogram
    const_sol=const.c.value*units.meter/units.second

    return 2*const_G*const_solarmass*m_in_solarms/const_sol**2

radius=Cal_R_Sch(1).magnitude
print("The Schwarzchild radius of the Sun is "+"%.2f"%radius+"m.")