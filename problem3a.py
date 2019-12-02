""" find the Schwarzschild radius of the Sun in m using pint"""

import pint

class Sun:
    """ Class to describe a star based on its mass in terms of solar masses """
    def __init__(self, mass):
        self.ureg = pint.UnitRegistry()
        self.ureg.define("Msolar = 1.98855*10**30 * kilogram")
        self.mass = mass * self.ureg.Msolar

    def schwarz(self):
        """ Find the Schwarzchild radius for the class """
        g_newt = self.ureg.newtonian_constant_of_gravitation
        msun = self.mass
        r_sch = 2 * g_newt * msun / self.ureg.speed_of_light**2
        return r_sch.to_base_units()

def schwarz_rad(mass):
    """ Given a mass, find the Schwarzschild radius """
    star = Sun(mass)
    radius = star.schwarz()
    return radius

if __name__ == "__main__":
    MASS = 1.0
    RAD = schwarz_rad(MASS)
    print(RAD)
