import pint
from pint import UnitRegistry

class sun:
    def schwarzchild(mass):
        ureg = pint.UnitRegistry()
        G = 1.0 * ureg.newtonian_constant_of_gravitation
        c = 1.0 * ureg.speed_of_light
        mass = mass * ureg.kilogram
        G_si = G.to_base_units()
        c_si = c.to_base_units()
        sch_r = G_si*mass/pow(c_si,2)
        return sch_r

    mass = 2 * pow(10,30)
    r = schwarzchild(mass)
    print(r)






