import pint

def calculate_schwarzchild_radius(mass):
    ureg = pint.UnitRegistry()
    ureg.define("Msun = 1.98855*10**30 * kg")

    schwarzchild_radius = 2 * ureg.newtonian_constant_of_gravitation * ureg.Msun * mass / ureg.c**2
    return schwarzchild_radius.to_base_units()

if __name__ == "__main__":
    mass = 1 # in solar mass unit
    print(calculate_schwarzchild_radius(mass))