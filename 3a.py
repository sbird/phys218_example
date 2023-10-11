'''
Problem 1.3a in PHYS 206 to calculate the Schwarzschild radius of the sun in metres
'''
import pint


u = pint.UnitRegistry()


G = (6.67e-11 * u.meter**3 / u.second**2 / u.kilogram)  # Gravitational constant in SI units
c = 3e8 * u.meter / u.second  # Speed of light in SI units


def get_sch_rad(m):
    """
    This function returns the Schwarzschlid radius of a black hole in metres

    Args:
    m: Mass of the black hole in Msun

    Returns:
    r: Schwarzchild radius for a given mass of the black hole in metres
    """

    r = 2 * G * m / c**2
    r = r.to(u.meter)  # This is to confirm if the final units are in metres

    return r


msun = 2e30 * u.kilogram  # mass of the sun

print(f"Schwarzschild radius of the sun is {round(get_sch_rad(msun), 2)}")
