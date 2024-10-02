#HW1, problem 1.3.a
#Writer: Aryana Haghjoo
#Oct 1, 2024

import pint
ureg = pint.UnitRegistry()
G = 6.674 * ureg.m**3 / (ureg.kg * ureg.s**2)  # Newtonian Gravitational constant in SI system
M = 1.989 * 10**30 * ureg.kg  # Mass of the Sun in SI system
c = 3.00 * 10**8 * ureg.m / ureg.s  # Speed of light in SI system

# Schwarzschild radius
r_sch = (2 * G * M / c**2).to(ureg.m)  # Calculate the Schwartzchild radius and ensure conversion to meters

# Print the result in scientific notation
print(f"The Schwarzschild radius of the Sun is {r_sch:.2e}.")