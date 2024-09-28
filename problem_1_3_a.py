import pint

# Create a unit registry
ureg = pint.UnitRegistry()

# Constants
G = 6.674 * ureg.m**3 / (ureg.kg * ureg.s**2)  # Gravitational constant
M = 1.989 * 10**30 * ureg.kg  # Mass of the Sun
c = 3.00 * 10**8 * ureg.m / ureg.s  # Speed of light

# Calculate the Schwarzschild radius
r_s = (2 * G * M / c**2).to(ureg.m)  # Convert to meters

# Print the result
print(f"The Schwarzschild radius of the Sun is {r_s:.2f}.")