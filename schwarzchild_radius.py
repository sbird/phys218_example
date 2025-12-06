import pint


ureg = pint.UnitRegistry()

# Define constants
# https://physics.nist.gov/cgi-bin/cuu/Value?bg
G = 6.67430E-11 * ureg.meter**3 / ureg.kilogram / ureg.second**2 # Newton's constant
# https://web.archive.org/web/20140811195806/http://www.bipm.org/en/si/new_si/explicit_constant.html
c = 299792458 * ureg.meter / ureg.second # speed of light
# https://www.mccc.edu/~dornemam/Planet_Walk/Sun/the_sun.htm#:~:text=The%20sun%20has%20a%20mass%20of%201.9891x1030,it%20could%20hold%20over%201%20million%20Earths.
M_sun = 1.9891E30 * ureg.kilogram # solar mass


if __name__ == "__main__":
    # Kutner 2003, p. 148
    # https://archive.org/details/astronomyphysica00kutn/mode/2up
    print(2*G*M_sun/(c**2)) # Schwarzchild radius formula