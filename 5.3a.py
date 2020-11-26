from pint import UnitRegistry
ureg=UnitRegistry()
#speed of light
c=ureg.speed_of_light
C=(c*1).to(ureg.m/ureg.s)
#G
g=ureg.gravitational_constant
G=(g*1).to(ureg.m**3/(ureg.kg * ureg.s**2))
#Solar Mass
SM= 2e+30 * ureg.kg

rs=(2*SM*G)/(C**2)
print("The Schwarzchild radius of the Sun is:")
print(rs)
