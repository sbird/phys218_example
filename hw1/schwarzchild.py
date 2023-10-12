# Calculate the Schwarzchild radius in meters using pint


import pint

mass = float(input("Enter mass in kg:"))

ureg = pint.UnitRegistry()
def Schwarzschild(m):
    """ Calculates Schwarzchild radius in Meters. Inpute mass in kg"""
    m = m * ureg.gram * 1000 
    G = 6.674e-11  * (ureg.meter)**3/(ureg.gram*1000 * (ureg.second**2)) 
    c = 3e8 * ureg.meter/ureg.second
    r = 2*G*m/(c**2)
    # assert r.units==ureg.meter
    assert r.magnitude>0
    
    print(r)
    print(r.to(ureg.km))
Schwarzschild(mass)