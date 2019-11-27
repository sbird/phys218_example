#!/usr/bin/env python
"""Program to find the schwarzchild radius of the sun"""
import numpy as np
import pint

#Units and constants
ureg = pint.UnitRegistry()
ureg.define("Msolar = 1.98855*10**30 * kilogram")
M_sun = 1*ureg.Msolar
c = ureg.c
G = ureg.newtonian_constant_of_gravitation


def get_Rs(m):
	"""
	Computes the Schwarzchild Radius of an object given its mass in units on meters
	"""
	return (2*G*m/c**2).to(ureg.m)

Rs_sun = get_Rs(M_sun)
print("Schwarzchild Radius of the Sun: {0:.3f}s".format(Rs_sun))