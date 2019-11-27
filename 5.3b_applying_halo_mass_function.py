#!/usr/bin/env python

from halo_mass_function import HaloMassFunction,Overdensities
import numpy as np
import matplotlib.pyplot as plt

save_fig=True #save the figure or not?
show_fig=False #show the figure or not?
fontsize=15

#Units and constants
M_array = np.linspace(1e12,1e20,1000)
dM = (1e20-1e12)/1000
z=0

#Calculation of number of halos
N_array = HaloMassFunction(redshift=z).dndm(mass=M_array)#*M_array
N_total = np.sum(N_array*dM)

#total number of halos above 1.e12
print("The total number of halos above 10**12 solar masses: {:.3f} halos".format(N_total))

# #Masking out non-positive number so I can take logarithm
# mask = N_halos==0
# N_halos[mask] = np.nan

# #Plotting routine
# fig,ax = plt.subplots(1)
# ax.plot(np.log10(M_array),np.log10(N_halos))
# ax.set_xlabel(r"$log_{10}\left[M_{\odot}\right]$",fontsize=fontsize)
# ax.set_ylabel(r"$log_{10}\left[N_{halo}\right]$",fontsize=fontsize)

# #Saving and or showing plot
# if save_fig:
# 	plt.savefig(fname="Number-of-halos_vs_halo-mass.pdf")
# if show_fig:
# 	plt.show()