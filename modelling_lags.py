#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 15:18:56 2020

Code to determine expected lags for my sample

@author: ppxee
"""

import time
start = time.time()
print(start)

### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
from matplotlib import cm
import matplotlib.colors as colors
from astropy.io import fits #for handling fits
from astropy.table import Table #for handling tables
import numpy as np #for handling arrays
#import math
#from astropy.stats import median_absolute_deviation
import vari_funcs #my module to help run code neatly
from astropy.constants import b_wien, G, M_sun, c
import astropy.units as u

def get_radius(wavelength, acc_rate, M_bh):
    term_1 = (b_wien/(6.3e5*u.K*wavelength*u.m))**(-4/3)
    term_2 = (acc_rate)**(1/3) # make sure acc rate is defined as % eddington
    term_3 = (M_bh*M_sun/(1e8*M_sun))**(-1/3) #make sure mass is in unit M_sun
    R_s = (2*G*M_bh*M_sun)/(c**2)#make sure mass is in unit M_sun
    
#    term_1 = (b_wien.value/(6.3e5*wavelength))**(-4/3)
#    term_2 = (acc_rate)**(1/3) # make sure acc rate is defined as % eddington
#    term_3 = (M_bh/(1e8))**(-1/3) #make sure mass is in unit M_sun
#    R_s = (2*G.value*M_bh)/(c.value**2)
    
    radius = term_1 * term_2 * term_3 * R_s
    
    return radius

def get_lag(acc_rate, M_bh, z):
    rad_j = get_radius(1.2e-6/(1+z), acc_rate, M_bh)
    rad_k = get_radius(2.2e-6/(1+z), acc_rate, M_bh)
    
    rad_diff = rad_k - rad_j
    
    t_light_s = rad_diff/c
    t_light_mon = t_light_s.value/(60*60*24*30) #min - hr - day - month
    
    return t_light_mon

t_lag = get_lag(0.1,1e9,1)
acc_rates = np.linspace(0.05,1,100)
bh_masses = np.logspace(6,10,100)

### create grid to put contours on ###
x, y = np.meshgrid(acc_rates, bh_masses)

### get lags at each point ###
lags = np.empty(np.shape(x))
for n, acc_rate in enumerate(acc_rates):
    for m, M_bh in enumerate(bh_masses):
        lags[m,n] = get_lag(acc_rate, M_bh, 1)
        
### Define contour levels ###
levels = np.array([0.25,0.5,0.75,1,1.25,1.5,1.75,2,3,4,5,6,7])

#%% Plot ###
fig = plt.figure(figsize=[9,6])
CS = plt.contour(x, y, lags, levels=levels, cmap=cm.viridis_r)
plt.clabel(CS, inline=1, fontsize=10)
im = plt.imshow(lags, origin='lower', extent=(0.05, 1, 1e6, 1e10),
               cmap='gray',norm=colors.LogNorm(vmin=lags.min(), vmax=lags.max()))
plt.xlabel('Accretion Rate ($\dot{M_{EDD}}$)')
plt.ylabel('Black Hole Mass ($M_\odot$)')
plt.yscale('log')
cbar = plt.colorbar(im)#, orientation='horizontal',  fraction=0.03)
cbar.set_label('Lag (months)')
cbar2 = plt.colorbar(CS)
cbar2.set_label('Lag (months)')
plt.tight_layout()

end = time.time()
print(end-start)




























