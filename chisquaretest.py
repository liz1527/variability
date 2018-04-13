#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  7 10:44:58 2018

@author: ppxee
"""

### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
#from astropy.table import Table #for handling tables
import numpy as np #for handling arrays
#import math
from astropy.stats import median_absolute_deviation
from scipy.stats import chisquare
from matplotlib.mlab import griddata

import vari_funcs #my module to help run code neatly
plt.close('all') #close any open plots

### Open the fits files and get data ###
combined = fits.open('mag_flux_table_best.fits')
tbdata = combined[1].data
chandra = fits.open('xray_mag_flux_table_best.fits')
chandata = chandra[1].data
stars = fits.open('stars_mag_flux_table.fits')
sdata = stars[1].data

## Create arrays of flux values ###
fluxn = vari_funcs.mag5_stacks(tbdata)
fluxchann = vari_funcs.mag5_stacks(chandata) # for chandra non-stellar objects
sfluxn = vari_funcs.mag5_stacks(sdata)

### remove values that are +/-99 ###
flux, tbdata = vari_funcs.no99(fluxn, tbdata)
fluxchann, chandata = vari_funcs.no99(fluxchann, chandata)
sfluxn, sdata = vari_funcs.no99(sfluxn, sdata)

## Create arrays of flux values ###
fluxn = vari_funcs.flux5_stacks(tbdata)
fluxchann = vari_funcs.flux5_stacks(chandata) # for chandra non-stellar objects
sfluxn = vari_funcs.flux5_stacks(sdata)

avgflux = np.mean(fluxn, axis=1)
expected = np.column_stack([avgflux, avgflux, avgflux, avgflux, avgflux, avgflux, avgflux, avgflux])#np.ones(np.shape(flux))
### normalise ###
#fluxn = vari_funcs.normalise_mag(flux)
[chisqn, pn] = chisquare(fluxn, expected, axis=1)
[chisq, p] = chisquare(flux, axis=1)

plt.hist(pn, label='normalised')
plt.legend()

#plt.figure()
#plt.plot(avgflux, chisq, 'b+')
#plt.yscale('log')

#fig = vari_funcs.flux_variability_plot(flux, fluxchann, 'chisq', 
#                                       normalised=True, stars=True,
#                                       starflux = sfluxn)
#
### plot pvalue as contours ###
#xi = np.linspace(min(avgflux), max(avgflux), 1000)
#yi = np.logspace(-4, 1, 5000)
#
#
#zi = griddata(avgflux, chisqn, pn, xi, yi, interp='linear')
#plt.contour(xi, yi, zi)#, [0,4,5,6,8,10], zorder=3)
#cbar = plt.colorbar()
#cbar.set_label('p-value')