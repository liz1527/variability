#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 20 17:40:28 2018

@author: ppxee
"""


### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
from astropy.table import Table #for handling tables
import numpy as np #for handling arrays
#import math
from astropy.stats import median_absolute_deviation
import vari_funcs #my module to help run code neatly
#from matplotlib.mlab import griddata
from scipy.interpolate import griddata
plt.close('all') #close any open plots

### Open the fits files and get data ###
combined = fits.open('mag_flux_tables/mag_flux_table_best.fits')
tbdata = combined[1].data
chandra = fits.open('mag_flux_tables/xray_mag_flux_table_best.fits')
chandata = chandra[1].data
stars = fits.open('mag_flux_tables/stars_mag_flux_table.fits')
sdata = stars[1].data

### Restrict objects to those in the Chandra field ###
#tbdata = vari_funcs.chandra_only(tbdata)

## Create arrays of flux values but without 06B ###
flux = vari_funcs.flux5_stacks(tbdata)
fluxchan = vari_funcs.flux5_stacks(chandata) # for chandra non-stellar objects
sflux = vari_funcs.flux5_stacks(sdata)


### remove values that are +/-99 ###
flux, tbdata = vari_funcs.noneg(flux, tbdata)
fluxchan, chandata = vari_funcs.noneg(fluxchan, chandata)
sflux, sdata = vari_funcs.noneg(sflux, sdata)

### Normalise arrays ###
#fluxn = vari_funcs.normalise_flux(flux)
#fluxchann = vari_funcs.normalise_flux(fluxchan)
#sfluxn = vari_funcs.normalise_flux(sflux)
#
#
#### Multiply all flux values in a yearstack by the correct constant ###
#fluxcorrn = vari_funcs.psf_correct(fluxn, fluxn, 'median') 
#fluxchancorrn = vari_funcs.psf_correct(fluxn, fluxchann, 'median') 
#fluxcorr = vari_funcs.psf_correct(flux, flux, 'median') 
#fluxchancorr = vari_funcs.psf_correct(flux, fluxchan, 'median') 

fig,_ = vari_funcs.flux_variability_plot(flux, fluxchan, 'mad', starflux=sflux,
                                      stars=True, normalised=True)

fig.canvas.mpl_connect('pick_event', vari_funcs.onpick)
# Find outliers ###
bins = np.array([13, 15])
bins = np.append(bins, np.arange(16,24,0.2))
bins = np.append(bins, [24, 25, 26])
bins = 10**((30-bins)/2.5)
bins = np.flip(bins, axis=0)

#outliers, tb, modz = vari_funcs.find_outliers(fluxn, tbdata, bins)
#varys = tb[outliers]
#varyflux = vari_funcs.mag5_stacks(varys)
#varymad = median_absolute_deviation(varyflux, axis=1) 
#varymean = np.mean(varyflux, axis=1)
#plt.figure(1)
#plt.plot(varymean, varymad, 'md', mfc='none', markersize=10)


outliers2, tb2, modz2 = vari_funcs.find_outliers(flux, tbdata, bins, threshold=5)
tb2['X-ray'][tb2['X-ray']==70] = False 
tb2['X-ray'][tb2['X-ray']==84] = True
### Find plotting values for the new flux table ###
flux2 = vari_funcs.flux5_stacks(tb2) 
flux2, tb2 = vari_funcs.noneg(flux2, tb2)
avgfluxperob2 = np.nanmean(flux2, axis=1) #for UDS
flux2 = vari_funcs.normalise_flux(flux2)
vary2 = median_absolute_deviation(flux2, axis=1)


### create table of varying ###
magmask = avgfluxperob2 < 10**((30-21)/2.5)
outliers2 = outliers2*magmask
varydata = tb2[outliers2]
#cols = fits.ColDefs(varydata)
#hdu = fits.BinTableHDU.from_columns(cols)
#hdu.writeto('variable_tables/variable_mag_flux_table_with06B_mag21.fits')

varyfluxcorr = flux2[outliers2]
varymadcorr = vary2[outliers2]
varymeancorr = avgfluxperob2[outliers2]
varymodz = modz2[outliers2]
#plt.plot(varymeancorr, varymadcorr, 'kd', mfc='none', markersize=10)

#%% Find the ones that used to be variable but not aren't
#oldvary = fits.open('variable_mag_flux_table.fits')[1].data
#nums = ~np.isin(oldvary['NUMBER_05B'], varydata['NUMBER_05B'])
#oldvary = oldvary[nums]
#ind = np.isin(tbdata['NUMBER_05B'], oldvary['NUMBER_05B'])
#oldvary = tbdata[ind]

#varyflux = vari_funcs.mag5_stacks(oldvary)
#varymean = np.mean(varyflux, axis=1)
#varymad = median_absolute_deviation(varyflux, axis=1)
#plt.plot(varymean, varymad, 'gd', mfc='none', markersize=10)

#plt.figure()
#plt.scatter(varymeancorr, varymadcorr, c=varymodz, marker='+')
#plt.colorbar()

## Plot mod z in scatter ###
plt.figure()
plt.scatter(avgfluxperob2, vary2, c=modz2, marker='+', vmax=6)
cbar = plt.colorbar()
cbar.set_label('Modified z-score')
plt.yscale('log')
plt.xscale('log')
plt.ylim(1e-4, 1e1)
plt.xlabel('Mean Magnitude')
plt.ylabel('MAD Magnitude')

#plt.plot(varymeancorr, varymadcorr, 'md', mfc='none', markersize=10)

### plot mod z score as contours ###
xi = np.logspace(2, 7, 1000)
yi = np.logspace(-4, 1, 1000)

xi, yi = np.meshgrid(xi,yi)
#plt.figure()
#plt.plot(xi,yi,'o')
#plt.xscale('log')
#plt.yscale('log')

#zi = griddata(avgfluxperob2, vary2, modz2, xi, yi)
zi = griddata((avgfluxperob2, vary2), modz2, (xi, yi))
#plt.plot(avgfluxperob, vary, 'b+')
plt.figure(1)
plt.contour(xi, yi, zi, [0,4,5,6,8,10], zorder=3, linewidths=2.5)
cbar = plt.colorbar()
cbar.set_label('Modified z-score')
ax = fig.add_subplot(111)
for axis in ['top','bottom','left','right']:
  ax.spines[axis].set_linewidth(1.5)
#interest = tb2['NUMBER_05B'] == 1
#xint = avgfluxperob2[interest]
#yint = vary2[interest]
#plt.plot(xint, yint, 'ks', markersize=10)
  

  
  
#  
#  
#  ########## code that worked for magnitudes
#  
#  ### plot mod z score as contours ###
#xi = np.linspace(min(avgfluxperob2), max(avgfluxperob2), 1000)
#yi = np.logspace(-4, 1, 5000)
#
#
#zi = griddata(avgfluxperob2, vary2, modz2, xi, yi, interp='linear')
##plt.plot(avgfluxperob, vary, 'b+')
#plt.figure(1)
#plt.contour(xi, yi, zi, [0,4,5,6,8,10], zorder=3, linewidths=2.5)
#cbar = plt.colorbar()
#cbar.set_label('Modified z-score')
#ax = fig.add_subplot(111)
#for axis in ['top','bottom','left','right']:
#  ax.spines[axis].set_linewidth(1.5)