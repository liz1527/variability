#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 10:58:14 2019

Code to check the Chi-squared of variables after most deviant point is removed

@author: ppxee
"""

### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
from astropy.table import Table #for handling tables
import numpy as np #for handling arrays
#import math
#from astropy.stats import median_absolute_deviation
from scipy import stats
import vari_funcs #my module to help run code neatly
#from scipy.stats import chisquare
plt.close('all') #close any open plots

### Read fits tables ###
varydata = fits.open('variable_tables/J_variables_chi40_noneg.fits')[1].data
sigtb = Table.read('sigma_tables/quad_epoch_sigma_table_extra_clean_2arcsec_J_noneg.fits')

### remove edges ###
varydata = vari_funcs.field_funcs.remove_edges(varydata)

### Create flux and error arrays ###
varyflux = vari_funcs.j_mag_flux.flux4_stacks(varydata)
varyflux, varyfluxerr, varydata = vari_funcs.j_mag_flux.create_quad_error_array_J(sigtb, varydata, aper=4)
meanflux = np.nanmean(varyflux, axis=1)

### Find difference from 1 ###
fluxn, fluxnerr = vari_funcs.flux_funcs.normalise_flux_and_errors(varyflux, varyfluxerr)
diff = abs(fluxn - 1)

### mask max difference ###
mask = np.argmax(diff, axis=1)
varyfluxnew = np.copy(varyflux)
varyfluxerrnew = np.copy(varyfluxerr)
varyfluxnew[np.arange(len(varyflux)),mask] = np.nan
varyfluxerrnew[np.arange(len(varyflux)),mask] = np.nan


### Calculate new chi^2 values ###
chisq = vari_funcs.vary_stats.my_chisquare_err(varyflux, varyfluxerr)
chisqnew = vari_funcs.vary_stats.my_chisquare_err(varyfluxnew, varyfluxerrnew)

#
### remove those drop below the 99.9% line ###
mask2 = chisqnew < 16.68
cutchisq = chisq[~mask2]
cutmeanflux = meanflux[~mask2]
cutchisqnew = chisqnew[~mask2]

#### Plot hist of differences ###
#plt.figure()
#plt.hist(diff, bins=np.logspace(0,4))
#plt.xscale('log')

### plot chisq v flux for both ###
plt.figure(figsize=[8,8])
plt.plot(meanflux, chisq, 'o')
plt.plot(meanflux, chisqnew, 'o')
for i in range(0, len(meanflux)):
    plt.plot([meanflux[i],meanflux[i]], [chisq[i],chisqnew[i]], 'r-')
for i in range(0, len(cutmeanflux)):
    plt.plot([cutmeanflux[i],cutmeanflux[i]], [cutchisq[i],cutchisqnew[i]], 'k-')
#plt.plot(cutmeanflux, cutchisq, 'rs', mfc='None', markersize=10)
#plt.plot(cutmeanflux, cutchisqnew, 'ks', mfc='None')
plt.hlines(30, 1e2,1e6, color='C0', label='Chi=30')
plt.hlines(16.68, 1e2,1e6, color='C1', label='Chi=16.68')

plt.xscale('log')
plt.yscale('log')

#### Plot hist of chi with p value marked ###
#plt.figure()
#plt.hist(chisqnew, bins=np.logspace(0,4))
#plt.vlines(27.83, 0, 50)
#plt.xscale('log')
#plt.xlabel(r'$\chi^{2}$')
#plt.ylabel('Number')
#
newvary = varydata[chisqnew<16.68] #99.9%
#
### plot positions ###
plt.figure(figsize=[8,7])
plt.scatter(varydata['X_IMAGE_05B'], varydata['Y_IMAGE_05B'])
plt.scatter(newvary['X_IMAGE_05B'], newvary['Y_IMAGE_05B'])
plt.xlabel('X')
plt.ylabel('Y')
plt.tight_layout()

### save checked table ###
#save30 = Table(newvary)
#save30.write('variable_tables/J_variables_chi40_noneg_deviant.fits')