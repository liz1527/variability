#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 11:35:11 2018

Code to isolate variable sources using the bootstraped error bars and chi^2

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

### Open the fits files and get data ###
tbdata = fits.open('mag_flux_tables/mag_flux_table_best_extra_clean_no06.fits')[1].data
chandata = fits.open('mag_flux_tables/xray_mag_flux_table_best_extra_clean_no06.fits')[1].data
sdata = fits.open('mag_flux_tables/stars_mag_flux_table_extra_clean_no06.fits')[1].data
sigtb = Table.read('sigma_tables/quad_epoch_sigma_table_extra_clean_no06.fits')

## Create arrays of flux values ###
flux = vari_funcs.flux5_stacks(tbdata)
fluxchan = vari_funcs.flux5_stacks(chandata) 
sflux = vari_funcs.flux5_stacks(sdata)

### remove values that are negative ###
flux, tbdata = vari_funcs.noneg(flux, tbdata)
fluxchan, chandata = vari_funcs.noneg(fluxchan, chandata)
sflux, sdata = vari_funcs.noneg(sflux, sdata)

### Get error arrays ###
flux, fluxerr, tbdata = vari_funcs.create_quad_error_array(sigtb, tbdata)
fluxchan, chanerr, chandata = vari_funcs.create_quad_error_array(sigtb, chandata)
sflux, serr, sdata = vari_funcs.create_quad_error_array(sigtb, sdata)

### Check chisq plot looks correct ###
fig = vari_funcs.flux_variability_plot(flux, fluxchan, 'chisq', 
                                       fluxerr=fluxerr, chanerr=chanerr,
                                       starflux=sflux, starfluxerr=serr,
                                       normalised=True, stars=True, scale='log')

### Calculate chi^2 values ###
chisq = vari_funcs.my_chisquare_err(flux, fluxerr)
chanchisq = vari_funcs.my_chisquare_err(fluxchan, chanerr)

### Select Variables as those with chisq > 24.322 and >50 ###
varydata24 = tbdata[chisq>24.322]
varydata30 = tbdata[chisq>30]
varydata40 = tbdata[chisq>40]
varydata50 = tbdata[chisq>50]

plt.hlines(24.322, 2e2, 1e7,zorder=4,label='Chi>24.3')
plt.hlines(30, 2e2, 1e7,'g', zorder=4,label='Chi>30')
plt.hlines(40, 2e2, 1e7,'y', zorder=4,label='Chi>40')
plt.hlines(50, 2e2, 1e7,'c', zorder=4,label='Chi>50')
plt.legend()

#### Save new tables ###
#save24 = Table(varydata24)
#save24.write('variable_tables/variables_chi24.fits')
#save30 = Table(varydata30)
#save30.write('variable_tables/variables_chi30.fits')
#save40 = Table(varydata40)
#save40.write('variable_tables/variables_chi40.fits')
#save50 = Table(varydata50)
#save50.write('variable_tables/variables_chi50.fits')

### Isolate top left (RA<34.45, Dec>-5.1) ###
mask1 = varydata24['ALPHA_J2000_05B'] < 34.45
mask2 = varydata24['DELTA_J2000_05B'] > -5.1
mask = mask1*mask2.astype('bool')
tlvary = varydata24[mask]

mask1 = tbdata['ALPHA_J2000_05B'] < 34.45
mask2 = tbdata['DELTA_J2000_05B'] > -5.1
mask = mask1*mask2.astype('bool')

tlflux = flux[mask]
tlmeanflux = np.nanmean(tlflux, axis=1)
tlchi = chisq[mask]
otherflux = flux[~mask]
othermeanflux = np.nanmean(otherflux, axis=1)
otherchi = chisq[~mask]

#plt.plot(othermeanflux, otherchi, 'mo')
plt.plot(tlmeanflux, tlchi, 'yo')

plt.figure()
plt.hist([otherchi, tlchi], bins=np.logspace(-2,5,50), histtype='step', 
         normed=True, label=['Other quadrants','Top left quadrant'])
plt.xscale('log')
plt.xlabel(r'$\chi^{2}$')
plt.ylabel('Normalised Counts')
plt.legend()
plt.tight_layout()


fig = vari_funcs.flux_variability_plot(flux, fluxchan, 'excess', 
                                       fluxerr=fluxerr, chanerr=chanerr,
                                       starflux=sflux, starfluxerr=serr,
                                       normalised=True, stars=True, scale='symlog')







