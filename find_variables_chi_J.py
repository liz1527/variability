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
tbdata = fits.open('mag_flux_tables/J/mag_flux_table_best_J_extra_clean.fits')[1].data
chandata = fits.open('mag_flux_tables/J/xray_mag_flux_table_best_J_extra_clean.fits')[1].data
sdata = fits.open('mag_flux_tables/J/stars_mag_flux_table_J_extra_clean.fits')[1].data
sigtb = Table.read('sigma_tables/quad_epoch_sigma_table_J_extra_clean_2arcsec_neg.fits')

### Remove edges ###
tbdata = vari_funcs.field_funcs.remove_edges(tbdata)
chandata = vari_funcs.field_funcs.remove_edges(chandata)
sdata = vari_funcs.field_funcs.remove_edges(sdata)

## Create arrays of flux values ###
flux = vari_funcs.j_mag_flux.flux_stacks(tbdata, aper=4)
fluxchan = vari_funcs.j_mag_flux.flux_stacks(chandata, aper=4) 
sflux = vari_funcs.j_mag_flux.flux_stacks(sdata, aper=4)
#
### remove values that are negative ###
#flux, tbdata = vari_funcs.flux_funcs.noneg(flux, tbdata)
#fluxchan, chandata = vari_funcs.flux_funcs.noneg(fluxchan, chandata)
#sflux, sdata = vari_funcs.flux_funcs.noneg(sflux, sdata)

### Get error arrays ###
flux, fluxerr, tbdata = vari_funcs.j_mag_flux.create_quad_error_array_J(sigtb, tbdata, aper=4)
fluxchan, chanerr, chandata = vari_funcs.j_mag_flux.create_quad_error_array_J(sigtb, chandata, aper=4)
sflux, serr, sdata = vari_funcs.j_mag_flux.create_quad_error_array_J(sigtb, sdata, aper=4)

### Check chisq plot looks correct ###
fig, chisq = vari_funcs.selection_plot_funcs.flux_variability_plot(flux, fluxchan, 'chisq', 
                                       fluxerr=fluxerr, chanerr=chanerr,
                                       starflux=sflux, starfluxerr=serr,
                                       normalised=True, stars=True, scale='log')

plt.ylim(3e-2,3e4)
plt.xlim(4e-1, 1e7)

#### Calculate chi^2 values ###
#chisq = vari_funcs.vary_stats.my_chisquare_err(flux, fluxerr)
#chanchisq = vari_funcs.vary_stats.my_chisquare_err(fluxchan, chanerr)

### Select Variables as those with chisq > 24.322 and >50 ###
varydata32 = tbdata[chisq>32.08]
varydata30 = tbdata[chisq>30]
varydata40 = tbdata[chisq>40]
varydata50 = tbdata[chisq>50]

plt.hlines(32.08, 4e-1, 1e7,zorder=4,label='Chi>32.08')
plt.hlines(30, 4e-1, 1e7,'g', zorder=4,label='Chi>30')
plt.hlines(40, 4e-1, 1e7,'y', zorder=4,label='Chi>40')
#plt.hlines(50, 4e-1, 1e7,'c', zorder=4,label='Chi>50')
plt.legend()

#### Save new tables ###
save32 = Table(varydata32)
save32.write('variable_tables/J/K_extraction/J_variables_chi32_neg.fits')
#save30 = Table(varydata30)
#save30.write('variable_tables/variables_chi30.fits')
#save40 = Table(varydata40)
#save40.write('variable_tables/J/K_extraction/J_variables_chi40.fits')
#save50 = Table(varydata50)
#save50.write('variable_tables/variables_chi50.fits')

#### Isolate top left (RA<34.45, Dec>-5.1) ###
#mask1 = varydata24['ALPHA_J2000_05B'] < 34.45
#mask2 = varydata24['DELTA_J2000_05B'] > -5.1
#mask = mask1*mask2.astype('bool')
#tlvary = varydata24[mask]
#
#mask1 = tbdata['ALPHA_J2000_05B'] < 34.45
#mask2 = tbdata['DELTA_J2000_05B'] > -5.1
#mask = mask1*mask2.astype('bool')
#
#tlflux = flux[mask]
#tlmeanflux = np.nanmean(tlflux, axis=1)
#tlchi = chisq[mask]
#otherflux = flux[~mask]
#othermeanflux = np.nanmean(otherflux, axis=1)
#otherchi = chisq[~mask]
#
##plt.plot(othermeanflux, otherchi, 'mo')
#plt.plot(tlmeanflux, tlchi, 'yo')
#
#plt.figure()
#plt.hist([otherchi, tlchi], bins=np.logspace(-2,5,50), histtype='step', 
#         normed=True, label=['Other quadrants','Top left quadrant'])
#plt.xscale('log')
#plt.xlabel(r'$\chi^{2}$')
#plt.ylabel('Normalised Counts')
#plt.legend()
#plt.tight_layout()
#
#
#fig = vari_funcs.flux_variability_plot(flux, fluxchan, 'excess', 
#                                       fluxerr=fluxerr, chanerr=chanerr,
#                                       starflux=sflux, starfluxerr=serr,
#                                       normalised=True, stars=True, scale='symlog')







