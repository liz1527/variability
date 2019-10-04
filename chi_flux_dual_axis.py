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
sigtb = Table.read('sigma_tables/quad_epoch_sigma_table_extra_clean_no06_2arcsec.fits')

### Remove edges ###
tbdata = vari_funcs.remove_edges(tbdata)
chandata = vari_funcs.remove_edges(chandata)
sdata = vari_funcs.remove_edges(sdata)

## Create arrays of flux values ###
flux = vari_funcs.flux4_stacks(tbdata)
fluxchan = vari_funcs.flux4_stacks(chandata) 
sflux = vari_funcs.flux4_stacks(sdata)

### remove values that are negative ###
flux, tbdata = vari_funcs.noneg(flux, tbdata)
fluxchan, chandata = vari_funcs.noneg(fluxchan, chandata)
sflux, sdata = vari_funcs.noneg(sflux, sdata)

### Get error arrays ###
flux, fluxerr, tbdata = vari_funcs.create_quad_error_array(sigtb, tbdata, aper=4)
fluxchan, chanerr, chandata = vari_funcs.create_quad_error_array(sigtb, chandata,aper=4)
sflux, serr, sdata = vari_funcs.create_quad_error_array(sigtb, sdata, aper=4)

### reset X-ray column as messed up by stacking ###
tbdata['X-ray'][tbdata['X-ray']==70] = False 
tbdata['X-ray'][tbdata['X-ray']==84] = True

### Check chisq plot looks correct ###
#fig,_ = vari_funcs.flux_variability_plot(flux, fluxchan, 'chisq', 
#                                       fluxerr=fluxerr, chanerr=chanerr,
#                                       starflux=sflux, starfluxerr=serr,
#                                       #normalised=True, 
#                                       stars=True, scale='log')
fig = plt.figure(figsize=[8,8])
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()

avgfluxperob = np.nanmean(flux, axis=1) #for UDS
avgfluxchanperob = np.nanmean(fluxchan, axis=1) #for non-stellar chandra
savgfluxperob = np.nanmean(sflux, axis=1) #for stars

### plot chi squared ###
vary = vari_funcs.my_chisquare_err(flux, fluxerr)
varychan = vari_funcs.my_chisquare_err(fluxchan, chanerr)
varystar = vari_funcs.my_chisquare_err(sflux, serr)

ax1.set_ylabel('Chi Squared')
### Plot the variability v mean as appropriate ###
plt.plot(savgfluxperob, varystar, 'm*', mfc = 'none', markersize = 10,
             label='DR11 Star') 
line, = plt.plot(avgfluxperob, vary, 'b+', label='Galaxy', picker=2)
plt.plot(avgfluxchanperob, varychan, 'ro', mfc = 'none', markersize = 10,
         label='X-ray detected') #no picker as will be selected in the UDS point

    
### Apply required plot charateristics ###
ax1.set_xscale('log')
ax1.set_yscale('log')
    
ax1.set_ylim(3e-2,3e4)
ax1.set_xlim(8e1, 1e7)
#    plt.xlim(13,26)
ax1.set_xlabel('Mean Flux')

ticks = np.array([1e2, 1e3, 1e4, 1e5, 1e6, 1e7])
new_ticks = 30 - 2.5*np.log10(ticks)

ax2.set_xlim(ax1.get_xlim())
ax2.set_xscale('log')
ax2.set_xticks(ticks)
ax2.set_xticklabels(new_ticks)
ax2.minorticks_off()
ax1.minorticks_on()
ax2.set_xlabel('Mean $K$-band Magnitude (AB)')


plt.legend()

fig.canvas.mpl_connect('pick_event', vari_funcs.onpickflux)

### Calculate chi^2 values ###
chisq = vari_funcs.my_chisquare_err(flux, fluxerr)
chanchisq = vari_funcs.my_chisquare_err(fluxchan, chanerr)

### Select Variables as those with chisq > 22.458 and >50 ###
varydata24 = tbdata[chisq>22.458]
varydata30 = tbdata[chisq>30]
varydata40 = tbdata[chisq>40]
varydata50 = tbdata[chisq>50]

plt.hlines(30, 8e1, 1e7,'g', zorder=4,label='Chi=30')

plt.legend()
plt.tight_layout()
