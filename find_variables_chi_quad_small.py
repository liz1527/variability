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
fig,_ = vari_funcs.flux_variability_plot(flux, fluxchan, 'chisq', 
                                       fluxerr=fluxerr, chanerr=chanerr,
                                       starflux=sflux, starfluxerr=serr,
                                       #normalised=True, 
                                       stars=True, scale='log')
fig.canvas.mpl_connect('pick_event', vari_funcs.onpickflux)

### Calculate chi^2 values ###
chisq = vari_funcs.my_chisquare_err(flux, fluxerr)
chanchisq = vari_funcs.my_chisquare_err(fluxchan, chanerr)

### Select Variables as those with chisq > 22.458 and >50 ###
varydata24 = tbdata[chisq>22.458]
varydata30 = tbdata[chisq>30]
varydata40 = tbdata[chisq>40]
varydata50 = tbdata[chisq>50]

#plt.hlines(22.458, 8e1, 1e7,zorder=4,label='Chi>22.5')
plt.hlines(30, 8e1, 1e7,'g', zorder=4,label='Chi=30')
#plt.hlines(40, 8e1, 1e7,'y', zorder=4,label='Chi>40')
#plt.hlines(50, 8e1, 1e7,'c', zorder=4,label='Chi>50')
plt.legend()
plt.tight_layout()

#### Check if variables vanish if the most deviant semester is removed ###
#varyflux = vari_funcs.flux5_stacks(varydata40)
#varyflux, varyfluxerr, varydata40 = vari_funcs.create_quad_error_array(sigtb, varydata40)
#fluxn, fluxnerr = vari_funcs.normalise_flux_and_errors(varyflux, varyfluxerr)
#diff = abs(fluxn - 1)
#mask = np.argmax(diff, axis=1)
#varyfluxnew = np.copy(varyflux)
#varyfluxerrnew = np.copy(varyfluxerr)
#varyfluxnew[np.arange(len(varyflux)),mask] = np.nan
#varyfluxerrnew[np.arange(len(varyflux)),mask] = np.nan
#
#varychi = vari_funcs.my_chisquare_err(varyflux, varyfluxerr)
#varychinew = vari_funcs.my_chisquare_err(varyfluxnew, varyfluxerrnew)
#chidiff = varychi - varychinew
### Save new tables ###
#save24 = Table(varydata24)
#save24.write('variable_tables/no06_variables_chi22.fits')
#save30 = Table(varydata30)
#save30.write('variable_tables/no06_variables_chi30_2arcsec.fits')
#save40 = Table(varydata40)
#save40.write('variable_tables/no06_variables_chi40.fits')
#save50 = Table(varydata50)
#save50.write('variable_tables/no06_variables_chi50.fits')
