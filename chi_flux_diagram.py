#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 11:35:11 2018

Code to creat chi-flux plot using the bootstraped error bars for 3 arcsec K band

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
tbdata = fits.open('mag_flux_tables/K/mag_flux_table_best_extra_clean_no06.fits')[1].data
chandata = fits.open('mag_flux_tables/K/xray_mag_flux_table_best_extra_clean_no06.fits')[1].data
sdata = fits.open('mag_flux_tables/K/stars_mag_flux_table_extra_clean_no06.fits')[1].data
sigtb = Table.read('sigma_tables/quad_epoch_sigma_table_extra_clean_no06.fits')

### Remove edges ###
tbdata = vari_funcs.field_funcs.remove_edges(tbdata)
chandata = vari_funcs.field_funcs.remove_edges(chandata)
sdata = vari_funcs.field_funcs.remove_edges(sdata)

#### Limit to Chandra region ###
#tbdata = vari_funcs.field_funcs.chandra_only(tbdata)
#chandata = vari_funcs.field_funcs.chandra_only(chandata)
#sdata = vari_funcs.field_funcs.chandra_only(sdata)

## Create arrays of flux values ###
flux = vari_funcs.k_mag_flux.flux5_stacks(tbdata)
fluxchan = vari_funcs.k_mag_flux.flux5_stacks(chandata) 
sflux = vari_funcs.k_mag_flux.flux5_stacks(sdata)

### remove values that are negative ###
flux, tbdata = vari_funcs.flux_funcs.noneg(flux, tbdata)
fluxchan, chandata = vari_funcs.flux_funcs.noneg(fluxchan, chandata)
sflux, sdata = vari_funcs.flux_funcs.noneg(sflux, sdata)

### Get error arrays ###
flux, fluxerr, tbdata = vari_funcs.k_mag_flux.create_quad_error_array(sigtb, tbdata)
fluxchan, chanerr, chandata = vari_funcs.k_mag_flux.create_quad_error_array(sigtb, chandata)
sflux, serr, sdata = vari_funcs.k_mag_flux.create_quad_error_array(sigtb, sdata)

### reset X-ray column as messed up by stacking ###
tbdata['X-ray'][tbdata['X-ray']==70] = False 
tbdata['X-ray'][tbdata['X-ray']==84] = True

### Check chisq plot looks correct ###
fig,_ = vari_funcs.selection_plot_funcs.flux_variability_plot(flux, fluxchan, 'chisq', 
                                       fluxerr=fluxerr, chanerr=chanerr,
                                       starflux=sflux, starfluxerr=serr,
                                       #normalised=True, 
                                       stars=True, scale='log')
fig.canvas.mpl_connect('pick_event', vari_funcs.selection_plot_funcs.onpickflux)

### Calculate chi^2 values ###
chisq = vari_funcs.vary_stats.my_chisquare_err(flux, fluxerr)
chanchisq = vari_funcs.vary_stats.my_chisquare_err(fluxchan, chanerr)

#plt.hlines(22.458, 8e1, 1e7,zorder=4,label='Chi>22.5')
plt.hlines(30, 8e1, 1e7,'g', zorder=4,label='Chi=30')
#plt.hlines(40, 8e1, 1e7,'y', zorder=4,label='Chi>40')
#plt.hlines(50, 8e1, 1e7,'c', zorder=4,label='Chi>50')
plt.legend()
plt.tight_layout()

#### Plot on those that remain after deviant point check ###
#varys = fits.open('variable_tables/no06_variables_chi30_deviant.fits')[1].data
#flux = vari_funcs.k_mag_flux.flux5_stacks(varys)
#flux, varys = vari_funcs.noneg(flux, varys)
#flux, fluxerr, newvarys = vari_funcs.k_mag_flux.create_quad_error_array(sigtb, varys)
#meanflux = np.nanmean(flux, axis=1)
#chisq = vari_funcs.vary_stats.my_chisquare_err(flux, fluxerr)
#
#plt.plot(meanflux, chisq, 'ks', mfc='None', markersize=10)





















