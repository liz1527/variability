#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  1 13:30:42 2018

Code to compute maximum likelihood variability

@author: ppxee
"""
import time
start = time.time()
#print(start)
print('Code was started '+time.ctime())

### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
from astropy.table import Table #for handling tables
import numpy as np #for handling arrays
#import math
#from astropy.stats import median_absolute_deviation
import vari_funcs #my module to help run code neatly
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
plt.close('all') #close any open plots
#from numpy.lib.recfunctions import append_fields

#%% Open the fits files and get data ###
#tbdata = fits.open('mag_flux_tables/mag_flux_table_best_extra_clean_no06.fits')[1].data
tbdata = fits.open('variable_tables/no06_variables_chi30_2arcsec_DR11data_restframe.fits')[1].data
#fullxray = fits.open('mag_flux_tables/chanDR11data_restframe_mag_flux_table_best_extra_clean_no06.fits')[1].data
#chandata = fits.open('variable_tables/no06_variables_chi30_chandata_DR11data_restframe.fits')[1].data
sigtb = Table.read('sigma_tables/quad_epoch_sigma_table_extra_clean_no06_2arcsec.fits')

### Limit to Chandra region for simplicity ###
#tbdata = vari_funcs.chandra_only(tbdata)
#fullxray = vari_funcs.chandra_only(fullxray)
#chandata = vari_funcs.chandra_only(chandata)

# Extract magnitude table and error table
flux = vari_funcs.flux4_stacks(tbdata)
flux, tbdata = vari_funcs.noneg(flux, tbdata)
flux, fluxerr, tbdata = vari_funcs.create_quad_error_array(sigtb, tbdata, aper=4)
#chanflux = vari_funcs.flux5_stacks(chandata)
#chanflux, chandata = vari_funcs.noneg(chanflux, chandata)
#chanflux, chanerr, chandata = vari_funcs.create_quad_error_array(sigtb, chandata, aper=5)
#fullflux = vari_funcs.flux5_stacks(fullxray)
#fullflux, fulldata = vari_funcs.noneg(fullflux, fullxray)
#fullflux, fullerr, fullxray = vari_funcs.create_quad_error_array(sigtb, fullxray, aper=5)

### Normalise ###
fluxnorm, fluxerrnorm = vari_funcs.normalise_flux_and_errors(flux, fluxerr)
#chanfluxnorm, chanerrnorm = vari_funcs.normalise_flux_and_errors(chanflux, chanerr)
#fullfluxnorm, fullerrnorm = vari_funcs.normalise_flux_and_errors(fullflux, fullerr)
#%% All points
posvar = np.linspace(0,2,5000)
#start = time.time()

numobs = np.shape(fluxnorm)[0]
meanflux = np.nanmean(fluxnorm, axis=1)
out = np.array([vari_funcs.maximum_likelihood(fluxnorm[n,:], fluxerrnorm[n,:], 
                                              meanflux[n], posvar, n=n, 
                                              printn=1000) for n in range(numobs)])

#numobs = np.shape(chanfluxnorm)[0]
#meanchan = np.nanmean(chanfluxnorm, axis=1)
#chanout = np.array([vari_funcs.maximum_likelihood(chanfluxnorm[n,:], chanerrnorm[n,:], meanchan[n], posvar) for n in range(numobs)])
#
#numobs = np.shape(fullfluxnorm)[0]
#meanfull = np.nanmean(fullfluxnorm, axis=1)
#fullout = np.array([vari_funcs.maximum_likelihood(fullfluxnorm[n,:], fullerrnorm[n,:], meanfull[n], posvar) for n in range(numobs)])

#%% Plots
mask=tbdata['X-ray'].astype(bool)
plt.figure(figsize=[7,7])
meanflux = np.nanmean(flux, axis=1)
plt.plot(meanflux, out[:,0], 'bo')
plt.errorbar(meanflux, out[:,0],yerr=out[:,1],fmt='+', color='tab:grey',alpha=0.5, zorder=0)
#plt.plot(meanflux[mask], out[mask,0], 'ro')#, markersize=10, mfc='None')
plt.xscale('log')
#plt.yscale('log')
plt.ylabel(r'$\sigma$')
plt.xlabel('Mean Flux')

plt.figure(figsize=[7,7])
meanflux = np.nanmean(flux, axis=1)
plt.plot(meanflux, out[:,0], 'bo')
plt.errorbar(meanflux, out[:,0],yerr=out[:,1],fmt='+', color='tab:grey',alpha=0.5, zorder=0)
#plt.plot(meanflux[mask], out[mask,0], 'ro')#, markersize=10, mfc='None')
plt.xscale('log')
plt.yscale('log')
plt.ylabel(r'$\sigma$')
plt.xlabel('Mean Flux')

end = time.time()
print(end-start)