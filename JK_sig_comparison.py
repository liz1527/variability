#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  1 13:30:42 2018

Code to compute maximum likelihood variability

@author: ppxee
"""
import time
start = time.time()
print(start)

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
#plt.close('all') #close any open plots
#from numpy.lib.recfunctions import append_fields

#%% Open the fits files and get data ###
tbdata = fits.open('variable_tables/no06_variables_chi30_2arcsec_DR11data_restframe.fits')[1].data
#chandata = fits.open('variable_tables/no06_variables_chi30_chandata_DR11data_restframe.fits')[1].data
sigtb = Table.read('sigma_tables/quad_epoch_sigma_table_extra_clean_no06_2arcsec.fits')
Jdata = fits.open('variable_tables/no06_variables_chi30_2arcsec_DR11data_J_best.fits')[1].data
Jsigtb = Table.read('sigma_tables/quad_epoch_sigma_table_extra_clean_2arcsec_J.fits')

#Jxraydata = Jdata[Jdata['X-ray']==True]

#### Limit to Chandra region for simplicity ###
#tbdata = vari_funcs.chandra_only(tbdata)
#chandata = vari_funcs.chandra_only(chandata)
#Jdata = vari_funcs.chandra_only(Jdata)
#Jxraydata = vari_funcs.chandra_only(Jxraydata)

# Extract magnitude table and error table
#flux = vari_funcs.flux4_stacks(tbdata)
#flux, tbdata = vari_funcs.noneg(flux, tbdata)
flux, fluxerr, tbdata = vari_funcs.create_quad_error_array(sigtb, tbdata, aper=4)
#chanflux = vari_funcs.flux4_stacks(chandata)
#chanflux, chandata = vari_funcs.noneg(chanflux, chandata)
#chanflux, chanerr, chandata = vari_funcs.create_quad_error_array(sigtb, chandata, aper=4)
#Jflux = vari_funcs.jflux4_stacks(Jdata)
#Jfluxerr = vari_funcs.jfluxerr4_stacks(Jdata)
#Jxrayflux = vari_funcs.flux4_stacks(Jxraydata)
#Jxrayfluxerr = vari_funcs.fluxerr4_stacks(Jxraydata)
Jflux, Jfluxerr, Jdata = vari_funcs.create_quad_error_array_J(Jsigtb, Jdata, aper=4)

### Normalise ###
fluxnorm, fluxerrnorm = vari_funcs.normalise_flux_and_errors(flux, fluxerr)
#chanfluxnorm, chanerrnorm = vari_funcs.normalise_flux_and_errors(chanflux, chanerr)
Jfluxnorm, Jfluxerrnorm = vari_funcs.normalise_flux_and_errors(Jflux, Jfluxerr)
#Jxrayfluxnorm, Jxrayfluxerrnorm = vari_funcs.normalise_flux_and_errors(Jxrayflux, Jxrayfluxerr)

#%% All points
posvar = np.linspace(0,2,5000)
meanflux = np.nanmean(fluxnorm, axis=1)
Jmeanflux = np.nanmean(Jfluxnorm, axis=1)

#start = time.time()
plt.figure()
for n in range(len(tbdata)):
    obnum = tbdata['NUMBER_1'][n]
    Jmask = np.isin(Jdata['NUMBER_1'], obnum)
    if ~np.any(Jmask):
        continue
    
    Kout = vari_funcs.maximum_likelihood(fluxnorm[n,:], fluxerrnorm[n,:], meanflux[n], posvar)
    
    Jmeanflux = np.nanmean(Jfluxnorm, axis=1)
    Jout = vari_funcs.maximum_likelihood(Jfluxnorm[Jmask,:].reshape(8), Jfluxerrnorm[Jmask,:].reshape(8), Jmeanflux[Jmask], posvar)

    plt.errorbar(Kout[0], Jout[0], xerr=Kout[1], yerr=Jout[1], fmt='.', color='tab:grey', zorder=0)
    if Jdata['X-ray'][Jmask] == True:
        plt.plot(Kout[0], Jout[0], 'rs', zorder=2)
    else:
        plt.plot(Kout[0], Jout[0], 'bo', zorder=1)
#
plt.xlabel('$\sigma_{K}$')
plt.ylabel('$\sigma_{J}$')
#plt.xscale('log')
#plt.yscale('log')
#numobs = np.shape(chanfluxnorm)[0]
#meanchan = np.nanmean(chanfluxnorm, axis=1)
#chanout = np.array([vari_funcs.maximum_likelihood(chanfluxnorm[n,:], chanerrnorm[n,:], meanchan[n], posvar) for n in range(numobs)])
#
#numobs = np.shape(Jfluxnorm)[0]
#Jmeanflux = np.nanmean(Jfluxnorm, axis=1)
#out = np.array([vari_funcs.maximum_likelihood(Jfluxnorm[n,:], Jfluxerrnorm[n,:], Jmeanflux[n], posvar) for n in range(numobs)])
#
#numobs = np.shape(Jxrayfluxnorm)[0]
#Jxraymeanflux = np.nanmean(Jxrayfluxnorm, axis=1)
#out = np.array([vari_funcs.maximum_likelihood(Jxrayfluxnorm[n,:], Jxrayfluxerrnorm[n,:], Jxraymeanflux[n], posvar) for n in range(numobs)])
#


x = np.linspace(0,2,10)
y = x
plt.plot(x,y,'k')
plt.tight_layout()

end = time.time()
print(end-start)















