#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 14:24:41 2020

Code to look at month variability properties

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
plt.close('all') #close any open plots
#from numpy.lib.recfunctions import append_fields

#%% Open the fits files and get data ###
### Import matched sample ###
varydata = Table.read('variable_tables/J_and_K_variables_month_varystats_DR11data.fits')
xvarydata = varydata[varydata['X-ray']==True]
noxvarydata = varydata[varydata['X-ray']==False]

#### Find stellarity for full tbdata ##
#phot = Table.read('UDS_catalogues/DR11-2arcsec-Jun-30-2019_best_photometry.fits')
#small = phot['MAG_APER'][:,0] # 0.7 arcsec
#big = phot['MAG_APER'][:,3] # 2 arcsec
#stell = big - small
#
#### set values where mag == 99 to nan ### 
#small[small==99] = np.nan
#big[big==99] = np.nan
#
#mask = np.isin(phot['ID'], varydata['ID'])
#
#varystell = stell[mask]
#xstell = varystell[varydata['X-ray']==True]
#noxstell = varystell[varydata['X-ray']==False]


### Get sig data ###
#monthout = varydata['Month_sig_K']
#monthouterr = varydata['Month_sig_K_err']
semout = varydata['sig_K']
semouterr = varydata['sig_K_err']
#xmonthout = xvarydata['Month_sig_K']
#xmonthouterr = xvarydata['Month_sig_K_err']
#noxmonthout = noxvarydata['Month_sig_K']
#noxmonthouterr = noxvarydata['Month_sig_K_err']

### get flux data ###
#monthflux = varydata['Month_Flux_K']
#monthfluxerr = varydata['Month_Fluxerr_K']
semflux = varydata['Flux_K']
semfluxerr = varydata['Fluxerr_K']
#xmonthflux = xvarydata['Month_Flux_K']
#xmonthfluxerr = xvarydata['Month_Fluxerr_K']
#noxmonthflux = noxvarydata['Month_Flux_K']
#noxmonthfluxerr = noxvarydata['Month_Fluxerr_K']

### calculate mean flux ###
#monthmeanflux = np.nanmean(monthflux, axis=1)
semmeanflux = np.nanmean(semflux, axis=1)
#xmeanflux = np.nanmean(xmonthflux, axis=1)
#noxmeanflux = np.nanmean(noxmonthflux, axis=1)

### calculate median flux ###
#monthmedianflux = np.nanmedian(monthflux, axis=1)
semmedianflux = np.nanmedian(semflux, axis=1)
#xmedianflux = np.nanmedian(xmonthflux, axis=1)
#noxmedianflux = np.nanmedian(noxmonthflux, axis=1)

posvar = np.linspace(0,40000,1000)
#posvar = np.linspace(0,4,100)

### Set up arrays for K selected ###
allsemout = np.array([])
allsemouterr = np.array([])

semfluxn, semfluxerrn = vari_funcs.flux_funcs.normalise_flux_and_errors(semflux, semfluxerr)
semmeanfluxn = np.nanmean(semfluxn, axis=1)

for n in range(len(varydata)): #loop over the selection band
    
    ### Get maximum likelihoods in K for that object ###
    KKout = vari_funcs.vary_stats.maximum_likelihood(semflux[n,:], 
                                                     semfluxerr[n,:], 
                                                     semmeanflux[n], posvar)
#    KKout = vari_funcs.vary_stats.maximum_likelihood(semfluxn[n,:], 
#                                                     semfluxerrn[n,:], 
#                                                     semmeanfluxn[n], posvar)
    

    ### save output into x-ray and band specfic arrays ###
    allsemout = np.append(allsemout, KKout[0])
    allsemouterr = np.append(allsemouterr, KKout[1])
    
### calculate sig_excess ###
#semfluxn, semfluxerrn = vari_funcs.flux_funcs.normalise_flux_and_errors(semflux, semfluxerr)
#allsemout = np.sqrt(vari_funcs.vary_stats.normsigmasq(semfluxn, semfluxerrn))

### plot sig_maxl vs sig_excess ###
plt.figure()
plt.errorbar(semout, allsemout, xerr=semouterr, fmt='o', color='tab:grey', alpha=0.5)
plt.plot(semout, allsemout, 'o')
plt.xscale('log')
#plt.xlim(xmin=1e4)
#plt.yscale('log')
#plt.xlabel('$sig_{max likely}$')
#plt.ylabel('$sig_{excess}$')
plt.xlabel('$sig_{norm}$')
plt.ylabel('$sig_{notnorm}$')
plt.tight_layout()


plt.figure()
plt.plot(semmeanflux, allsemout, 'o')
plt.xscale('log')
#plt.xlim(xmin=1e4)
plt.yscale('log')
plt.xlabel('Mean K band flux')
plt.ylabel('$sig_{not norm}$')
plt.tight_layout()

plt.figure()
plt.plot(semmeanflux, semout, 'o')
plt.xscale('log')
#plt.xlim(xmin=1e4)
plt.yscale('log')
plt.xlabel('Mean K band flux')
plt.ylabel('$sig_{norm}$')
plt.tight_layout()

### calculate chi ###
allchi = np.sqrt(vari_funcs.vary_stats.my_chisquare_err(semfluxn, semfluxerrn))
chi = varydata['Chi_K']

plt.figure()
plt.plot(chi, allchi, 'o')
plt.xscale('log')
#plt.xlim(xmin=1e4)
plt.yscale('log')
plt.xlabel('$\chi_{norm}$')
plt.ylabel('$\chi_{not norm}$')
plt.tight_layout()


plt.figure()
plt.plot(semmeanflux, allchi, 'o')
plt.xscale('log')
#plt.xlim(xmin=1e4)
plt.yscale('log')
plt.xlabel('Mean K band flux')
plt.ylabel('$\chi_{not norm}$')
plt.tight_layout()

plt.figure()
plt.plot(semmeanflux, chi, 'o')
plt.xscale('log')
#plt.xlim(xmin=1e4)
plt.yscale('log')
plt.xlabel('Mean K band flux')
plt.ylabel('$\chi_{norm}$')
plt.tight_layout()


#### plot mean fluxes ###
#x = np.linspace(np.min(semmeanflux), np.max(semmeanflux),10)
#y = x
#plt.figure()
#plt.plot(x,y,'k')
#plt.plot(semmeanflux, monthmeanflux, 'o')
#plt.xscale('log')
##plt.xlim(xmin=1e4)
#plt.yscale('log')
#plt.xlabel('Mean K Band Semester Flux')
#plt.ylabel('Mean K Band Month Flux')
#plt.tight_layout()
#
#
#### plot median fluxes ###
#x = np.linspace(np.min(semmedianflux), np.max(semmedianflux),10)
#y = x
#plt.figure()
#plt.plot(x,y,'k')
#plt.plot(semmedianflux, monthmedianflux, 'o')
#plt.xscale('log')
##plt.xlim(xmin=1e4)
#plt.yscale('log')
#plt.xlabel('Mean K Band Semester Flux')
#plt.ylabel('Mean K Band Month Flux')
#plt.tight_layout()
#
#
#### plot weighted mean fluxes ###
#monthwmeanflux = np.average(monthflux, weights=monthfluxerr, axis=1)
#semwmeanflux = np.average(semflux, weights=semfluxerr, axis=1)
#x = np.linspace(np.min(semwmeanflux), np.max(semwmeanflux),10)
#y = x
#plt.figure()
#plt.plot(x,y,'k')
##plt.plot(semwmeanflux, monthwmeanflux, 'o')
##plt.plot(semwmeanflux, semmeanflux, 'o')
#plt.xscale('log')
##plt.xlim(xmin=1e4)
#plt.yscale('log')
#plt.xlabel('Mean K Band Semester Flux')
#plt.ylabel('Mean K Band Month Flux')
#plt.tight_layout()







