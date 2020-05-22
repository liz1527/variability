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
monthout = varydata['Month_sig_K']
monthouterr = varydata['Month_sig_K_err']
semout = varydata['sig_K']
semouterr = varydata['sig_K_err']
xmonthout = xvarydata['Month_sig_K']
xmonthouterr = xvarydata['Month_sig_K_err']
noxmonthout = noxvarydata['Month_sig_K']
noxmonthouterr = noxvarydata['Month_sig_K_err']

### get flux data ###
monthflux = varydata['Month_Flux_K']
monthfluxerr = varydata['Month_Fluxerr_K']
semflux = varydata['Flux_K']
semfluxerr = varydata['Fluxerr_K']
xmonthflux = xvarydata['Month_Flux_K']
xmonthfluxerr = xvarydata['Month_Fluxerr_K']
noxmonthflux = noxvarydata['Month_Flux_K']
noxmonthfluxerr = noxvarydata['Month_Fluxerr_K']

### calculate mean flux ###
monthmeanflux = np.nanmean(monthflux, axis=1)
semmeanflux = np.nanmean(semflux, axis=1)
xmeanflux = np.nanmean(xmonthflux, axis=1)
noxmeanflux = np.nanmean(noxmonthflux, axis=1)

### calculate median flux ###
monthmedianflux = np.nanmedian(monthflux, axis=1)
semmedianflux = np.nanmedian(semflux, axis=1)
xmedianflux = np.nanmedian(xmonthflux, axis=1)
noxmedianflux = np.nanmedian(noxmonthflux, axis=1)

### plot sig vs flux ###
plt.figure()
plt.errorbar(noxmeanflux, noxmonthout, 
             yerr=noxmonthouterr, fmt='bo')
plt.errorbar(xmeanflux, xmonthout, 
             yerr=xmonthouterr, fmt='rs')
#plt.errorbar(noxmeanflux[noxmeanflux>1e4], noxmonthout[noxmeanflux>1e4], 
#             yerr=noxmonthouterr[noxmeanflux>1e4], fmt='bo')
#plt.errorbar(xmeanflux[xmeanflux>1e4], xmonthout[xmeanflux>1e4], 
#yerr=xmonthouterr[xmeanflux>1e4], fmt='rs')
plt.xscale('log')
#plt.xlim(xmin=1e4)
plt.yscale('log')
plt.xlabel('Mean K Band Flux')
plt.ylabel('$sig_{month}$')
plt.tight_layout()

### Plot mon and sem sig vs flux ###
plt.figure()
plt.plot(semmeanflux, semout, 'ko')
plt.plot(monthmeanflux, monthout, 'ms')
#plt.errorbar(semmeanflux, semout,'ko')
#plt.errorbar(monthmeanflux, monthout, 
#             yerr=monthouterr, fmt='ms')
#plt.errorbar(noxmeanflux[noxmeanflux>1e4], noxmonthout[noxmeanflux>1e4], 
#             yerr=noxmonthouterr[noxmeanflux>1e4], fmt='bo')
#plt.errorbar(xmeanflux[xmeanflux>1e4], xmonthout[xmeanflux>1e4], 
#yerr=xmonthouterr[xmeanflux>1e4], fmt='rs')
plt.xscale('log')
#plt.xlim(xmin=1e4)
plt.yscale('log')
plt.xlabel('Mean K Band Flux')
plt.ylabel('$sig_{month}$')
plt.tight_layout()


### Get chi data ###
#monthout = varydata['Month_sig_K']
#monthouterr = varydata['Month_sig_K_err']
xmonthout = xvarydata['Month_Chi_K']
noxmonthout = noxvarydata['Month_Chi_K']

end = time.time()
print(end-start)

plt.figure()
plt.errorbar(noxmeanflux, noxmonthout, 
             yerr=noxmonthouterr, fmt='bo')
plt.errorbar(xmeanflux, xmonthout, 
yerr=xmonthouterr, fmt='rs')
#plt.errorbar(noxmeanflux[noxmeanflux>1e4], noxmonthout[noxmeanflux>1e4], 
#             yerr=noxmonthouterr[noxmeanflux>1e4], fmt='bo')
#plt.errorbar(xmeanflux[xmeanflux>1e4], xmonthout[xmeanflux>1e4], 
#yerr=xmonthouterr[xmeanflux>1e4], fmt='rs')
plt.xscale('log')
#plt.xlim(xmin=1e4)
plt.yscale('log')
plt.xlabel('Mean K Band Flux')
plt.ylabel('$\chi^{2}_{month}$')
plt.tight_layout()

### plot mean fluxes ###
x = np.linspace(np.min(semmeanflux), np.max(semmeanflux),10)
y = x
plt.figure()
plt.plot(x,y,'k')
plt.plot(semmeanflux, monthmeanflux, 'o')
plt.xscale('log')
#plt.xlim(xmin=1e4)
plt.yscale('log')
plt.xlabel('Mean K Band Semester Flux')
plt.ylabel('Mean K Band Month Flux')
plt.tight_layout()


### plot median fluxes ###
x = np.linspace(np.min(semmedianflux), np.max(semmedianflux),10)
y = x
plt.figure()
plt.plot(x,y,'k')
plt.plot(semmedianflux, monthmedianflux, 'o')
plt.xscale('log')
#plt.xlim(xmin=1e4)
plt.yscale('log')
plt.xlabel('Median K Band Semester Flux')
plt.ylabel('Median K Band Month Flux')
plt.tight_layout()


#### plot weigth






