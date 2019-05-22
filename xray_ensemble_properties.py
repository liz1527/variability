#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 29 12:13:52 2018

Code to look at ensemble variability of objects compared to X-ray luminosity,
with a split based on hardness ratio

@author: ppxee
"""

import time
start = time.time()
#print(start)

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

### Define cosmology ###
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

#%% Open the fits files and get data ###
chandata = fits.open('variable_tables/no06_variables_chi30_chandata_DR11data_restframe.fits')[1].data
xmmdata = fits.open('variable_tables/no06_variables_chi30_xmmdata_DR11data_restframe.fits')[1].data
fullxray = fits.open('mag_flux_tables/novarys_chanDR11data_restframe_mag_flux_table_best_extra_clean_no06.fits')[1].data
sigtb = Table.read('sigma_tables/quad_epoch_sigma_table_extra_clean_no06.fits')


#%% create function to do things 
def get_luminosity_and_flux(tbdata, xmm=False):
#    
    ### Extract magnitude table and error table ###
    flux = vari_funcs.flux5_stacks(tbdata)
    flux, tbdata = vari_funcs.noneg(flux, tbdata)
    flux, fluxerr, tbdata = vari_funcs.create_quad_error_array(sigtb, tbdata)
    
    ### Normalise ###
    fluxnorm, fluxerrnorm = vari_funcs.normalise_flux_and_errors(flux, fluxerr)
    
    ### Find luminosity distance ###
    z = tbdata['z_spec']#[mask]
    z[z==-1] = tbdata['z_p'][z==-1]
    DL = cosmo.luminosity_distance(z)
    DL = DL.to(u.cm)
    
    ### Calculate luminosity ###
    if xmm == True:
        xrayF = tbdata['CR(S)'] * 0.171 * (10**(-14)) #conversion into ergs/s/cm2
    else:
        xrayF = tbdata['Soft_flux'] #no conversion required in chandra
    xrayL = xrayF*4*np.pi*(DL.value**2)
    
    return tbdata, xrayL, fluxnorm, fluxerrnorm

def make_ensemble(tbdata, xrayL, binedge):
    
    ### Extract magnitude table and error table ###
    flux = vari_funcs.flux5_stacks(tbdata)
    flux, tbdata = vari_funcs.noneg(flux, tbdata)
    flux, fluxerr, tbdata = vari_funcs.create_quad_error_array(sigtb, tbdata)
    
    ### Normalise ###
    fluxnorm, fluxerrnorm = vari_funcs.normalise_flux_and_errors(flux, fluxerr)
    
    ### Create dicts to save data into ###
    enflux = {}
    enfluxerr = {}
    enXrayL = {}
    
    ### loop over bins ###
    for m, enmin in enumerate(binedge):
        ### Isolate data needed ###
        mask1 = xrayL >= enmin
        if m != size-1:
            mask2 = xrayL < binedge[m+1]
        else:
            mask2 = np.ones(len(mask1))
        enmask = mask1*mask2.astype(bool)
        enflux[m] = fluxnorm[enmask]
        enfluxerr[m] = fluxerrnorm[enmask]
        enXrayL[m] = xrayL[enmask]
        
    return enflux, enfluxerr, enXrayL

def get_ensemble_sig(enflux, enfluxerr, enXrayL, posvar):
    ### Create dicts to save data into ###
    size = len(enflux)
    sig = np.empty(size)
    sigerr = np.empty(size)
    meanxrayL = np.empty(size)
    for m in enflux:
        ### Combine into one flux curve per bin ###
        enfluxcurve = np.ravel(enflux[m])
        enfluxcurveerr = np.ravel(enfluxerr[m])
        
        ### Find max likelihood sig of curve ###
        [sig[m],sigerr[m]] = vari_funcs.maximum_likelihood(enfluxcurve, enfluxcurveerr, 1, posvar)

        ### find mean xrayL ###
        meanxrayL[m] = np.nanmean(enXrayL[m])
    
    return sig, sigerr, meanxrayL
        
#%%
chandata, xrayL, chanfluxnorm, chanerrnorm= get_luminosity_and_flux(chandata)
fullxray, fullxrayL, fullfluxnorm, fullerrnorm = get_luminosity_and_flux(fullxray)
xmmdata, xmmxrayL, xmmfluxnorm, xmmerrnorm = get_luminosity_and_flux(xmmdata, xmm=True)

#%% Split into bins according to X-ray luminosity ###
allxrayL = np.append(xrayL, xmmxrayL)
sortedxrayL = np.sort(allxrayL)
posvar = np.linspace(0,2,5000)

binsize = 14
binedge = np.empty(int(len(allxrayL)/binsize))
size = len(binedge)
for n in range(size):
    ### Define bin edges ###
    enXrayL = sortedxrayL[n*binsize:(n*binsize)+binsize]
    enmin = np.nanmin(enXrayL)
    binedge[n] = enmin

enflux, enfluxerr, enXrayL = make_ensemble(chandata, xrayL, binedge)
fenflux, fenfluxerr, fenXrayL = make_ensemble(fullxray, fullxrayL, binedge)
xmmenflux, xmmenfluxerr, xmmenXrayL = make_ensemble(xmmdata, xmmxrayL, binedge)

### Combine chan and xmm ##
xenflux = {}
xenfluxerr = {}
xenXrayL = {}
for m, enmin in enumerate(binedge):
    xenflux[m] = np.vstack([enflux[m],xmmenflux[m]])
    xenfluxerr[m] = np.vstack([enfluxerr[m],xmmenfluxerr[m]])
    xenXrayL[m] = np.append(enXrayL[m],xmmenXrayL[m])

### Get ensemble maximum likelihood ###
xsig, xsigerr, xmeanxrayL = get_ensemble_sig(xenflux, xenfluxerr, xenXrayL, posvar)
fsig, fsigerr, fmeanxrayL = get_ensemble_sig(fenflux, fenfluxerr, fenXrayL, posvar)

### Get non ensemble results ###
numobs = np.shape(chanfluxnorm)[0]
meanchan = np.nanmean(chanfluxnorm, axis=1)
chanout = np.array([vari_funcs.maximum_likelihood(chanfluxnorm[n,:], chanerrnorm[n,:], meanchan[n], posvar) for n in range(numobs)])

numobs = np.shape(fullfluxnorm)[0]
meanfull = np.nanmean(fullfluxnorm, axis=1)
fullout = np.array([vari_funcs.maximum_likelihood(fullfluxnorm[n,:], fullerrnorm[n,:], meanfull[n], posvar) for n in range(numobs)])

numobs = np.shape(xmmfluxnorm)[0]
meanxmm = np.nanmean(xmmfluxnorm, axis=1)
xmmout = np.array([vari_funcs.maximum_likelihood(xmmfluxnorm[n,:], xmmerrnorm[n,:], meanxmm[n], posvar) for n in range(numobs)])

#%% Plot results ###
binupper = np.append(binedge[1:],np.max(sortedxrayL))
fxlow = fmeanxrayL-binedge
fxhigh = binupper - fmeanxrayL
xxlow = xmeanxrayL-binedge
xxhigh = binupper - xmeanxrayL
plt.figure(figsize=[12,7])
plt.subplot(121)
plt.errorbar(xrayL, chanout[:,0],yerr=chanout[:,1],fmt='o', alpha=0.25,
             color='r',zorder=2)#,label='Variable X-ray Source')
plt.errorbar(xmmxrayL, xmmout[:,0],yerr=xmmout[:,1],fmt='o', alpha=0.25,
             color='r',zorder=2)
plt.errorbar(fullxrayL, fullout[:,0],yerr=fullout[:,1],fmt='o', alpha=0.25,
             color='tab:grey',zorder=1)#, label='Non variable X-ray source')
plt.scatter(xrayL, chanout[:,0],
             c='r',label='Variable X-ray Source',zorder=2)
plt.scatter(xmmxrayL, xmmout[:,0], c='r',zorder=2)
plt.scatter(fullxrayL, fullout[:,0],
             c='tab:grey', label='Non variable X-ray source',zorder=1)
plt.xscale('log')
plt.xlabel('X-ray Luminosity (ergs/s)')
plt.ylabel(r'$\sigma$')
plt.ylim(ymin=-0.05,ymax=0.5)
plt.xlim(xmin=5e41)
plt.legend()

plt.subplot(122)
plt.errorbar(fmeanxrayL, fsig, xerr=[fxlow, fxhigh], yerr=fsigerr, fmt='o',
             color='k', label='Ensemble Non Variable X-ray Sources')
plt.errorbar(xmeanxrayL, xsig, xerr=[xxlow, xxhigh], yerr=xsigerr, fmt='o',
             color='r', label='Ensemble Variable X-ray Sources')
plt.xscale('log')
plt.xlabel('X-ray Luminosity (ergs/s)')
plt.ylabel(r'$\sigma$')
plt.ylim(ymin=-0.05,ymax=0.5)
plt.xlim(xmin=5e41)
plt.legend()

plt.tight_layout()




















end = time.time()
print(end-start)

