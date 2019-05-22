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

posvar = np.linspace(0,2,5000)

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
        if m != len(binedge)-1:
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
        
#%% Split by hardness ratio ###
xmmhr = xmmdata['HR2']
medhr = np.median(xmmhr)
xmmhard = xmmdata[xmmhr >= medhr]
xmmsoft = xmmdata[xmmhr < medhr]
#%%
#chandata, xrayL, chanfluxnorm, chanerrnorm= get_luminosity_and_flux(chandata)
#fullxray, fullxrayL, fullfluxnorm, fullerrnorm = get_luminosity_and_flux(fullxray)
xmmhard, hardxrayL, hardfluxnorm, harderrnorm = get_luminosity_and_flux(xmmhard, xmm=True)
xmmsoft, softxrayL, softfluxnorm, softerrnorm = get_luminosity_and_flux(xmmsoft, xmm=True)

#%% Split into bins according to X-ray luminosity ###
#allxrayL = np.append(hardxrayL, softxrayL)
sortedhardxrayL = np.sort(hardxrayL)
sortedsoftxrayL = np.sort(softxrayL)
#sortedxrayL = np.sort(fullxrayL)

binsize = 10
hardbinedge = np.empty(int(len(softxrayL)/binsize))
softbinedge = np.empty(int(len(softxrayL)/binsize))
#binedge = np.empty(int(len(fullxrayL)/binsize))
size = len(hardbinedge)
for n in range(size):
    ### Define bin edges ###
    enXrayL = sortedhardxrayL[n*binsize:(n*binsize)+binsize]
    enmin = np.nanmin(enXrayL)
    hardbinedge[n] = enmin
size = len(softbinedge)
for n in range(size):
    ### Define bin edges ###
    enXrayL = sortedsoftxrayL[n*binsize:(n*binsize)+binsize]
    enmin = np.nanmin(enXrayL)
    softbinedge[n] = enmin
#size = len(binedge)
#for n in range(size):
#    ### Define bin edges ###
#    enXrayL = sortedxrayL[n*binsize:(n*binsize)+binsize]
#    enmin = np.nanmin(enXrayL)
#    binedge[n] = enmin

#enflux, enfluxerr, enXrayL = make_ensemble(chandata, xrayL, binedge)
#fenflux, fenfluxerr, fenXrayL = make_ensemble(fullxray, fullxrayL, binedge)
hardenflux, hardenfluxerr, hardenXrayL = make_ensemble(xmmhard, hardxrayL, hardbinedge)
softenflux, softenfluxerr, softenXrayL = make_ensemble(xmmsoft, softxrayL, softbinedge)

#
##### Combine chan and xmm ##
##xenflux = {}
##xenfluxerr = {}
##xenXrayL = {}
##for m, enmin in enumerate(binedge):
##    xenflux[m] = np.vstack([enflux[m],xmmenflux[m]])
##    xenfluxerr[m] = np.vstack([enfluxerr[m],xmmenfluxerr[m]])
##    xenXrayL[m] = np.append(enXrayL[m],xmmenXrayL[m])
#
### Get ensemble maximum likelihood ###
hardsig, hardsigerr, hardmeanxrayL = get_ensemble_sig(hardenflux, hardenfluxerr, hardenXrayL, posvar)
softsig, softsigerr, softmeanxrayL = get_ensemble_sig(softenflux, softenfluxerr, softenXrayL, posvar)
#fsig, fsigerr, fmeanxrayL = get_ensemble_sig(fenflux, fenfluxerr, fenXrayL, posvar)

#%% Plot results ###
hardbinupper = np.append(hardbinedge[1:],np.max(sortedhardxrayL))
softbinupper = np.append(softbinedge[1:],np.max(sortedsoftxrayL))
#binupper = np.append(binedge[1:],np.max(sortedxrayL))
#fxlow = fmeanxrayL-binedge
#fxhigh = binupper - fmeanxrayL
hardxlow = hardmeanxrayL-hardbinedge
hardxhigh = hardbinupper - hardmeanxrayL
softxlow = softmeanxrayL-softbinedge
softxhigh = softbinupper - softmeanxrayL

plt.figure(figsize=[8,7])
#plt.errorbar(fmeanxrayL, fsig, xerr=[fxlow, fxhigh], yerr=fsigerr, fmt='o',
#             color='k', label='Ensemble Non Variable X-ray Sources')
plt.errorbar(softmeanxrayL, softsig, xerr=[softxlow, softxhigh], yerr=softsigerr, fmt='o',
             color='b', label='HR2 >= -0.625')
plt.errorbar(hardmeanxrayL, hardsig, xerr=[hardxlow, hardxhigh], yerr=hardsigerr, fmt='o',
             color='g', label='HR2 < -0.625')
plt.xscale('log')
plt.xlabel('X-ray Luminosity (ergs/s)')
plt.ylabel(r'$\sigma$')
plt.ylim(ymin=0,ymax=0.15)
#plt.xlim(xmin=1e41)
plt.legend()





















end = time.time()
print(end-start)

