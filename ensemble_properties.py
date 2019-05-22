#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 29 12:13:52 2018

Code to look at ensemble variability of objects

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
tbdata = fits.open('variable_tables/no06_variables_chi30_DR11data_restframe.fits')[1].data
fullxray = fits.open('mag_flux_tables/novarys_chanDR11data_restframe_mag_flux_table_best_extra_clean_no06.fits')[1].data
sigtb = Table.read('sigma_tables/quad_epoch_sigma_table_extra_clean_no06.fits')

#### Limit to Chandra region for simplicity ###
#tbdata = vari_funcs.chandra_only(tbdata)
#fullxray = vari_funcs.chandra_only(fullxray)
#chandata = vari_funcs.chandra_only(chandata)

# Extract magnitude table and error table
flux = vari_funcs.flux5_stacks(tbdata)
flux, tbdata = vari_funcs.noneg(flux, tbdata)
flux, fluxerr, tbdata = vari_funcs.create_quad_error_array(sigtb, tbdata)
chanflux = vari_funcs.flux5_stacks(chandata)
chanflux, chandata = vari_funcs.noneg(chanflux, chandata)
chanflux, chanerr, chandata = vari_funcs.create_quad_error_array(sigtb, chandata)
fullflux = vari_funcs.flux5_stacks(fullxray)
fullflux, fulldata = vari_funcs.noneg(fullflux, fullxray)
fullflux, fullerr, fullxray = vari_funcs.create_quad_error_array(sigtb, fullxray)
xmmflux = vari_funcs.flux5_stacks(xmmdata)
xmmflux, xmmdata = vari_funcs.noneg(xmmflux, xmmdata)
xmmflux, xmmerr, xmmdata = vari_funcs.create_quad_error_array(sigtb, xmmdata)

### Normalise ###
fluxnorm, fluxerrnorm = vari_funcs.normalise_flux_and_errors(flux, fluxerr)
chanfluxnorm, chanerrnorm = vari_funcs.normalise_flux_and_errors(chanflux, chanerr)
fullfluxnorm, fullerrnorm = vari_funcs.normalise_flux_and_errors(fullflux, fullerr)
xmmfluxnorm, xmmerrnorm = vari_funcs.normalise_flux_and_errors(xmmflux, xmmerr)


#%% Find luminosity distance both for all and just for chandra sources ###
z = tbdata['z_spec']#[mask]
z[z==-1] = tbdata['z_p'][z==-1]
DL = cosmo.luminosity_distance(z)
DL = DL.to(u.cm)
chanz = chandata['z_spec']#[chanmask]
chanz[chanz==-1] = chandata['z_p'][chanz==-1]
chanDL = cosmo.luminosity_distance(chanz)
chanDL = chanDL.to(u.cm)
fullz = fullxray['z_spec']#[fullmask]
fullz[fullz==-1] = fullxray['z_p'][fullz==-1]
fullDL = cosmo.luminosity_distance(fullz)
fullDL = fullDL.to(u.cm)
xmmz = xmmdata['z_spec']#[chanmask]
xmmz[xmmz==-1] = xmmdata['z_p'][xmmz==-1]
xmmDL = cosmo.luminosity_distance(xmmz)
xmmDL = xmmDL.to(u.cm)

#%% Calculate the luminosity ###
xrayF = chandata['Soft_flux']#[chanmask]
xrayL = xrayF*4*np.pi*(chanDL.value**2)
fullxrayF = fullxray['Soft_flux']#[fullmask]
fullxrayL = fullxrayF*4*np.pi*(fullDL.value**2)
xmmxrayF = xmmdata['CR(S)'] * 0.171 * (10**(-14)) #conversion into ergs/s/cm2
xmmxrayL = xmmxrayF*4*np.pi*(xmmDL.value**2)

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

### Create dicts to save data into ###
enflux = {}
enfluxerr = {}
enXrayL = {}
sig = np.empty(size)
sigerr = np.empty(size)
meanxrayL = np.empty(size)
fenflux = {}
fenfluxerr = {}
fenXrayL = {}
fsig = np.empty(size)
fsigerr = np.empty(size)
fmeanxrayL = np.empty(size)
xmmenflux = {}
xmmenfluxerr = {}
xmmenXrayL = {}
xmmsig = np.empty(size)
xmmsigerr = np.empty(size)
xmmmeanxrayL = np.empty(size)
xenflux = {}
xenfluxerr = {}
xenXrayL = {}
xsig = np.empty(size)
xsigerr = np.empty(size)
xmeanxrayL = np.empty(size)
for m, enmin in enumerate(binedge):
    ### Isolate data needed ###
    mask1 = xrayL >= enmin
    fmask1 = fullxrayL >enmin
    xmmmask1 = xmmxrayL >= enmin
    if m != size-1:
        mask2 = xrayL < binedge[m+1]
        fmask2 = fullxrayL < binedge[m+1]
        xmmmask2 = xmmxrayL < binedge[m+1]
    else:
        mask2 = np.ones(len(mask1))
        fmask2 = np.ones(len(fmask1))
        xmmmask2 = np.ones(len(xmmmask1))
    
    enmask = mask1*mask2.astype(bool)
    fenmask = fmask1*fmask2.astype(bool)
    xmmenmask = xmmmask1*xmmmask2.astype(bool)

    enflux[m] = chanfluxnorm[enmask]
    enfluxerr[m] = chanerrnorm[enmask]
    enXrayL[m] = xrayL[enmask]
    
    fenflux[m] = fullfluxnorm[fenmask]
    fenfluxerr[m] = fullerrnorm[fenmask]
    fenXrayL[m] = fullxrayL[fenmask]
    
    xmmenflux[m] = xmmfluxnorm[xmmenmask]
    xmmenfluxerr[m] = xmmerrnorm[xmmenmask]
    xmmenXrayL[m] = xmmxrayL[xmmenmask]
    
    ### Combine chan and xmm ###
    xenflux[m] = np.vstack([enflux[m],xmmenflux[m]])
    xenfluxerr[m] = np.vstack([enfluxerr[m],xmmenfluxerr[m]])
    xenXrayL[m] = np.append(enXrayL[m],xmmenXrayL[m])
    
    ### Combine into one flux curve per bin ###
#    enfluxcurve = np.ravel(enflux[m])
#    enfluxcurveerr = np.ravel(enfluxerr[m])
    fenfluxcurve = np.ravel(fenflux[m])
    fenfluxcurveerr = np.ravel(fenfluxerr[m])
#    xmmenfluxcurve = np.ravel(xmmenflux[m])
#    xmmenfluxcurveerr = np.ravel(xmmenfluxerr[m])
    xenfluxcurve = np.ravel(xenflux[m])
    xenfluxcurveerr = np.ravel(xenfluxerr[m])
    
    
    ### Find max likelihood sig of curve ###
#    [sig[m],sigerr[m]] = vari_funcs.maximum_likelihood(enfluxcurve, enfluxcurveerr, 1, posvar)
    [fsig[m],fsigerr[m]] = vari_funcs.maximum_likelihood(fenfluxcurve, fenfluxcurveerr, 1, posvar)
#    [xmmsig[m],xmmsigerr[m]] = vari_funcs.maximum_likelihood(xmmenfluxcurve, xmmenfluxcurveerr, 1, posvar)
    [xsig[m],xsigerr[m]] = vari_funcs.maximum_likelihood(xenfluxcurve, xenfluxcurveerr, 1, posvar)

    ### find mean xrayL ###
#    meanxrayL[m] = np.nanmean(enXrayL[m])
    fmeanxrayL[m] = np.nanmean(fenXrayL[m])
#    xmmmeanxrayL[m] = np.nanmean(xmmenXrayL[m])
    xmeanxrayL[m] = np.nanmean(xenXrayL[m])

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
#xlow = meanxrayL-binedge
#xhigh = binupper - meanxrayL
fxlow = fmeanxrayL-binedge
fxhigh = binupper - fmeanxrayL
#xmmxlow = xmmmeanxrayL-binedge
#xmmxhigh = binupper - xmmmeanxrayL
xxlow = xmeanxrayL-binedge
xxhigh = binupper - xmeanxrayL
plt.figure(figsize=[8,7])
plt.errorbar(xrayL, chanout[:,0],yerr=chanout[:,1],fmt='o', alpha=0.25,
             color='r',label='Variable X-ray Source',zorder=2)
plt.errorbar(xmmxrayL, xmmout[:,0],yerr=xmmout[:,1],fmt='o', alpha=0.25,
             color='r',zorder=2)
plt.errorbar(fullxrayL, fullout[:,0],yerr=fullout[:,1],fmt='o', alpha=0.25,
             color='tab:grey', label='Non variable X-ray source',zorder=1)
#plt.errorbar(meanxrayL, sig, xerr=[xlow, xhigh], yerr=sigerr, fmt='o',
#             color='r', label='Ensemble Variable X-ray Sources')
plt.errorbar(fmeanxrayL, fsig, xerr=[fxlow, fxhigh], yerr=fsigerr, fmt='o',
             color='k', label='Ensemble Non Variable X-ray Sources')
#plt.errorbar(xmmmeanxrayL, xmmsig, xerr=[xmmxlow, xmmxhigh], yerr=xmmsigerr, fmt='o',
#             color='r', label='Ensemble Variable X-ray Sources')
plt.errorbar(xmeanxrayL, xsig, xerr=[xxlow, xxhigh], yerr=xsigerr, fmt='o',
             color='r', label='Ensemble Variable X-ray Sources')
plt.xscale('log')
plt.xlabel('X-ray Luminosity (ergs/s)')
plt.ylabel(r'$\sigma$')
plt.ylim(ymin=-0.1,ymax=0.5)
plt.xlim(xmin=5e41)
plt.legend()





















end = time.time()
print(end-start)

