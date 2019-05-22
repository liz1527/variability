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
tbdata = fits.open('variable_tables/no06_variables_chi30_DR11data_restframe.fits')[1].data
fullxray = fits.open('mag_flux_tables/chanDR11data_restframe_mag_flux_table_best_extra_clean_no06.fits')[1].data
sigtb = Table.read('quad_epoch_sigma_table_extra_clean_no06_067arcsec.fits')

### Limit to Chandra region for simplicity ###
tbdata = vari_funcs.chandra_only(tbdata)
fullxray = vari_funcs.chandra_only(fullxray)
chandata = vari_funcs.chandra_only(chandata)

# Extract magnitude table and error table
flux = vari_funcs.flux1_stacks(tbdata)
flux, tbdata = vari_funcs.noneg(flux, tbdata)
flux, fluxerr, tbdata = vari_funcs.create_quad_error_array(sigtb, tbdata,aper=1)
chanflux = vari_funcs.flux1_stacks(chandata)
chanflux, chandata = vari_funcs.noneg(chanflux, chandata)
chanflux, chanerr, chandata = vari_funcs.create_quad_error_array(sigtb, chandata,aper=1)
fullflux = vari_funcs.flux1_stacks(fullxray)
fullflux, fulldata = vari_funcs.noneg(fullflux, fullxray)
fullflux, fullerr, fullxray = vari_funcs.create_quad_error_array(sigtb, fullxray,aper=1)

### Normalise ###
fluxnorm, fluxerrnorm = vari_funcs.normalise_flux_and_errors(flux, fluxerr)
chanfluxnorm, chanerrnorm = vari_funcs.normalise_flux_and_errors(chanflux, chanerr)
fullfluxnorm, fullerrnorm = vari_funcs.normalise_flux_and_errors(fullflux, fullerr)


#%% Find z ###
z = tbdata['z_spec']#[mask]
z[z==-1] = tbdata['z_p'][z==-1]
chanz = chandata['z_spec']#[chanmask]
chanz[chanz==-1] = chandata['z_p'][chanz==-1]
fullz = fullxray['z_spec']#[fullmask]
fullz[fullz==-1] = fullxray['z_p'][fullz==-1]

#%% Split into bins according to z ###
sortedchanz = np.sort(chanz)
posvar = np.linspace(0,2,5000)

binsize = 10
binedge = np.empty(int(len(chanz)/binsize))
size = len(binedge)
for n in range(size):
    ### Define bin edges ###
    enz = sortedchanz[n*binsize:(n*binsize)+binsize]
    enmin = np.nanmin(enz)
    binedge[n] = enmin

### Create dicts to save data into ###
enflux = {}
enfluxerr = {}
enz = {}
sig = np.empty(size)
sigerr = np.empty(size)
meanz = np.empty(size)
cenflux = {}
cenfluxerr = {}
cenz = {}
csig = np.empty(size)
csigerr = np.empty(size)
cmeanz = np.empty(size)
fenflux = {}
fenfluxerr = {}
fenz = {}
fsig = np.empty(size)
fsigerr = np.empty(size)
fmeanz = np.empty(size)
for m, enmin in enumerate(binedge):
    ### Isolate data needed ###
    mask1 = z >= enmin
    cmask1 = chanz >= enmin
    fmask1 = fullz >= enmin
    if m != size-1:
        mask2 = z < binedge[m+1]
        cmask2 = chanz < binedge[m+1]
        fmask2 = fullz < binedge[m+1]
    else:
        mask2 = z < 4.5#np.ones(len(mask1))
        cmask2 = chanz < 4.5#np.ones(len(cmask1))
        fmask2 = fullz < 4.5#np.ones(len(fmask1))
    
    enmask = mask1*mask2.astype(bool)
    cenmask = cmask1*cmask2.astype(bool)
    fenmask = fmask1*fmask2.astype(bool)

    enflux[m] = fluxnorm[enmask]
    enfluxerr[m] = fluxerrnorm[enmask]
    enz[m] = z[enmask]
    
    cenflux[m] = chanfluxnorm[cenmask]
    cenfluxerr[m] = chanerrnorm[cenmask]
    cenz[m] = chanz[cenmask]
    
    fenflux[m] = fullfluxnorm[fenmask]
    fenfluxerr[m] = fullerrnorm[fenmask]
    fenz[m] = fullz[fenmask]
    
    ### Combine into one flux curve per bin ###
    enfluxcurve = np.ravel(enflux[m])
    enfluxcurveerr = np.ravel(enfluxerr[m])
    cenfluxcurve = np.ravel(cenflux[m])
    cenfluxcurveerr = np.ravel(cenfluxerr[m])
    fenfluxcurve = np.ravel(fenflux[m])
    fenfluxcurveerr = np.ravel(fenfluxerr[m])
    
    
    ### Find max likelihood sig of curve ###
    [sig[m],sigerr[m]] = vari_funcs.maximum_likelihood(enfluxcurve, enfluxcurveerr, 1, posvar)
    [csig[m],csigerr[m]] = vari_funcs.maximum_likelihood(cenfluxcurve, cenfluxcurveerr, 1, posvar)
    [fsig[m],fsigerr[m]] = vari_funcs.maximum_likelihood(fenfluxcurve, fenfluxcurveerr, 1, posvar)
    
    ### find mean z ###
    meanz[m] = np.nanmean(enz[m])
    cmeanz[m] = np.nanmean(cenz[m])
    fmeanz[m] = np.nanmean(fenz[m])

### Get non ensemble results ###
numobs = np.shape(fluxnorm)[0]
meanflux = np.nanmean(fluxnorm, axis=1)
out = np.array([vari_funcs.maximum_likelihood(fluxnorm[n,:], fluxerrnorm[n,:], meanflux[n], posvar) for n in range(numobs)])

numobs = np.shape(chanfluxnorm)[0]
meanchan = np.nanmean(chanfluxnorm, axis=1)
chanout = np.array([vari_funcs.maximum_likelihood(chanfluxnorm[n,:], chanerrnorm[n,:], meanchan[n], posvar) for n in range(numobs)])

numobs = np.shape(fullfluxnorm)[0]
meanfull = np.nanmean(fullfluxnorm, axis=1)
fullout = np.array([vari_funcs.maximum_likelihood(fullfluxnorm[n,:], fullerrnorm[n,:], meanfull[n], posvar) for n in range(numobs)])

#%% Plot results ###
binupper = np.append(binedge[1:],np.max(z))
cbinupper = np.append(binedge[1:],np.max(chanz))
fbinupper = np.append(binedge[1:],np.max(fullz))
xlow = meanz-binedge
xhigh = binupper - meanz
cxlow = cmeanz-binedge
cxhigh = cbinupper - cmeanz
fxlow = fmeanz-binedge
fxhigh = fbinupper - fmeanz
plt.figure(figsize=[8,7])
plt.errorbar(z, out[:,0],yerr=out[:,1], fmt='o',color='b', alpha=0.2,
             label='Variable Non X-ray Source')
plt.errorbar(chanz, chanout[:,0],yerr=chanout[:,1], fmt='o', alpha=0.2,
             color='r',label='Variable X-ray Source')
plt.errorbar(fullz, fullout[:,0],yerr=fullout[:,1],fmt='o', alpha=0.2,
             color='tab:grey', label='Non variable X-ray source',zorder=1)
plt.errorbar(meanz, sig, xerr=[xlow, xhigh], yerr=sigerr, fmt='o',
             color='b', label='Ensemble Variable Sources')
plt.errorbar(cmeanz, csig, xerr=[cxlow, cxhigh], yerr=csigerr, fmt='o',
             color='r', label='Ensemble Variable X-ray Sources')
plt.errorbar(fmeanz, fsig, xerr=[fxlow, fxhigh], yerr=fsigerr, fmt='o',
             color='k', label='Ensemble Non Variable X-ray Sources')
#plt.xscale('log')
plt.xlabel('z')
plt.ylabel(r'$\sigma$')
plt.ylim(ymin=-0.25,ymax=1.755)
plt.xlim(xmin=-0.25, xmax=4.5)
plt.legend()





















end = time.time()
print(end-start)

