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
tbdata = fits.open('variable_tables/no06_variables_chi30_2arcsec_DR11data_restframe.fits')[1].data
#tbdata = fits.open('variable_tables/no06_variables_chi30_2arcsec_not_deviant_DR11data_restframe.fits')[1].data
#tbdata = fits.open('variable_tables/no06_variables_chi30_2arcsec_noXray_DR11data_restframe.fits')[1].data
chandata = fits.open('variable_tables/no06_variables_chi30_2arcsec_chandata_DR11data_restframe.fits')[1].data
xmmdata = fits.open('variable_tables/no06_variables_chi30_2arcsec_xmmdata_DR11data_restframe.fits')[1].data
fullxray = fits.open('mag_flux_tables/novarys_chanDR11data_restframe_mag_flux_table_best_extra_clean_no06.fits')[1].data
sigtb = Table.read('sigma_tables/quad_epoch_sigma_table_extra_clean_no06_2arcsec.fits')

posvar = np.linspace(0,2,5000)

#%% create function to do things ###

def get_luminosity_and_flux(tbdata, xmm=False):
    
    ### Extract magnitude table and error table ###
    flux = vari_funcs.flux4_stacks(tbdata)
    flux, tbdata = vari_funcs.noneg(flux, tbdata)
#    tbdata = tbdata[np.nanmean(flux,axis=1)>1e4]
    flux, fluxerr, tbdata = vari_funcs.create_quad_error_array(sigtb, tbdata, aper=4)
    
    ### Normalise ###
    fluxnorm, fluxerrnorm = vari_funcs.normalise_flux_and_errors(flux, fluxerr)
    
    ### Get Luminosity ###
    L = tbdata['M_K_z_p']
    L[L == 99] = np.nan
    mask = ~np.isnan(L)
    tbdata = tbdata[mask]
    L = L[mask]
    fluxnorm = fluxnorm[mask]
    fluxerrnorm = fluxerrnorm[mask]
    
    return tbdata, L, fluxnorm, fluxerrnorm

def make_ensemble(tbdata, L, binedge, fluxnorm, fluxerrnorm):
    
    ### Create dicts to save data into ###
    enflux = {}
    enfluxerr = {}
    enL = {}
    
    ### loop over bins ###
    for m, enmin in enumerate(binedge):
        ### Isolate data needed ###
        mask1 = L >= enmin
        if m != len(binedge)-1:
            mask2 = L < binedge[m+1]
        else:
            mask2 = np.ones(len(mask1))
        enmask = mask1*mask2.astype(bool)
        enflux[m] = fluxnorm[enmask]
        enfluxerr[m] = fluxerrnorm[enmask]
        enL[m] = L[enmask]
        
    return enflux, enfluxerr, enL

def get_ensemble_sig(enflux, enfluxerr, enL, posvar):
    ### Create dicts to save data into ###
    size = len(enflux)
    sig = np.empty(size)
    sigerr = np.empty(size)
    meanL = np.empty(size)
    for m in enflux:
        ### Combine into one flux curve per bin ###
        enfluxcurve = np.ravel(enflux[m])
        enfluxcurveerr = np.ravel(enfluxerr[m])
        
        ### Find max likelihood sig of curve ###
        [sig[m],sigerr[m]] = vari_funcs.maximum_likelihood(enfluxcurve, enfluxcurveerr, 1, posvar)

        ### find mean L ###
        meanL[m] = np.nanmean(enL[m])
    
    return sig, sigerr, meanL
        

def get_binedge(binsize, sortedL):
    binedge = np.empty(int(len(sortedL)/binsize))
    size = len(binedge)
    for n in range(size):
        ### Define bin edges ###
        enL = sortedL[n*binsize:(n*binsize)+binsize]
        enmin = np.nanmin(enL)
        binedge[n] = enmin
    return binedge
    
def combine_Xray(binedge, xmmenflux, xmmenfluxerr, xmmenL, 
                 chanenflux, chanenfluxerr, chanenL):
    enflux = {}
    enfluxerr = {}
    enL = {}
    for m, enmin in enumerate(binedge):
        enflux[m] = np.vstack([xmmenflux[m],chanenflux[m]])
        enfluxerr[m] = np.vstack([xmmenfluxerr[m],chanenfluxerr[m]])
        enL[m] = np.append(xmmenL[m],chanenL[m])
    return enflux, enfluxerr, enL

#%% remove edges ###
    
tbdata = vari_funcs.remove_edges(tbdata)
chandata = vari_funcs.remove_edges(chandata)
xmmdata = vari_funcs.remove_edges(xmmdata)
fullxray = vari_funcs.remove_edges(fullxray)

#%% get_luminosity_and_flux ###

tbdata, L, fluxnorm, fluxerrnorm= get_luminosity_and_flux(tbdata)
chandata, chanL, chanfluxnorm, chanerrnorm= get_luminosity_and_flux(chandata)
xmmdata, xmmL, xmmfluxnorm, xmmerrnorm = get_luminosity_and_flux(xmmdata, xmm=True)
fullxray, fullL, fullfluxnorm, fullerrnorm= get_luminosity_and_flux(fullxray)

#%% Split into bins according to nir luminosity ###

### Join xmm and chan arrays ###
xrayL = np.append(xmmL, chanL)

### Sort arrays so can find bin edges ###
sortedL = np.sort(L)
sortedxrayL = np.sort(xrayL)
sortedfullL = np.sort(fullL)

### define bin size and find edge arrays ###
binsize = 8
binedge = get_binedge(binsize, sortedL)
xraybinedge = get_binedge(binsize, sortedxrayL)
fullbinedge = get_binedge(binsize, sortedfullL)

#%% create ensemble arrays ###

enflux, enfluxerr, enL = make_ensemble(tbdata, L, binedge, 
                                       fluxnorm, fluxerrnorm)
chanenflux, chanenfluxerr, chanenL = make_ensemble(chandata, chanL, xraybinedge, 
                                                   chanfluxnorm, chanerrnorm)
xmmenflux, xmmenfluxerr, xmmenL = make_ensemble(xmmdata, xmmL, xraybinedge,
                                                xmmfluxnorm, xmmerrnorm )
fullenflux, fullenfluxerr, fullenL = make_ensemble(fullxray, fullL, fullbinedge, 
                                                   fullfluxnorm,fullerrnorm)

#%% Combine chan and xmm ensembles ###

xrayenflux, xrayenfluxerr, xrayenL = combine_Xray(xraybinedge, xmmenflux,
                                                  xmmenfluxerr, xmmenL, 
                                                  chanenflux, chanenfluxerr, 
                                                  chanenL)

#%% Get ensemble maximum likelihood ###
sig, sigerr, meanL = get_ensemble_sig(enflux, enfluxerr, enL, posvar)
xraysig, xraysigerr, xraymeanL = get_ensemble_sig(xrayenflux, xrayenfluxerr, 
                                                  xrayenL, posvar)
fullsig, fullsigerr, fullmeanL = get_ensemble_sig(fullenflux, fullenfluxerr, 
                                                  fullenL, posvar)

#%% Plot results ###
def find_err_values(binedge, sortedL, meanL):
    binupper = np.append(binedge[1:],np.max(sortedL))
    xlow = meanL-binedge
    xhigh = binupper - meanL
    return xlow, xhigh
    
xerr = find_err_values(binedge, sortedL, meanL)
xrayxerr = find_err_values(xraybinedge, sortedxrayL, xraymeanL)
fullxerr = find_err_values(fullbinedge, sortedfullL, fullmeanL)


plt.figure(figsize=[12,7])
plt.errorbar(meanL, sig, xerr=xerr, yerr=sigerr, fmt='o',
             color='b', label='Variable')
#plt.errorbar(meanL, sig, xerr=xerr, yerr=sigerr, fmt='o',
#             color='b', label='Variable Non-X-ray')
#plt.errorbar(xraymeanL, xraysig, xerr=xrayxerr, yerr=xraysigerr, fmt='o',
#             color='r', label='Variable X-ray')
plt.errorbar(fullmeanL, fullsig, xerr=fullxerr, 
             yerr=fullsigerr, fmt='o',color='tab:grey', label='Non Variable X-ray')
plt.xlabel('K Band Absolute Magnitude')
plt.ylabel(r'$\sigma$')
plt.ylim(ymin=-0.05,ymax=1.5)
plt.xlim(xmin=-29,xmax=-2)
#plt.ylim(ymin=-0.05,ymax=0.5)
#plt.xlim(xmin=-29,xmax=-18)
plt.gca().invert_xaxis()
plt.legend()






















end = time.time()
print(end-start)

