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
#chandata = fits.open('variable_tables/no06_variables_chi30_2arcsec_chandata_DR11data_restframe.fits')[1].data
#xmmdata = fits.open('variable_tables/no06_variables_chi30_2arcsec_xmmdata_DR11data_restframe.fits')[1].data
#fullxray = fits.open('mag_flux_tables/novarys_chanDR11data_restframe_mag_flux_table_best_extra_clean_no06.fits')[1].data
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
    L = tbdata['KFLUX_20']
    L[L < 0] = np.nan # remove those with null values
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
        
def z_split(tbdata):
    z = tbdata['z_spec']#[mask]
    z[z==-1] = tbdata['z_p'][z==-1]
    tb1 = tbdata[z <= 0.5]
    mask = np.array(z>0.5)*np.array(z<=1).astype('bool')
    tb2 = tbdata[mask]
    mask = np.array(z>1)*np.array(z<=1.4).astype('bool')
    tb3 = tbdata[mask]
    mask = np.array(z>1.5)*np.array(z<=2).astype('bool')
    tb4 = tbdata[mask]
    tb5 = tbdata[z > 2]
#    tb2 = tbdata[z>2]
    return tb1, tb2, tb3, tb4, tb5


def get_binedge(binsize, sortedL):
    binedge = np.empty(int(len(sortedL)/binsize))
    size = len(binedge)
    for n in range(size):
        ### Define bin edges ###
        enL = sortedL[n*binsize:(n*binsize)+binsize]
        enmin = np.nanmin(enL)
        binedge[n] = enmin
    return binedge
    
def find_err_values(binedge, sortedL, meanL):
    binupper = np.append(binedge[1:],np.max(sortedL))
    xlow = meanL-binedge
    xhigh = binupper - meanL
    return xlow, xhigh

def get_ensemble_data(tbdata):
    
    ### get_luminosity_and_flux ###
    tbdata, L, fluxnorm, fluxerrnorm= get_luminosity_and_flux(tbdata)
    
    ### Split into bins according to nir luminosity ###
    ### Sort arrays so can find bin edges ###
    sortedL = np.sort(L)
    
    ### define bin size and find edge arrays ###
    binsize = 8
    binedge = get_binedge(binsize, sortedL)
    
    ### create ensemble arrays ###
    
    enflux, enfluxerr, enL = make_ensemble(tbdata, L, binedge,
                                                       fluxnorm, fluxerrnorm)
    
    ### Get ensemble maximum likelihood ###
    sig, sigerr, meanL = get_ensemble_sig(enflux, enfluxerr, enL, posvar)
    
    #### get errors ###
    xerr = find_err_values(binedge, sortedL, meanL)
    
    return tbdata, meanL, xerr, sig, sigerr
    
#%% remove edges and split by z ###
tbdata = vari_funcs.remove_edges(tbdata)
#fullxray = vari_funcs.remove_edges(fullxray)

tb1, tb2, tb3, tb4, tb5 = z_split(tbdata)
#fullhigh, fulllow = z_split(fullxray)

#%% run function ###
tb1, meanL1, xerr1, sig1, sigerr1 = get_ensemble_data(tb1)
tb2, meanL2, xerr2, sig2, sigerr2 = get_ensemble_data(tb2)
tb3, meanL3, xerr3, sig3, sigerr3 = get_ensemble_data(tb3)
tb4, meanL4, xerr4, sig4, sigerr4 = get_ensemble_data(tb4)
tb5, meanL5, xerr5, sig5, sigerr5 = get_ensemble_data(tb5)

#%% run plot ###
plt.figure(figsize=[12,7])
#plt.subplot(121)
plt.errorbar(meanL1, sig1, xerr=xerr1, yerr=sigerr1, fmt='o',
             color='g', label='0 < z <= 0.5')
plt.errorbar(meanL2, sig2, xerr=xerr2, yerr=sigerr2, fmt='o',
             color='y', label='0.5 < z <= 1')
plt.errorbar(meanL3, sig3, xerr=xerr3, yerr=sigerr3, fmt='o',
             color='b', label='1 < z <= 1.5')
plt.errorbar(meanL4, sig4, xerr=xerr4, yerr=sigerr4, fmt='o',
             color='m', label='1.5 < z <= 2')
plt.errorbar(meanL5, sig5, xerr=xerr5, yerr=sigerr5, fmt='o',
             color='r', label='z > 2')
plt.xlabel('K Band Flux')
plt.ylabel(r'$\sigma$')
plt.ylim(ymin=-0.05,ymax=1.3)
plt.xlim(xmin=4e1,xmax=1e5)
#plt.ylim(ymin=-0.05,ymax=0.5)
#plt.xlim(xmin=-29,xmax=-18)
plt.xscale('log')
#plt.gca().invert_xaxis()
plt.legend()
plt.tight_layout()





















end = time.time()
print(end-start)

