#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 29 12:13:52 2018

Code to look at how the fractional variability changes with flux

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
import matplotlib.colors as colors
from scipy.optimize import curve_fit
#from numpy.lib.recfunctions import append_fields
    
#%% Open the fits files and get data ###
chandata = fits.open('variable_tables/K/no06_variables_chi30_2arcsec_chandata_DR11data_restframe.fits')[1].data
xmmdata = fits.open('variable_tables/K/no06_variables_chi30_2arcsec_xmmdata_DR11data_restframe.fits')[1].data
fullxray = fits.open('mag_flux_tables/K/chanDR11data_restframe_mag_flux_table_best_extra_clean_no06.fits')[1].data
#fulluds = fits.open('mag_flux_tables/mag_flux_table_best_extra_clean_no06_DR11data.fits')[1].data
fulluds = fits.open('mag_flux_tables/K/mag_flux_table_extra_clean_no06_DR11data.fits')[1].data
sdata = fits.open('mag_flux_tables/K/stars_mag_flux_table_extra_clean_no06_DR11data.fits')[1].data
tbdata = fits.open('variable_tables/K/no06_variables_chi30_2arcsec_DR11data_restframe.fits')[1].data
sigtb = Table.read('sigma_tables/quad_epoch_sigma_table_extra_clean_no06_2arcsec.fits')

fulluds = fulluds[~fulluds['Stars-DR11']]
fulluds = fulluds[fulluds['Best']]

#posvar = np.linspace(0,2,5000)

#%% create function to do things 
def get_luminosity_and_flux(tbdata, xmm=False):
#    
    ### Extract magnitude table and error table ###
    flux = vari_funcs.k_mag_flux.flux4_stacks(tbdata)
    flux, tbdata = vari_funcs.flux_funcs.noneg(flux, tbdata)
#    tbdata = tbdata[np.nanmean(flux,axis=1)>1e4]
    flux, fluxerr, tbdata = vari_funcs.k_mag_flux.create_quad_error_array(sigtb, tbdata, aper=4)
    
    ### Normalise ###
    fluxnorm, fluxerrnorm = vari_funcs.flux_funcs.normalise_flux_and_errors(flux, fluxerr)
    
    ### Get Luminosity ###    
    L = tbdata['KFLUX_20']
    L[L < 0] = np.nan # remove those with null values
    mask = ~np.isnan(L)
    tbdata = tbdata[mask]
    L = L[mask]
    fluxnorm = fluxnorm[mask]
    fluxerrnorm = fluxerrnorm[mask]
    
    return tbdata, L, fluxnorm, fluxerrnorm


def run_max_likely(tbdata):
    
    ### Remove edges ###
    tbdata = vari_funcs.field_funcs.remove_edges(tbdata)
    
    ### Get luminosity and flux ###
    tbdata, flux, fluxnorm, fluxerrnorm = get_luminosity_and_flux(tbdata)
    
    ### read in max likelihood ###
    arr = np.load('np_arrays/fulluds_max_likelihood_values.npy')
    
    ### match the IDS in current tbdata to DR11 IDS in fulll uds ###
    mask = np.isin(arr[:,0], tbdata['NUMBER_1']) 
    out = np.empty([len(tbdata),2])
    out[:,0] = arr[mask,1]
    out[:,1] = arr[mask,2]
    
    return flux, out, tbdata
    
#%% run function ###
flux, out, tbdata = run_max_likely(tbdata)
chanflux, chanout, chandata = run_max_likely(chandata)
xmmflux, xmmout, xmmdata = run_max_likely(xmmdata)
fullflux, fullout, fullxray = run_max_likely(fullxray)
udsflux, udsout, udsdata = run_max_likely(fulluds)
sflux, sout, sdata = run_max_likely(sdata)


#%% Combine xray into single data set ###
xrayflux = np.append(chanflux, xmmflux)
xraysig = np.append(chanout[:,0],xmmout[:,0])
xraysigerr = np.append(chanout[:,1],xmmout[:,1])

#%% Plot results ###

#plt.figure(figsize=[12,7])
##plt.plot(fullL, fullout[:,0],'o', color='tab:grey', label='X-ray non variable', zorder=1)
##plt.errorbar(fullL, fullout[:,0], yerr=fullout[:,1], fmt='o', color='tab:grey', zorder=0, alpha=0.2)
##plt.plot(L, out[:,0], 'bo', label='Non-X-ray Variable', zorder=1)
#plt.plot(L, out[:,0], 'bo', label='Variable', zorder=1)
#plt.errorbar(L, out[:,0], yerr=out[:,1], fmt='ko', zorder=0, alpha=0.2)
##plt.plot(xrayL, xraysig, 'ro', label='X-ray variable', zorder=1)
##plt.errorbar(xrayL, xraysig, yerr=xraysigerr, fmt='ro', zorder=0, alpha=0.2)
##plt.plot(snL, snout[:,0], 'm*', markersize=10,mfc='None', label='Potential SN', zorder=3)
##plt.plot(radioL, radioout[:,0], 'rx', markersize=7,mfc='None', label='bad 07B', zorder=3)
##plt.plot(radioL, radioout[:,0], 'ks', markersize=10,mfc='None', label='Radio Source', zorder=3)
##plt.xscale('log')
#plt.xlabel('K Band Absolute Magnitude')
#plt.ylabel(r'$\sigma$')
##plt.ylim(ymin=-0.1,ymax=1.5)
##plt.xlim(xmin=1e32,xmax=2e40)
#plt.ylim(ymin=-0.05,ymax=1.5)
#plt.xlim(xmin=-29,xmax=-2)
##plt.ylim(ymin=-0.05,ymax=0.5)
##plt.xlim(xmin=-29,xmax=-18)
#plt.gca().invert_xaxis()
#plt.legend()

plt.figure(figsize=[12,7])
#plt.errorbar(flux, out[:,0], yerr=out[:,1], fmt='bo', zorder=2)
#plt.errorbar(udsflux, udsout[:,0], yerr=udsout[:,1], fmt='.', color='tab:grey', 
#             markersize=1, alpha=0.5, zorder=1)
plt.plot(flux, out[:,0], 'bo', zorder=3, label='Variable Non-X-ray AGN')
#plt.errorbar(chanflux, chanout[:,0], yerr=chanout[:,1], fmt='ro', zorder=2, label='Variable Chandra AGN')
#plt.plot(chanflux, chanout[:,0], 'ro', zorder=3, label='Variable Chandra AGN')
plt.plot(xrayflux, xraysig, 'ro', zorder=4, label='Variable X-ray AGN')
#plt.errorbar(fullflux, fullout[:,0], yerr=fullout[:,1], fmt='ko', zorder=1, label='Chandra AGN')
plt.plot(udsflux, udsout[:,0], '.', color='tab:grey', 
             markersize=1, zorder=1, label='Galaxy')
plt.plot(sflux, sout[:,0], 'm*', zorder=0, mfc='None', label='Star', markersize=10)
plt.plot(fullflux, fullout, 'ks', zorder=2, markersize=4,  label='Non-Variable X-ray AGN')
#plt.errorbar(fullL, fullout[:,0], yerr=fullout[:,1], fmt='ro',markersize=10, mfc='None', zorder=0, alpha=0.5)
#plt.plot(fullL, fullout[:,0], 'ro',markersize=10, mfc='None', zorder=0, alpha=0.5)
#plt.ylim(ymin=-0.05,ymax=1.4)
#plt.xlim(xmin=4e1,xmax=1e5)
plt.xscale('log')
#plt.yscale('log')
plt.xlabel(r'$K$-Band Flux')
plt.ylabel(r'$\sigma$')
#cbar = plt.colorbar()
plt.legend()
#cbar.set_label('z')
plt.tight_layout()

##%% find median of the full data set ###
#bins = np.array([13, 15])
#bins = np.append(bins, np.arange(16,24,0.2))
#bins = np.append(bins, [24])
#
#bins = 10**((30-bins)/2.5)
#bins = np.flip(bins, axis=0)
##bins = bins[16:44] #because of flux limit
#
#### Bin data ###
#allmedsig = np.array([])
#plt.figure(2, figsize=[12,7])
#for n, binedge in enumerate(bins):
##    print(binedge)
#    if n==np.size(bins)-1:
#        break
#    binflux, bindata = vari_funcs.fluxbin(binedge, bins[n+1], udsflux, tbdata) #bindata
#    
#    ### read in max likelihood ###
#    arr = np.load('fulluds_max_likelihood_values.npy')
#    
#    ### match the IDS in current tbdata to DR11 IDS in fulll uds ###
#    mask = np.isin(arr[:,0], bindata['NUMBER_1']) 
#    binsigs = arr[mask,1] ## get just the sig value
#    medsig = np.nanmedian(binsigs)
#    allmedsig = np.append(allmedsig, medsig)
#    
#    ### find difference from median for all points in bin
#    sig_diff = binsigs - medsig
#    sig_diff_sig = sig_diff/binsigs
#    
#    ### plot points ###
#    plt.plot(binflux, sig_diff_sig,'bo')
#    
#plt.figure(1)
#plt.plot(bins[0:42], allmedsig, 'k--')
#
#
##plt.figure(2, figsize=[12,7])
##plt.plot(flux, sig_diff_sig, 'bo')
#plt.xscale('log')
#plt.xlabel('K Band Flux')
#plt.ylabel(r'$\frac{\delta\sigma}{\sigma}$')
#cbar = plt.colorbar()
#plt.legend()
#cbar.set_label('z')
#plt.tight_layout()
#


end = time.time()
print(end-start)

