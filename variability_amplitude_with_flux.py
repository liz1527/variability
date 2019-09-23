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

### Define cosmology ###
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    
#%% Open the fits files and get data ###
chandata = fits.open('variable_tables/no06_variables_chi30_2arcsec_chandata_DR11data_restframe.fits')[1].data
xmmdata = fits.open('variable_tables/no06_variables_chi30_2arcsec_xmmdata_DR11data_restframe.fits')[1].data
fullxray = fits.open('mag_flux_tables/novarys_chanDR11data_restframe_mag_flux_table_best_extra_clean_no06.fits')[1].data
fulluds = fits.open('mag_flux_tables/mag_flux_table_extra_clean_no06_DR11data.fits')[1].data
tbdata = fits.open('variable_tables/no06_variables_chi30_2arcsec_noXray_DR11data_restframe.fits')[1].data
sigtb = Table.read('sigma_tables/quad_epoch_sigma_table_extra_clean_no06_2arcsec.fits')

posvar = np.linspace(0,2,5000)

#%% create function to do things 
def get_luminosity_and_flux(tbdata, xmm=False):
#    
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


def run_max_likely(tbdata):
    posvar = np.linspace(0,2,5000)
    ### Remove edges ###
    tbdata = vari_funcs.remove_edges(tbdata)
    
    ### Get luminosity and flux ###
    tbdata, L, fluxnorm, fluxerrnorm= get_luminosity_and_flux(tbdata)
    
    ### Get sig values ###
    numobs = np.shape(fluxnorm)[0]
    meanflux = np.nanmean(fluxnorm, axis=1)
    out = np.array([vari_funcs.maximum_likelihood(fluxnorm[n,:], 
                                                  fluxerrnorm[n,:], meanflux[n], 
                                                  posvar, n=n, printn=100) for n in range(numobs)])
    return L, out, tbdata
    
#%% run function ###
L, out, tbdata = run_max_likely(tbdata)
chanL, chanout, chandata = run_max_likely(chandata)
xmmL, xmmout, xmmdata = run_max_likely(xmmdata)
fullL, fullout, fullxray = run_max_likely(fullxray)
udsL, udsout, udsdata = run_max_likely(fulluds)


#%% Combine xray into single data set ###
xrayL = np.append(chanL, xmmL)
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
#plt.errorbar(L, out[:,0], yerr=out[:,1], fmt='bo', zorder=2)
#plt.errorbar(udsL, udsout[:,0], yerr=udsout[:,1], fmt='ko', zorder=1)
plt.plot(L, out[:,0], 'bo', zorder=2)
plt.plot(udsL, udsout[:,0], 'ko', zorder=1)
#plt.errorbar(fullL, fullout[:,0], yerr=fullout[:,1], fmt='ro',markersize=10, mfc='None', zorder=0, alpha=0.5)
#plt.plot(fullL, fullout[:,0], 'ro',markersize=10, mfc='None', zorder=0, alpha=0.5)
plt.ylim(ymin=-0.05,ymax=1.5)
#plt.xlim(xmin=4e1,xmax=1e5)
plt.xscale('log')
plt.xlabel('K Band Flux')
plt.ylabel(r'$\sigma$')
#cbar = plt.colorbar()
plt.legend()
#cbar.set_label('z')
plt.tight_layout()

#%% Save maximum likelihood in a np array so it can be read in, not calculated ###
ids = udsdata['NUMBER_1']
arr = np.array([ids, udsout[:,0], udsout[:,1]])
arr = arr.T
np.save('fulluds_max_likelihood_values', arr)

#%% curve fit to faint end ###
def func_pol(x,a,b):
    return a*+b*(x**-1)
def func_exp(x,a,b):
#    return a*np.exp(b*x)
    return a*(x**(b))

Lfit1 = L<1e3
Lfit2 = L>1e2
Lfit = Lfit1*Lfit2.astype(bool)
popt, pcov = curve_fit(func_pol, L[Lfit], out[Lfit,0], sigma=out[Lfit,1])
x = np.logspace(2,3)
y = func_pol(x, popt[0], popt[1])
popt2, pcov2 = curve_fit(func_exp, L[Lfit], out[Lfit,0], sigma=out[Lfit,1])
y2 = func_exp(x, popt2[0], popt2[1])

perr = np.sqrt(np.diag(pcov))
perr2 = np.sqrt(np.diag(pcov2))
#yerr1 = func_pol(x, popt[0]+perr[0], popt[1]-perr[1])
#yerr2 = func_pol(x, popt[0]-perr[0], popt[1]+perr[1])
yerr1 = func_exp(x, popt2[0]+perr2[0], popt2[1]-perr2[1])
yerr2 = func_exp(x, popt2[0]-perr2[0], popt2[1]+perr2[1])
#plt.plot(x,y)
plt.plot(x,y2)
#plt.plot(x,yerr1,'r--')
#plt.plot(x,yerr2,'r--')

#%% find distance from that fit for each point ###
sig_fit = func_exp(L[Lfit], popt2[0]+perr2[0], popt2[1]-perr2[1])
sig_diff = out[Lfit,0] - sig_fit
sig_diff_sig = sig_diff/out[Lfit,0]

plt.figure(figsize=[12,7])
plt.plot(L[Lfit], sig_diff_sig, 'bo')
plt.xscale('log')
plt.xlabel('K Band Flux')
plt.ylabel(r'$\frac{\delta\sigma}{\sigma}$')
cbar = plt.colorbar()
plt.legend()
cbar.set_label('z')
plt.tight_layout()



end = time.time()
print(end-start)

