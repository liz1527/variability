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
#def remove_edges(tbdata):
#    
#    ### Set X limits ###
#    x = tbdata['X_IMAGE_05B']
#    xmask1 = x > 1000
#    xmask2 = x < 24000
#    xmask = xmask1 * xmask2.astype(bool)
#    tbdata = tbdata[xmask]
#    
#    ### Set Y limits ###
#    y = tbdata['Y_IMAGE_05B']
#    ymask1 = y > 1000
#    ymask2 = y < 24000
#    ymask = ymask1 * ymask2.astype(bool)
#    tbdata = tbdata[ymask]
#    
#    return tbdata
    
#%% Open the fits files and get data ###
chandata = fits.open('variable_tables/no06_variables_chi30_chandata_DR11data_restframe.fits')[1].data
xmmdata = fits.open('variable_tables/no06_variables_chi30_xmmdata_DR11data_restframe.fits')[1].data
fullxray = fits.open('mag_flux_tables/novarys_chanDR11data_restframe_mag_flux_table_best_extra_clean_no06.fits')[1].data
tbdata = fits.open('variable_tables/no06_variables_chi30_DR11data_restframe.fits')[1].data
radiodata = fits.open('variable_tables/no06_variables_chi30_radiodata_DR11data_restframe.fits')[1].data
sigtb = Table.read('sigma_tables/quad_epoch_sigma_table_extra_clean_no06.fits')

posvar = np.linspace(0,2,5000)

#%% create function to do things 
def get_luminosity_and_flux(tbdata, xmm=False):
#    
    ### Extract magnitude table and error table ###
    flux = vari_funcs.flux5_stacks(tbdata)
    flux, tbdata = vari_funcs.noneg(flux, tbdata)
    tbdata = tbdata[np.nanmean(flux,axis=1)>1e4]
    flux, fluxerr, tbdata = vari_funcs.create_quad_error_array(sigtb, tbdata)
    
    ### Normalise ###
    fluxnorm, fluxerrnorm = vari_funcs.normalise_flux_and_errors(flux, fluxerr)
    
    ### Get absolute magnitude ###
    L = tbdata['M_K_z_p']
    
    return tbdata, L, fluxnorm, fluxerrnorm

def z_split(tbdata):
    z = tbdata['z_spec']#[mask]
    z[z==-1] = tbdata['z_p'][z==-1]
    tbhigh = tbdata[z > 1.6]
    tblow = tbdata[z <= 1.6]
    return tbhigh, tblow
#%% Remove edges ###
tbdata = vari_funcs.remove_edges(tbdata)
chandata = vari_funcs.remove_edges(chandata)
xmmdata = vari_funcs.remove_edges(xmmdata)
fullxray = vari_funcs.remove_edges(fullxray)
radiodata = vari_funcs.remove_edges(radiodata)

#%% Split by z ###
tbhigh, tblow = z_split(tbdata)
xmmhigh, xmmlow = z_split(xmmdata)
chanhigh, chanlow = z_split(chandata)
fullhigh, fulllow = z_split(fullxray)
radiohigh, radiolow = z_split(radiodata)

#xmmz = xmmdata['z_spec']#[mask]
#xmmz[xmmz==-1] = xmmdata['z_p'][xmmz==-1]
#xmmhigh = xmmdata[xmmz > 1.6]
#xmmlow = xmmdata[xmmz <= 1.6]
#chanz = chandata['z_spec']#[mask]
#chanz[chanz==-1] = chandata['z_p'][chanz==-1]
#chanhigh = chandata[chanz > 1.6]
#chanlow = chandata[chanz <= 1.6]
#fullz = fullxray['z_spec']#[mask]
#fullz[fullz==-1] = fullxray['z_p'][fullz==-1]
#fullhigh = fullxray[fullz > 1.6]
#fulllow = fullxray[fullz <= 1.6]

#%%
high, highL, highfluxnorm, higherrnorm= get_luminosity_and_flux(tbhigh)
low, lowL, lowfluxnorm, lowerrnorm= get_luminosity_and_flux(tblow)
fullhigh, fullhighL, fullhighfluxnorm, fullhigherrnorm= get_luminosity_and_flux(fullhigh)
fulllow, fulllowL, fulllowfluxnorm, fulllowerrnorm= get_luminosity_and_flux(fulllow)
chanhigh, highchanL, highchanfluxnorm, highchanerrnorm= get_luminosity_and_flux(chanhigh)
chanlow, lowchanL, lowchanfluxnorm, lowchanerrnorm= get_luminosity_and_flux(chanlow)
xmmhigh, highxmmL, highxmmfluxnorm, highxmmerrnorm = get_luminosity_and_flux(xmmhigh, xmm=True)
xmmlow, lowxmmL, lowxmmfluxnorm, lowxmmerrnorm = get_luminosity_and_flux(xmmlow, xmm=True)
radiohigh, highradioL, highradiofluxnorm, highradioerrnorm= get_luminosity_and_flux(radiohigh)
radiolow, lowradioL, lowradiofluxnorm, lowradioerrnorm= get_luminosity_and_flux(radiolow)

#%% get sig values
numobs = np.shape(highfluxnorm)[0]
meanflux = np.nanmean(highfluxnorm, axis=1)
highout = np.array([vari_funcs.maximum_likelihood(highfluxnorm[n,:], 
                                              higherrnorm[n,:], meanflux[n], 
                                              posvar) for n in range(numobs)])

numobs = np.shape(lowfluxnorm)[0]
meanflux = np.nanmean(lowfluxnorm, axis=1)
lowout = np.array([vari_funcs.maximum_likelihood(lowfluxnorm[n,:], 
                                                  lowerrnorm[n,:], meanflux[n], 
                                                  posvar) for n in range(numobs)])
numobs = np.shape(fullhighfluxnorm)[0]
meanflux = np.nanmean(fullhighfluxnorm, axis=1)
fullhighout = np.array([vari_funcs.maximum_likelihood(fullhighfluxnorm[n,:], 
                                              fullhigherrnorm[n,:], meanflux[n], 
                                              posvar) for n in range(numobs)])

numobs = np.shape(fulllowfluxnorm)[0]
meanflux = np.nanmean(fulllowfluxnorm, axis=1)
fulllowout = np.array([vari_funcs.maximum_likelihood(fulllowfluxnorm[n,:], 
                                                  fulllowerrnorm[n,:], meanflux[n], 
                                                  posvar) for n in range(numobs)])
numobs = np.shape(highchanfluxnorm)[0]
meanflux = np.nanmean(highchanfluxnorm, axis=1)
highchanout = np.array([vari_funcs.maximum_likelihood(highchanfluxnorm[n,:], 
                                              highchanerrnorm[n,:], meanflux[n], 
                                              posvar) for n in range(numobs)])

numobs = np.shape(lowchanfluxnorm)[0]
meanchan = np.nanmean(lowchanfluxnorm, axis=1)
lowchanout = np.array([vari_funcs.maximum_likelihood(lowchanfluxnorm[n,:], 
                                                  lowchanerrnorm[n,:], meanchan[n], 
                                                  posvar) for n in range(numobs)])

numobs = np.shape(highxmmfluxnorm)[0]
meanfull = np.nanmean(highxmmfluxnorm, axis=1)
highxmmout = np.array([vari_funcs.maximum_likelihood(highxmmfluxnorm[n,:], 
                                                  highxmmerrnorm[n,:], meanfull[n], 
                                                  posvar) for n in range(numobs)])

numobs = np.shape(lowxmmfluxnorm)[0]
meanfull = np.nanmean(lowxmmfluxnorm, axis=1)
lowxmmout = np.array([vari_funcs.maximum_likelihood(lowxmmfluxnorm[n,:], 
                                                  lowxmmerrnorm[n,:], meanfull[n], 
                                                  posvar) for n in range(numobs)])
#numobs = np.shape(highradiofluxnorm)[0]
#meanflux = np.nanmean(highradiofluxnorm, axis=1)
#highradioout = np.array([vari_funcs.maximum_likelihood(highradiofluxnorm[n,:], 
#                                              highradioerrnorm[n,:], meanflux[n], 
#                                              posvar) for n in range(numobs)])
#
#numobs = np.shape(lowradiofluxnorm)[0]
#meanflux = np.nanmean(lowradiofluxnorm, axis=1)
#lowradioout = np.array([vari_funcs.maximum_likelihood(lowradiofluxnorm[n,:], 
#                                                  lowradioerrnorm[n,:], meanflux[n], 
#                                                  posvar) for n in range(numobs)])

#%% Combine into single data set ###
hL = np.append(highchanL, highxmmL)
hsig = np.append(highchanout[:,0],highxmmout[:,0])
hsigerr = np.append(highchanout[:,1],highxmmout[:,1])
lL = np.append(lowchanL, lowxmmL)
lsig = np.append(lowchanout[:,0], lowxmmout[:,0])
lsigerr = np.append(lowchanout[:,1], lowxmmout[:,1])

#%% Plot results ###
plt.figure(figsize=[12,7])
plt.subplot(121)
plt.plot(fulllowL, fulllowout[:,0],'o', color='tab:grey', label='X-ray non variable', zorder=1)
plt.errorbar(fulllowL, fulllowout[:,0], yerr=fulllowout[:,1], fmt='o', color='tab:grey', zorder=0, alpha=0.2)
plt.plot(lowL, lowout[:,0], 'bo', label='z <= 1.6', zorder=1)
plt.errorbar(lowL, lowout[:,0], yerr=lowout[:,1], fmt='bo', zorder=0, alpha=0.2)
#plt.plot(lL, lsig, 'ro', label='X-ray variable', zorder=1)
#plt.errorbar(lL, lsig, yerr=lsigerr, fmt='ro', zorder=0, alpha=0.2)
#plt.plot(lowradioL, lowradioout[:,0], 'ks', markersize=10,mfc='None', label='Radio Source', zorder=3)
#plt.xscale('log')
plt.xlabel('K Band Absolute Magnitude')
plt.ylabel(r'$\sigma$')
#plt.ylim(ymin=-0.1,ymax=1.5)
#plt.xlim(xmin=-29,xmax=-12)
plt.ylim(ymin=-0.05,ymax=0.35)
plt.xlim(xmin=-29,xmax=-15)
plt.gca().invert_xaxis()
plt.legend()

plt.subplot(122)
plt.plot(fullhighL, fullhighout[:,0], 'o', color='tab:grey', label='X-ray non variable', zorder=1)
plt.errorbar(fullhighL, fullhighout[:,0], yerr=fullhighout[:,1], fmt='o', color='tab:grey', zorder=0, alpha=0.2)
plt.plot(highL, highout[:,0], 'go', label='z > 1.6', zorder=1)
plt.errorbar(highL, highout[:,0], yerr=highout[:,1], fmt='go', zorder=0, alpha=0.2)
#plt.plot(hL, hsig, 'ro', label='X-ray variable',zorder=2)
#plt.errorbar(hL, hsig, yerr=hsigerr, fmt='ro', zorder=0, alpha=0.2)
#plt.plot(highradioL, highradioout[:,0], 'ks', markersize=10,mfc='None', label='Radio Source', zorder=3)
#plt.xscale('log')
plt.xlabel('K Band Absolute Magnitude')
plt.ylabel(r'$\sigma$')
#plt.ylim(ymin=-0.1,ymax=1.5)
#plt.xlim(xmin=-29,xmax=-12)
plt.ylim(ymin=-0.05,ymax=0.35)
plt.xlim(xmin=-29,xmax=-15)
plt.gca().invert_xaxis()
plt.legend()

plt.tight_layout()
#
#plt.figure(figsize=[12,7])
#plt.subplot(121)
#plt.plot(lL, lsig, 'bo', label='z <= 1.6', zorder=1)
#plt.errorbar(lL, lsig, yerr=lsigerr, fmt='bo', zorder=0, alpha=0.2)
#plt.xscale('log')
#plt.xlabel('K Band Luminosity (W)')
#plt.ylabel(r'$\sigma$')
#plt.ylim(ymin=0,ymax=0.25)
#plt.xlim(xmin=1e37,xmax=1e40)
#plt.legend()
#
#plt.subplot(122)
#plt.plot(hL, hsig, 'go', label='z > 1.6',zorder=1)
#plt.errorbar(hL, hsig, yerr=hsigerr, fmt='go', zorder=0, alpha=0.2)
#plt.xscale('log')
#plt.xlabel('K Band Luminosity (W)')
#plt.ylabel(r'$\sigma$')
#plt.ylim(ymin=0,ymax=0.25)
#plt.xlim(xmin=1e37,xmax=1e40)
#plt.legend()
#plt.tight_layout()


end = time.time()
print(end-start)

