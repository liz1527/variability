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
chandata = fits.open('variable_tables/no06_variables_chi30_2arcsec_chandata_DR11data_restframe.fits')[1].data
xmmdata = fits.open('variable_tables/no06_variables_chi30_2arcsec_xmmdata_DR11data_restframe.fits')[1].data
fullxray = fits.open('mag_flux_tables/novarys_chanDR11data_restframe_mag_flux_table_best_extra_clean_no06.fits')[1].data
tbdata = fits.open('variable_tables/no06_variables_chi30_2arcsec_not_deviant_DR11data_restframe.fits')[1].data
sndata = fits.open('variable_tables/no06_variables_chi30_2arcsec_DR11data_restframe_SN.fits')[1].data
radiodata = fits.open('variable_tables/no06_variables_chi30_2arcsec_DR11data_restframe_07B.fits')[1].data
#radiodata = fits.open('variable_tables/no06_variables_chi30_radiodata_DR11data_restframe.fits')[1].data
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
    
    ### Find luminosity distance ###
    z = tbdata['z_spec']#[mask]
    z[z==-1] = tbdata['z_p'][z==-1]
    DL = cosmo.luminosity_distance(z)
    DL = DL.to(u.m)
    
    ### Get AB magnitude ###
#    abmag = tbdata['KMAG_20']
    
    ### Convert to luminosity using formula worked out in lab book ###
#    L = 10 ** (-abmag/2.5) * 3631 * 1e-26 * (3e8/0.34e-6) * 4 * np.pi * (DL.value**2) # gives luminosity in W
    L = tbdata['M_K_z_p']
    
    return tbdata, L, fluxnorm, fluxerrnorm

def z_split(tbdata):
    z = tbdata['z_spec']#[mask]
    z[z==-1] = tbdata['z_p'][z==-1]
    tbhigh = tbdata[z > 1.6]
    tblow = tbdata[z <= 1.6]
    return tbhigh, tblow

def z_split_max_likely(tbdata):
    posvar = np.linspace(0,2,5000)
    ### Remove edges ###
    tbdata = vari_funcs.remove_edges(tbdata)
    
    ### Split by z ###
    tbhigh, tblow = z_split(tbdata)
    
    ### Get luminosity and flux ###
    tbhigh, highL, highfluxnorm, higherrnorm= get_luminosity_and_flux(tbhigh)
    tblow, lowL, lowfluxnorm, lowerrnorm= get_luminosity_and_flux(tblow)
    
    ### Get sig values ###
    numobs = np.shape(highfluxnorm)[0]
    meanflux = np.nanmean(highfluxnorm, axis=1)
    highout = np.array([vari_funcs.maximum_likelihood(highfluxnorm[n,:], 
                                                  higherrnorm[n,:], meanflux[n], 
                                                  posvar, n=n, printn=100) for n in range(numobs)])
    
    numobs = np.shape(lowfluxnorm)[0]
    meanflux = np.nanmean(lowfluxnorm, axis=1)
    lowout = np.array([vari_funcs.maximum_likelihood(lowfluxnorm[n,:], 
                                                      lowerrnorm[n,:], meanflux[n], 
                                                      posvar, n=n, printn=100) for n in range(numobs)])
    return highL, highout, lowL, lowout, tbdata
    
#%% run function ###
highL, highout, lowL, lowout, tbdata = z_split_max_likely(tbdata)
highchanL, highchanout, lowchanL, lowchanout, chandata = z_split_max_likely(chandata)
highxmmL, highxmmout, lowxmmL, lowxmmout, xmmdata = z_split_max_likely(xmmdata)
highfullL, highfullout, lowfullL, lowfullout, fullxray = z_split_max_likely(fullxray)
highsnL, highsnout, lowsnL, lowsnout, sndata = z_split_max_likely(sndata)
highradioL, highradioout, lowradioL, lowradioout, radiodata = z_split_max_likely(radiodata)

#%% Combine xray into single data set ###
hL = np.append(highchanL, highxmmL)
hsig = np.append(highchanout[:,0],highxmmout[:,0])
hsigerr = np.append(highchanout[:,1],highxmmout[:,1])
lL = np.append(lowchanL, lowxmmL)
lsig = np.append(lowchanout[:,0], lowxmmout[:,0])
lsigerr = np.append(lowchanout[:,1], lowxmmout[:,1])

#%% Plot results ###
plt.figure(figsize=[12,7])
#plt.subplot(121)
#plt.plot(lowfullL, lowfullout[:,0],'o', color='tab:grey', label='X-ray non variable', zorder=1)
#plt.errorbar(lowfullL, lowfullout[:,0], yerr=lowfullout[:,1], fmt='o', color='tab:grey', zorder=0, alpha=0.2)
plt.plot(lowL, lowout[:,0], 'bo', label='z <= 1.6', zorder=1)
plt.errorbar(lowL, lowout[:,0], yerr=lowout[:,1], fmt='bo', zorder=0, alpha=0.2)
#plt.plot(lL, lsig, 'ro', label='X-ray variable', zorder=1)
#plt.errorbar(lL, lsig, yerr=lsigerr, fmt='ro', zorder=0, alpha=0.2)
#plt.plot(lowsnL, lowsnout[:,0], 'm*', markersize=10,mfc='None', label='Potential SN', zorder=3)
#plt.plot(lowradioL, lowradioout[:,0], 'rx', markersize=7,mfc='None', label='bad 07B', zorder=3)
#plt.plot(lowradioL, lowradioout[:,0], 'ks', markersize=10,mfc='None', label='Radio Source', zorder=3)
#plt.xscale('log')
plt.xlabel('K Band Absolute Magnitude')
plt.ylabel(r'$\sigma$')
#plt.ylim(ymin=-0.1,ymax=1.5)
#plt.xlim(xmin=1e32,xmax=2e40)
plt.ylim(ymin=-0.05,ymax=1.3)
plt.xlim(xmin=-29,xmax=-2)
#plt.ylim(ymin=-0.05,ymax=0.5)
#plt.xlim(xmin=-29,xmax=-18)
plt.gca().invert_xaxis()
plt.legend()

#plt.subplot(122)
#plt.plot(highfullL, highfullout[:,0], 'o', color='tab:grey', label='X-ray non variable', zorder=1)
#plt.errorbar(highfullL, highfullout[:,0], yerr=highfullout[:,1], fmt='o', color='tab:grey', zorder=0, alpha=0.2)
plt.plot(highL, highout[:,0], 'go', label='z > 1.6', zorder=1)
plt.errorbar(highL, highout[:,0], yerr=highout[:,1], fmt='go', zorder=0, alpha=0.2)
#plt.plot(hL, hsig, 'ro', label='X-ray variable',zorder=2)
#plt.errorbar(hL, hsig, yerr=hsigerr, fmt='ro', zorder=0, alpha=0.2)
#plt.plot(highsnL, highsnout[:,0], 'm*', markersize=10,mfc='None', label='Potential SN', zorder=3)
#plt.plot(highradioL, highradioout[:,0], 'rx', markersize=7,mfc='None', label='bad 07B', zorder=3)
#plt.plot(highradioL, highradioout[:,0], 'ks', markersize=10,mfc='None', label='Radio Source', zorder=3)
#plt.xscale('log')
plt.xlabel('K Band Absolute Magnitude')
plt.ylabel(r'$\sigma$')
#plt.ylim(ymin=-0.1,ymax=1.5)
#plt.xlim(xmin=1e32,xmax=2e40)
plt.ylim(ymin=-0.05,ymax=1.3)
plt.xlim(xmin=-29,xmax=-2)
#plt.ylim(ymin=-0.05,ymax=0.5)
#plt.xlim(xmin=-29,xmax=-18)
plt.gca().invert_xaxis()
plt.legend()
plt.tight_layout()

##
##plt.figure(figsize=[12,7])
##plt.subplot(121)
##plt.plot(lL, lsig, 'bo', label='z <= 1.6', zorder=1)
##plt.errorbar(lL, lsig, yerr=lsigerr, fmt='bo', zorder=0, alpha=0.2)
##plt.xscale('log')
##plt.xlabel('K Band Luminosity (W)')
##plt.ylabel(r'$\sigma$')
##plt.ylim(ymin=0,ymax=0.25)
##plt.xlim(xmin=1e37,xmax=1e40)
##plt.legend()
##
##plt.subplot(122)
##plt.plot(hL, hsig, 'go', label='z > 1.6',zorder=1)
##plt.errorbar(hL, hsig, yerr=hsigerr, fmt='go', zorder=0, alpha=0.2)
##plt.xscale('log')
##plt.xlabel('K Band Luminosity (W)')
##plt.ylabel(r'$\sigma$')
##plt.ylim(ymin=0,ymax=0.25)
##plt.xlim(xmin=1e37,xmax=1e40)
##plt.legend()
##plt.tight_layout()
#
##%% Plot the split agains Mstar ###
#highmstar = tbhigh['Mstar_z_p']
#lowmstar = tblow['Mstar_z_p']
#highchanmstar = chanhigh['Mstar_z_p']
#lowchanmstar = chanlow['Mstar_z_p']
#highxmmmstar = xmmhigh['Mstar_z_p']
#lowxmmmstar = xmmlow['Mstar_z_p']
#fullhighmstar = fullhigh['Mstar_z_p']
#fulllowmstar = fulllow['Mstar_z_p']
#highradiomstar = radiohigh['Mstar_z_p']
#lowradiomstar = radiolow['Mstar_z_p']
#
#### combine data ###
#hmstar = np.append(highchanmstar, highxmmmstar)
#lmstar = np.append(lowchanmstar, lowxmmmstar)
#
#plt.figure(figsize=[12,7])
#plt.subplot(121)
#plt.plot(fulllowmstar, fulllowout[:,0],'o', color='tab:grey', label='X-ray non variable', zorder=1)
#plt.errorbar(fulllowmstar, fulllowout[:,0], yerr=fulllowout[:,1], fmt='o', color='tab:grey', zorder=0, alpha=0.2)
#plt.plot(lowmstar, lowout[:,0], 'bo', label='z <= 1.6', zorder=1)
#plt.errorbar(lowmstar, lowout[:,0], yerr=lowout[:,1], fmt='bo', zorder=0, alpha=0.2)
#plt.plot(lmstar, lsig, 'ro', label='X-ray variable', zorder=1)
#plt.errorbar(lmstar, lsig, yerr=lsigerr, fmt='ro', zorder=0, alpha=0.2)
#plt.plot(lowradiomstar, lowradioout[:,0], 'ks', markersize=10,mfc='None', label='Radio Source', zorder=3)
#plt.xscale('log')
#plt.xlabel(r'$M_{stellar}$')
#plt.ylabel(r'$\sigma$')
#plt.ylim(ymin=-0.1,ymax=1.5)
#plt.xlim(xmin=1e4,xmax=2e12)
##plt.ylim(ymin=-0.05,ymax=0.5)
##plt.xlim(xmin=1e35,xmax=2e40)
#plt.legend()
#
#plt.subplot(122)
#plt.plot(fullhighmstar, fullhighout[:,0], 'o', color='tab:grey', label='X-ray non variable', zorder=1)
#plt.errorbar(fullhighmstar, fullhighout[:,0], yerr=fullhighout[:,1], fmt='o', color='tab:grey', zorder=0, alpha=0.2)
#plt.plot(highmstar, highout[:,0], 'go', label='z > 1.6', zorder=1)
#plt.errorbar(highmstar, highout[:,0], yerr=highout[:,1], fmt='go', zorder=0, alpha=0.2)
#plt.plot(hmstar, hsig, 'ro', label='X-ray variable',zorder=2)
#plt.errorbar(hmstar, hsig, yerr=hsigerr, fmt='ro', zorder=0, alpha=0.2)
#plt.plot(highradiomstar, highradioout[:,0], 'ks', markersize=10,mfc='None', label='Radio Source', zorder=3)
#plt.xscale('log')
#plt.xlabel(r'$M_{stellar}$')
#plt.ylabel(r'$\sigma$')
#plt.ylim(ymin=-0.1,ymax=1.5)
#plt.xlim(xmin=1e4,xmax=2e12)
##plt.ylim(ymin=-0.05,ymax=0.5)
##plt.xlim(xmin=1e35,xmax=2e40)
#plt.legend()
#
#end = time.time()
#print(end-start)
#
