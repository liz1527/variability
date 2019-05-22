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
        
def z_split(tbdata):
    z = tbdata['z_spec']#[mask]
    z[z==-1] = tbdata['z_p'][z==-1]
    tbhigh = tbdata[z > 2]
    tblow = tbdata[z <= 2]
    return tbhigh, tblow

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

#%% Split by z ###

tbhigh, tblow = z_split(tbdata)
xmmhigh, xmmlow = z_split(xmmdata)
chanhigh, chanlow = z_split(chandata)
fullhigh, fulllow = z_split(fullxray)

#%% get_luminosity_and_flux ###

tbhigh, highL, highfluxnorm, highfluxerrnorm= get_luminosity_and_flux(tbhigh)
tblow, lowL, lowfluxnorm, lowfluxerrnorm= get_luminosity_and_flux(tblow)
chanhigh, highchanL, highchanfluxnorm, highchanerrnorm= get_luminosity_and_flux(chanhigh)
chanlow, lowchanL, lowchanfluxnorm, lowchanerrnorm= get_luminosity_and_flux(chanlow)
xmmhigh, highxmmL, highxmmfluxnorm, highxmmerrnorm = get_luminosity_and_flux(xmmhigh, xmm=True)
xmmlow, lowxmmL, lowxmmfluxnorm, lowxmmerrnorm = get_luminosity_and_flux(xmmlow, xmm=True)
fullhigh, highfullL, highfullfluxnorm, highfullerrnorm= get_luminosity_and_flux(fullhigh)
fulllow, lowfullL, lowfullfluxnorm, lowfullerrnorm= get_luminosity_and_flux(fulllow)

#%% Split into bins according to nir luminosity ###

### Join xmm and chan arrays ###
hxrayL = np.append(highxmmL, highchanL)
lxrayL = np.append(lowxmmL, lowchanL)

### Sort arrays so can find bin edges ###
sortedhighL = np.sort(highL)
sortedlowL = np.sort(lowL)
sortedhighxrayL = np.sort(hxrayL)
sortedlowxrayL = np.sort(lxrayL)
sortedhighfullL = np.sort(highfullL)
sortedlowfullL = np.sort(lowfullL)

### define bin size and find edge arrays ###
binsize = 8
highbinedge = get_binedge(binsize, sortedhighL)
lowbinedge = get_binedge(binsize, sortedlowL)
highxraybinedge = get_binedge(binsize, sortedhighxrayL)
lowxraybinedge = get_binedge(binsize, sortedlowxrayL)
highfullbinedge = get_binedge(binsize, sortedhighfullL)
lowfullbinedge = get_binedge(binsize, sortedlowfullL)

#%% create ensemble arrays ###

highenflux, highenfluxerr, highenL = make_ensemble(tbhigh, highL, highbinedge,
                                                   highfluxnorm, highfluxerrnorm)
lowenflux, lowenfluxerr, lowenL = make_ensemble(tblow, lowL, lowbinedge, 
                                                lowfluxnorm, lowfluxerrnorm)
highchanenflux, highchanenfluxerr, highchanenL = make_ensemble(chanhigh, highchanL, 
                                                               highxraybinedge,
                                                               highchanfluxnorm, 
                                                               highchanerrnorm)
lowchanenflux, lowchanenfluxerr, lowchanenL = make_ensemble(chanlow, lowchanL, 
                                                            lowxraybinedge, 
                                                            lowchanfluxnorm, 
                                                            lowchanerrnorm)
highxmmenflux, highxmmenfluxerr, highxmmenL = make_ensemble(xmmhigh, highxmmL, 
                                                            highxraybinedge,
                                                            highxmmfluxnorm,
                                                            highxmmerrnorm )
lowxmmenflux, lowxmmenfluxerr, lowxmmenL = make_ensemble(xmmlow, lowxmmL, 
                                                         lowxraybinedge, 
                                                         lowxmmfluxnorm, 
                                                         lowxmmerrnorm)
highfullenflux, highfullenfluxerr, highfullenL = make_ensemble(fullhigh, 
                                                               highfullL, 
                                                               highfullbinedge,
                                                               highfullfluxnorm, 
                                                               highfullerrnorm)
lowfullenflux, lowfullenfluxerr, lowfullenL = make_ensemble(fulllow, lowfullL, 
                                                            lowfullbinedge, 
                                                            lowfullfluxnorm, 
                                                            lowfullerrnorm)


#%% Combine chan and xmm ensembles ###

henflux, henfluxerr, henL = combine_Xray(highxraybinedge, highxmmenflux,
                                      highxmmenfluxerr, highxmmenL, 
                                      highchanenflux, highchanenfluxerr, 
                                      highchanenL)
lenflux, lenfluxerr, lenL = combine_Xray(lowxraybinedge, lowxmmenflux,
                                      lowxmmenfluxerr, lowxmmenL, 
                                      lowchanenflux, lowchanenfluxerr, 
                                      lowchanenL)

#%% Get ensemble maximum likelihood ###
highsig, highsigerr, highmeanL = get_ensemble_sig(highenflux, highenfluxerr, highenL, posvar)
lowsig, lowsigerr, lowmeanL = get_ensemble_sig(lowenflux, lowenfluxerr, lowenL, posvar)
highxraysig, highxraysigerr, highxraymeanL = get_ensemble_sig(henflux, henfluxerr, henL, posvar)
lowxraysig, lowxraysigerr, lowxraymeanL = get_ensemble_sig(lenflux, lenfluxerr, lenL, posvar)
highfullsig, highfullsigerr, highfullmeanL = get_ensemble_sig(highfullenflux, highfullenfluxerr, highfullenL, posvar)
lowfullsig, lowfullsigerr, lowfullmeanL = get_ensemble_sig(lowfullenflux, lowfullenfluxerr, lowfullenL, posvar)

#%% Plot results ###
def find_err_values(binedge, sortedL, meanL):
    binupper = np.append(binedge[1:],np.max(sortedL))
    xlow = meanL-binedge
    xhigh = binupper - meanL
    return xlow, xhigh
    
highxerr = find_err_values(highbinedge, sortedhighL, highmeanL)
lowxerr = find_err_values(lowbinedge, sortedlowL, lowmeanL)
highxrayxerr = find_err_values(highxraybinedge, sortedhighxrayL, highxraymeanL)
lowxrayxerr = find_err_values(lowxraybinedge, sortedlowxrayL, lowxraymeanL)
highfullxerr = find_err_values(highfullbinedge, sortedhighfullL, highfullmeanL)
lowfullxerr = find_err_values(lowfullbinedge, sortedlowfullL, lowfullmeanL)


plt.figure(figsize=[12,7])
#plt.subplot(121)
plt.errorbar(lowmeanL, lowsig, xerr=lowxerr, yerr=lowsigerr, fmt='o',
             color='b', label='z <= 1.6')
#plt.errorbar(lowxraymeanL, lowxraysig, xerr=lowxrayxerr, yerr=lowxraysigerr, fmt='o',
#             color='r', label='Variable X-ray')
#plt.errorbar(lowfullmeanL, lowfullsig, xerr=lowfullxerr, 
#             yerr=lowfullsigerr, fmt='o',color='tab:grey', label='Non Variable X-ray')
plt.xlabel('K Band Absolute Magnitude')
plt.ylabel(r'$\sigma$')
plt.ylim(ymin=-0.05,ymax=1.5)
plt.xlim(xmin=-29,xmax=-2)
#plt.ylim(ymin=-0.05,ymax=0.5)
#plt.xlim(xmin=-29,xmax=-18)
plt.gca().invert_xaxis()
plt.legend()

#plt.subplot(122)
plt.errorbar(highmeanL, highsig, xerr=highxerr, yerr=highsigerr, fmt='o',
             color='g', label='z > 1.6')
#plt.errorbar(highxraymeanL, highxraysig, xerr=highxrayxerr, yerr=highxraysigerr, fmt='o',
#             color='r', label='Variable X-ray')
#plt.errorbar(highfullmeanL, highfullsig, xerr=highfullxerr, 
#             yerr=highfullsigerr, fmt='o',color='tab:grey', label='Non Variable X-ray')
plt.xlabel('K Band Absolute Magnitude')
plt.ylabel(r'$\sigma$')
plt.ylim(ymin=-0.05,ymax=1.3)
plt.xlim(xmin=-29,xmax=-2)
#plt.ylim(ymin=-0.05,ymax=0.5)
#plt.xlim(xmin=-29,xmax=-18)
plt.gca().invert_xaxis()
plt.legend()
plt.tight_layout()





















end = time.time()
print(end-start)

