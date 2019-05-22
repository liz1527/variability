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

#%% Split by z ###
z = xmmdata['z_spec']#[mask]
z[z==-1] = xmmdata['z_p'][z==-1]
xmmhigh = xmmdata[z > 1.6]
xmmlow = xmmdata[z <= 1.6]
chanz = chandata['z_spec']#[mask]
chanz[chanz==-1] = chandata['z_p'][chanz==-1]
chanhigh = chandata[chanz > 1.6]
chanlow = chandata[chanz <= 1.6]
fullz = fullxray['z_spec']#[mask]
fullz[fullz==-1] = fullxray['z_p'][fullz==-1]
fullhigh = fullxray[fullz > 1.6]
fulllow = fullxray[fullz <= 1.6]

#%%
chanhigh, highchanxrayL, highchanfluxnorm, highchanerrnorm= get_luminosity_and_flux(chanhigh)
chanlow, lowchanxrayL, lowchanfluxnorm, lowchanerrnorm= get_luminosity_and_flux(chanlow)
fullhigh, fullhighxrayL, fullhighfluxnorm, fullhigherrnorm= get_luminosity_and_flux(fullhigh)
fulllow, fulllowxrayL, fulllowfluxnorm, fulllowerrnorm= get_luminosity_and_flux(fulllow)
xmmhigh, highxmmxrayL, highxmmfluxnorm, highxmmerrnorm = get_luminosity_and_flux(xmmhigh, xmm=True)
xmmlow, lowxmmxrayL, lowxmmfluxnorm, lowxmmerrnorm = get_luminosity_and_flux(xmmlow, xmm=True)

#%% get sig values
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

#%% Combine into single data set ###
hxrayL = np.append(highchanxrayL, highxmmxrayL)
hsig = np.append(highchanout[:,0],highxmmout[:,0])
hsigerr = np.append(highchanout[:,1],highxmmout[:,1])
lxrayL = np.append(lowchanxrayL, lowxmmxrayL)
lsig = np.append(lowchanout[:,0], lowxmmout[:,0])
lsigerr = np.append(lowchanout[:,1], lowxmmout[:,1])

#%% Plot results ###
plt.figure(figsize=[12,7])
plt.subplot(121)
plt.plot(lxrayL, lsig, 'bo', label='z <= 1.6', zorder=1)
plt.errorbar(lxrayL, lsig, yerr=lsigerr, fmt='bo', zorder=0, alpha=0.2)
plt.plot(fulllowxrayL, fulllowout[:,0],'o', color='tab:grey', label='X-ray non variable', zorder=0, alpha=0.7)
plt.errorbar(fulllowxrayL, fulllowout[:,0], yerr=fulllowout[:,1], fmt='o', color='tab:grey', zorder=0, alpha=0.2)
plt.xscale('log')
plt.xlabel('X-ray Luminosity (ergs/s)')
plt.ylabel(r'$\sigma$')
plt.ylim(ymin=-0.05,ymax=0.4)
plt.xlim(xmin=7e40,xmax=7e45)
plt.legend()

plt.subplot(122)
plt.errorbar(hxrayL, hsig, yerr=hsigerr, fmt='go',zorder=0, alpha=0.2)
plt.plot(hxrayL, hsig, 'go', label='z > 1.6', zorder=1)
plt.errorbar(hxrayL, hsig, yerr=hsigerr, fmt='go', zorder=0, alpha=0.2)
plt.plot(fullhighxrayL, fullhighout[:,0], 'o', color='tab:grey', label='X-ray non variable', zorder=0, alpha=0.7)
plt.errorbar(fullhighxrayL, fullhighout[:,0], yerr=fullhighout[:,1], fmt='o', color='tab:grey', zorder=0, alpha=0.2)
plt.xscale('log')
plt.xlabel('X-ray Luminosity (ergs/s)')
plt.ylabel(r'$\sigma$')
plt.ylim(ymin=-0.05,ymax=0.4)
plt.xlim(xmin=7e40,xmax=7e45)
plt.legend()
plt.tight_layout()



#%% Plot against eddington ratio (~Lx/Mstar)

### Get stellar mass ###
highchanmstar = chanhigh['Mstar_z_p']
lowchanmstar = chanlow['Mstar_z_p']
highxmmmstar = xmmhigh['Mstar_z_p']
lowxmmmstar = xmmlow['Mstar_z_p']
fullhighmstar = fullhigh['Mstar_z_p']
fulllowmstar = fulllow['Mstar_z_p']

### combine data ###
hmstar = np.append(highchanmstar, highxmmmstar)
lmstar = np.append(lowchanmstar, lowxmmmstar)

hLedd = 1.26e31 * (1e-3*hmstar)
heddrat = np.log(10*hxrayL*1e-7) - np.log(hLedd)
lLedd = 1.26e31 * (1e-3*lmstar)
leddrat = np.log(10*lxrayL*1e-7) - np.log(lLedd)
fullhighLedd = 1.26e31 * (1e-3*fullhighmstar)
fullhigheddrat = np.log(10*fullhighxrayL*1e-7) - np.log(fullhighLedd)
fulllowLedd = 1.26e31 * (1e-3*fulllowmstar)
fullloweddrat = np.log(10*fulllowxrayL*1e-7) - np.log(fulllowLedd)

### Plot results ###
plt.figure(figsize=[12,7])
plt.subplot(121)
plt.plot(leddrat, lsig, 'bo', label='z <= 1.6', zorder=1)
plt.errorbar(leddrat, lsig, yerr=lsigerr, fmt='bo', zorder=0, alpha=0.2)
plt.plot(fullloweddrat, fulllowout[:,0],'o', color='tab:grey', label='X-ray non variable', zorder=0, alpha=0.7)
plt.errorbar(fullloweddrat, fulllowout[:,0], yerr=fulllowout[:,1], fmt='o', color='tab:grey', zorder=0, alpha=0.2)
#plt.xscale('log')
plt.xlabel(r'$log(\lambda_{edd}) = log(10 L_{X}) - log(L_{edd})$')
plt.ylabel(r'$\sigma$')
plt.ylim(ymin=-0.05,ymax=0.4)
plt.xlim(xmin=-11,xmax=2)
plt.legend()

plt.subplot(122)
plt.errorbar(heddrat, hsig, yerr=hsigerr, fmt='go',zorder=0, alpha=0.2)
plt.plot(heddrat, hsig, 'go', label='z > 1.6', zorder=1)
plt.errorbar(heddrat, hsig, yerr=hsigerr, fmt='go', zorder=0, alpha=0.2)
plt.plot(fullhigheddrat, fullhighout[:,0], 'o', color='tab:grey', label='X-ray non variable', zorder=0, alpha=0.7)
plt.errorbar(fullhigheddrat, fullhighout[:,0], yerr=fullhighout[:,1], fmt='o', color='tab:grey', zorder=0, alpha=0.2)
#plt.xscale('log')
plt.xlabel(r'$log(\lambda_{edd}) = log(10 L_{X}) - log(L_{edd})$')
plt.ylabel(r'$\sigma$')
plt.ylim(ymin=-0.05,ymax=0.4)
plt.xlim(xmin=-11,xmax=2)
plt.legend()
plt.tight_layout()

end = time.time()
print(end-start)

