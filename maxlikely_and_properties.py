#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  1 13:30:42 2018

Code to compute maximum likelihood variability

@author: ppxee
"""
import time
start = time.time()
#print(start)
print('Code was started at '+time.ctime())

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

#%% Open the fits files and get data ###
tbdata = fits.open('variable_tables/no06_variables_chi30_DR11data_restframe.fits')[1].data
fullxray = fits.open('mag_flux_tables/chanDR11data_restframe_mag_flux_table_best_extra_clean_no06.fits')[1].data
chandata = fits.open('variable_tables/no06_variables_chi30_chandata_DR11data_restframe.fits')[1].data
sigtb = Table.read('sigma_tables/quad_epoch_sigma_table_extra_clean_no06.fits')

### Limit to Chandra region for simplicity ###
tbdata = vari_funcs.chandra_only(tbdata)
fullxray = vari_funcs.chandra_only(fullxray)
chandata = vari_funcs.chandra_only(chandata)

# Extract magnitude table and error table
flux = vari_funcs.flux5_stacks(tbdata)
flux, tbdata = vari_funcs.noneg(flux, tbdata)
flux, fluxerr, tbdata = vari_funcs.create_quad_error_array(sigtb, tbdata, aper=5)
chanflux = vari_funcs.flux5_stacks(chandata)
chanflux, chandata = vari_funcs.noneg(chanflux, chandata)
chanflux, chanerr, chandata = vari_funcs.create_quad_error_array(sigtb, chandata, aper=5)
fullflux = vari_funcs.flux5_stacks(fullxray)
fullflux, fulldata = vari_funcs.noneg(fullflux, fullxray)
fullflux, fullerr, fullxray = vari_funcs.create_quad_error_array(sigtb, fullxray, aper=5)

### Normalise ###
fluxnorm, fluxerrnorm = vari_funcs.normalise_flux_and_errors(flux, fluxerr)
chanfluxnorm, chanerrnorm = vari_funcs.normalise_flux_and_errors(chanflux, chanerr)
fullfluxnorm, fullerrnorm = vari_funcs.normalise_flux_and_errors(fullflux, fullerr)
#%% All points
posvar = np.linspace(0,2,5000)
#start = time.time()

numobs = np.shape(fluxnorm)[0]
meanflux = np.nanmean(fluxnorm, axis=1)
out = np.array([vari_funcs.maximum_likelihood(fluxnorm[n,:], fluxerrnorm[n,:], meanflux[n], posvar) for n in range(numobs)])

numobs = np.shape(chanfluxnorm)[0]
meanchan = np.nanmean(chanfluxnorm, axis=1)
chanout = np.array([vari_funcs.maximum_likelihood(chanfluxnorm[n,:], chanerrnorm[n,:], meanchan[n], posvar) for n in range(numobs)])

numobs = np.shape(fullfluxnorm)[0]
meanfull = np.nanmean(fullfluxnorm, axis=1)
fullout = np.array([vari_funcs.maximum_likelihood(fullfluxnorm[n,:], fullerrnorm[n,:], meanfull[n], posvar) for n in range(numobs)])

#%% Plots
#out2 = out
#out2[out2[:,0]==0] = np.nan
meanflux = np.nanmean(flux, axis=1)
meanchan = np.nanmean(chanflux, axis=1)
meanfull = np.nanmean(fullflux, axis=1)

#plt.figure()
#plt.errorbar(meanflux, out[:,0], out[:,1], fmt='x', zorder=0)
#plt.ylabel('Maximum Likelihood')
#plt.xlabel('Mean Flux')
#plt.xscale('log')
##plt.figure()
#plt.scatter(meanchan, chanout[:,0], s=80, c='none', marker='o', zorder=3,edgecolors='r')
#plt.scatter(meanflux, out[:,0],c='b', marker='x', zorder=2)
#plt.ylabel('Maximum Likelihood')
#plt.xlabel('Mean Flux')

### Plot against X-ray Luminosity ###
### Define cosmology ###
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

### Find luminosity distance both for all and just for chandra sources ###
#mask = tbdata['z_spec']>=0
z = tbdata['z_spec']#[mask]
z[z==-1] = tbdata['z_p'][z==-1]
DL = cosmo.luminosity_distance(z)
DL = DL.to(u.cm)
#chanmask = chandata['z_spec']>=0
chanz = chandata['z_spec']#[chanmask]
chanz[chanz==-1] = chandata['z_p'][chanz==-1]
chanDL = cosmo.luminosity_distance(chanz)
chanDL = chanDL.to(u.cm)
#fullmask = fullxray['z_spec']>=0
fullz = fullxray['z_spec']#[fullmask]
fullz[fullz==-1] = fullxray['z_p'][fullz==-1]
fullDL = cosmo.luminosity_distance(fullz)
fullDL = fullDL.to(u.cm)

### Calculate the luminosity ###
xrayF = chandata['Full_flux']#[chanmask]
xrayL = xrayF*4*np.pi*chanDL.value**2
fullxrayF = fullxray['Full_flux']#[fullmask]
fullxrayL = fullxrayF*4*np.pi*fullDL.value**2

### Plot against max likelihood ###
plt.figure(figsize=[7,7])
plt.errorbar(fullxrayL, fullout[:,0],yerr=fullout[:,1],fmt='o',
             color='tab:grey', label='Non variable X-ray source',zorder=1)
plt.errorbar(xrayL, chanout[:,0],yerr=chanout[:,1],fmt='o',
             color='r',label='Variable X-ray Source',zorder=2)
plt.xlabel('X-ray Luminosity (ergs/s)')
plt.ylabel(r'$\sigma$')
plt.xlim(xmin=1e40, xmax=8e45)
plt.ylim(ymin=-0.3, ymax=0.7)
plt.xscale('log')
plt.legend()
plt.title('Amplitude of variability against x-ray luminosity')
plt.tight_layout()

### Plot redshift against max likelihood ###
plt.figure(figsize=[7,7])
plt.errorbar(fullz, fullout[:,0],yerr=fullout[:,1], fmt='o', 
             color='tab:gray',label='Non Variable X-ray Source',zorder=1)
plt.errorbar(z, out[:,0],yerr=out[:,1], fmt='o',color='b',
             label='Variable Non X-ray Source',zorder=2)
plt.errorbar(chanz, chanout[:,0],yerr=chanout[:,1], fmt='o', 
             color='r',label='Variable X-ray Source',zorder=3)
plt.xlabel('z')
plt.xlim(xmin=-0.1, xmax=4.2)
plt.ylim(ymin=-0.3, ymax=1.75)
plt.ylabel(r'$\sigma$')
plt.title('Amplitude of variability against redshift')
plt.legend()
plt.tight_layout()

#plt.figure(figsize=[7,7])
#plt.hist([out[:,0], chanout[:,0], fullout[:,0]], bins=50, 
#         color = ['b','r','tab:grey'], normed=True)

#%% Plot maximum likelihood against stellar mass ###
mstar = tbdata['Mstar_z_p']
chanmstar = chandata['Mstar_z_p']
fullmstar = fullxray['Mstar_z_p']

plt.figure(figsize=[8,7])
plt.errorbar(fullmstar, fullout[:,0],yerr=fullout[:,1], fmt='o', 
             color='tab:gray',label='Non Variable X-ray Source',zorder=1)
plt.errorbar(mstar, out[:,0],yerr=out[:,1], fmt='o',color='b',
             label='Variable Non X-ray Source',zorder=2)
plt.errorbar(chanmstar, chanout[:,0],yerr=chanout[:,1], fmt='o', 
             color='r',label='Variable X-ray Source',zorder=3)
plt.xscale('log')
plt.xlabel(r'$M_{stellar}$')
plt.ylabel(r'$\sigma$')
plt.xlim(xmin=1e8, xmax=6e11)
plt.ylim(ymin=-0.3, ymax=1.75)
plt.legend()
plt.tight_layout()

#%% Plot against eddington ratio (~Lx/Mstar)
chanLedd = 1.26e31 * (0.1*chanmstar)
chaneddrat = np.log(10*xrayL) - np.log(chanLedd)
fullLedd = 1.26e31 * (0.1*fullmstar)
fulleddrat = np.log(10*fullxrayL) - np.log(fullLedd)

plt.figure(figsize=[8,7])
plt.errorbar(fulleddrat, fullout[:,0],yerr=fullout[:,1], fmt='o', 
             color='tab:gray',label='Non Variable X-ray Source',zorder=1)
plt.errorbar(chaneddrat, chanout[:,0],yerr=chanout[:,1], fmt='o', 
             color='r',label='Variable X-ray Source',zorder=3)
plt.xlabel(r'$log(\lambda_{edd}) = log(10 L_{X}) - log(L_{edd})$')
plt.ylabel(r'$\sigma$')
#plt.xlim(xmin=73)
#plt.ylim(ymin=-0.3, ymax=1.75)
plt.legend()
plt.tight_layout()

##%% Plot against LD
#plt.figure(figsize=[8,7])
#plt.errorbar(fullDL.value, fullout[:,0],yerr=fullout[:,1], fmt='o', 
#             color='tab:gray',label='Non Variable X-ray Source',zorder=1)
#plt.errorbar(DL.value, out[:,0],yerr=out[:,1], fmt='o',color='b',
#             label='Variable Non X-ray Source',zorder=2)
#plt.errorbar(chanDL.value, chanout[:,0],yerr=chanout[:,1], fmt='o', 
#             color='r',label='Variable X-ray Source',zorder=3)
#plt.xlabel(r'$D_{Luminosity}$')
#plt.ylabel(r'$\sigma$')
#plt.legend()
#plt.tight_layout()

end = time.time()
print(end-start)