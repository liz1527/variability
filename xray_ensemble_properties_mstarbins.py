#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 29 12:13:52 2018

Code to look at ensemble variability of objects in variety of aperture sizes

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
#sigtb1 = Table.read('quad_epoch_sigma_table_extra_clean_no06_067arcsec.fits')
#sigtb2 = Table.read('quad_epoch_sigma_table_extra_clean_no06_1arcsec.fits')
#sigtb3 = Table.read('quad_epoch_sigma_table_extra_clean_no06_1-5arcsec.fits')
#sigtb4 = Table.read('quad_epoch_sigma_table_extra_clean_no06_2arcsec.fits')
sigtb5 = Table.read('quad_epoch_sigma_table_extra_clean_no06.fits')

### Limit to Chandra region for simplicity ###
tbdata = vari_funcs.chandra_only(tbdata)
fullxray = vari_funcs.chandra_only(fullxray)
chandata = vari_funcs.chandra_only(chandata)

def get_and_normalise_flux(tbdata, sigtb, aper=5):
    # Extract magnitude table and error table
    flux = vari_funcs.flux_stacks(tbdata,aper)
    flux, tbdata = vari_funcs.noneg(flux, tbdata)
    flux, fluxerr, tbdata = vari_funcs.create_quad_error_array(sigtb, tbdata,aper)
    
    ### Normalise ###
    fluxnorm, fluxerrnorm = vari_funcs.normalise_flux_and_errors(flux, fluxerr)
    return fluxnorm, fluxerrnorm, tbdata
    
def get_ensemble_sig(tbdata, sigtb, binedge, posvar, aper=5):
    # Extract magnitude table and error table
    fluxnorm, fluxerrnorm, tbdata = get_and_normalise_flux(tbdata, sigtb, aper=5)
    
    ### Find luminosity distance ###
    z = tbdata['z_spec']#[mask]
    z[z==-1] = tbdata['z_p'][z==-1]
    DL = cosmo.luminosity_distance(z)
    DL = DL.to(u.cm)
    
    ### Calculate the luminosity ###
    xrayF = tbdata['Full_flux']#[chanmask]
    xrayL = xrayF*4*np.pi*(DL.value**2)
    
    ### Get stellar mass ###
    Mstar = tbdata['Mstar_z_p']
    
    ### Create dicts to save data into ###
    enflux = {}
    enfluxerr = {}
    enxrayL = {}
    sig = np.empty(size)
    sigerr = np.empty(size)
    meanxrayL = np.empty(size)
    for m, enmin in enumerate(binedge):
        ### Isolate data needed ###
        mask1 = Mstar >= enmin
        if m != size-1:
            mask2 = Mstar < binedge[m+1]
        else:
            mask2 = np.ones(len(mask1))
        
        enmask = mask1*mask2.astype(bool)
    
        enflux[m] = fluxnorm[enmask]
        enfluxerr[m] = fluxerrnorm[enmask]
        enxrayL[m] = xrayL[enmask]
        
        ### Combine into one flux curve per bin ###
        enfluxcurve = np.ravel(enflux[m])
        enfluxcurveerr = np.ravel(enfluxerr[m])
        
        
        ### Find max likelihood sig of curve ###
        [sig[m],sigerr[m]] = vari_funcs.maximum_likelihood(enfluxcurve, enfluxcurveerr, 1, posvar)
        
        ### find mean z ###
        meanxrayL[m] = np.nanmean(enxrayL[m])
        
    return fluxnorm, fluxerrnorm, sig, sigerr, xrayL, meanxrayL

def plot_sig(tbdata, sigtb, binedge, aper=5, format='o',colour='b', label='Ensemble'):
    posvar = np.linspace(0,2,5000)
    fluxnorm, fluxerrnorm, sig, sigerr, xrayL, meanxrayL = get_ensemble_sig(tbdata, sigtb, binedge,posvar, aper)

    for k, enmin in enumerate(binedge):
        if enmin == binedge[-1]:
            plt.errorbar(meanxrayL[k], sig[k], yerr=sigerr[k], fmt=format, 
                         label='{0:1.2e}'.format(enmin)+r'$<M_{star}$', markersize=10,
                         color='C'+str(k))
        else:                
            plt.errorbar(meanxrayL[k], sig[k], yerr=sigerr[k], fmt=format, 
                     label='{0:1.2e}'.format(enmin)+r'$<M_{star}<$'+'{0:1.2e}'.format(binedge[k+1]), 
                     markersize=10,color='C'+str(k))

#%% Find stellar mass for chandra sources and create bins ###
Mstar = chandata['Mstar_z_p']
sortedMstar = np.sort(Mstar)
posvar = np.linspace(0,2,5000)

binsize = 8
binedge = np.empty(int(len(Mstar)/binsize))
size = len(binedge)
plt.figure(figsize=[8,7])
for n in range(size):
    ### Define bin edges ###
    enXrayL = sortedMstar[n*binsize:(n*binsize)+binsize]
    enmin = np.nanmin(enXrayL)
    binedge[n] = enmin


for m, enmin in enumerate(binedge):
    ### Isolate data needed ###
    mask1 = Mstar >= enmin
    if m != size-1:
        mask2 = Mstar < binedge[m+1]
    else:
        mask2 = np.ones(len(mask1))
    
    enmask = mask1*mask2.astype(bool)
    tempdata = chandata[enmask]
    print(len(tempdata))
    #%% Plot the single sources colour coded by M_star bin ###
    chanfluxnorm, chanerrnorm, tempdata = get_and_normalise_flux(tempdata, sigtb5, aper=5)
    numobs = np.shape(chanfluxnorm)[0]
    meanchan = np.nanmean(chanfluxnorm, axis=1)
    chanout = np.array([vari_funcs.maximum_likelihood(chanfluxnorm[n,:], chanerrnorm[n,:], meanchan[n], posvar) for n in range(numobs)])
    
    ### Find luminosity distance for chandra sources ###
    chanz = tempdata['z_spec']#[chanmask]
    chanz[chanz==-1] = tempdata['z_p'][chanz==-1]
    chanDL = cosmo.luminosity_distance(chanz)
    chanDL = chanDL.to(u.cm)

    ### Calculate the luminosity ###
    xrayF = tempdata['Full_flux']#[chanmask]
    xrayL = xrayF*4*np.pi*(chanDL.value**2)
    
#    if enmin == binedge[-1]:
#        plt.errorbar(xrayL, chanout[:,0],yerr=chanout[:,1],fmt='o',
#                     label=str(enmin)+r'$<M_{star}$',zorder=0,alpha=0.5)
#    else:
#        plt.errorbar(xrayL, chanout[:,0],yerr=chanout[:,1],fmt='o',
#                     label=str(enmin)+r'$<M_{star}<$'+str(binedge[m+1]),
#                     zorder=0,alpha=0.5)
#    
    if enmin == binedge[-1]:
        plt.plot(xrayL, chanout[:,0],'o',
#                     label=str(enmin)+r'$<M_{star}$',
                     zorder=0,alpha=0.5)
    else:
        plt.plot(xrayL, chanout[:,0],'o',
#                     label=str(enmin)+r'$<M_{star}<$'+str(binedge[m+1]),
                     zorder=0,alpha=0.5)

### plot for 3" ###
#plot_sig(tbdata, sigtb5, binedge,aper=5,format = 'o',
#                 colour='C0', label='3"')#label='Ensemble Variable Sources')
plot_sig(chandata, sigtb5, binedge,aper=5,format = 's',
                 colour='C0', label='3"')#, label='Ensemble Variable X-ray Sources')
#plot_sig(fullxray, sigtb5, binedge,aper=5,format = 'o',
#                 colour='k', label='Ensemble Non Variable X-ray Sources')


plt.xscale('log')
plt.xlabel('X-ray Luminosity (ergs/s)')
plt.ylabel(r'$\sigma$')
#plt.ylim(ymin=-0.1,ymax=0.5)
#plt.xlim(xmin=1e43)
plt.legend()





















end = time.time()
print(end-start)

