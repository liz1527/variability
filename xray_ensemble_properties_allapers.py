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
sigtb1 = Table.read('quad_epoch_sigma_table_extra_clean_no06_067arcsec.fits')
sigtb2 = Table.read('quad_epoch_sigma_table_extra_clean_no06_1arcsec.fits')
sigtb3 = Table.read('quad_epoch_sigma_table_extra_clean_no06_1-5arcsec.fits')
sigtb4 = Table.read('quad_epoch_sigma_table_extra_clean_no06_2arcsec.fits')
sigtb5 = Table.read('quad_epoch_sigma_table_extra_clean_no06.fits')

### Limit to Chandra region for simplicity ###
tbdata = vari_funcs.chandra_only(tbdata)
fullxray = vari_funcs.chandra_only(fullxray)
chandata = vari_funcs.chandra_only(chandata)

def get_ensemble_sig(tbdata, sigtb, binedge, posvar, aper=5):
    # Extract magnitude table and error table
    flux = vari_funcs.flux_stacks(tbdata,aper)
    flux, tbdata = vari_funcs.noneg(flux, tbdata)
    flux, fluxerr, tbdata = vari_funcs.create_quad_error_array(sigtb, tbdata,aper)
    
    ### Normalise ###
    fluxnorm, fluxerrnorm = vari_funcs.normalise_flux_and_errors(flux, fluxerr)
    
    ### Find luminosity distance ###
    z = tbdata['z_spec']#[mask]
    z[z==-1] = tbdata['z_p'][z==-1]
    DL = cosmo.luminosity_distance(z)
    DL = DL.to(u.cm)
    
    ### Calculate the luminosity ###
    xrayF = tbdata['Full_flux']#[chanmask]
    xrayL = xrayF*4*np.pi*(DL.value**2)
    
    ### Create dicts to save data into ###
    enflux = {}
    enfluxerr = {}
    enxrayL = {}
    sig = np.empty(size)
    sigerr = np.empty(size)
    meanxrayL = np.empty(size)
    for m, enmin in enumerate(binedge):
        ### Isolate data needed ###
        mask1 = xrayL >= enmin
        if m != size-1:
            mask2 = xrayL < binedge[m+1]
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


#%% Find luminosity distance for chandra sources ###
chanz = chandata['z_spec']#[chanmask]
chanz[chanz==-1] = chandata['z_p'][chanz==-1]
chanDL = cosmo.luminosity_distance(chanz)
chanDL = chanDL.to(u.cm)

#%% Calculate the luminosity ###
xrayF = chandata['Full_flux']#[chanmask]
xrayL = xrayF*4*np.pi*(chanDL.value**2)

#%% Split into bins according to X-ray luminosity ###
sortedxrayL = np.sort(xrayL)
posvar = np.linspace(0,2,5000)

binsize = 10
binedge = np.empty(int(len(xrayL)/binsize))
size = len(binedge)
for n in range(size):
    ### Define bin edges ###
    enXrayL = sortedxrayL[n*binsize:(n*binsize)+binsize]
    enmin = np.nanmin(enXrayL)
    binedge[n] = enmin


plt.figure(figsize=[8,7])
def plot_sig(tbdata, sigtb, binedge, aper=5, format='o',colour='b', label='Ensemble'):
    posvar = np.linspace(0,2,5000)
    fluxnorm, fluxerrnorm, sig, sigerr, xrayL, meanxrayL = get_ensemble_sig(tbdata, sigtb, binedge,posvar, aper)
    binupper = np.append(binedge[1:],np.max(xrayL))
    xlow = meanxrayL - binedge
    xhigh = binupper - meanxrayL
    plt.errorbar(meanxrayL, sig, xerr=[xlow, xhigh], yerr=sigerr, fmt=format,
                 color=colour, label=label)

### plot for 3" ###
#plot_sig(tbdata, sigtb5, binedge,aper=5,format = 'o',
#                 colour='C0', label='3"')#label='Ensemble Variable Sources')
plot_sig(chandata, sigtb5, binedge,aper=5,format = 'o',
                 colour='C0', label='3"')#, label='Ensemble Variable X-ray Sources')
#plot_sig(fullxray, sigtb5, binedge,aper=5,format = 'o',
#                 colour='k', label='Ensemble Non Variable X-ray Sources')

### plot for 2" ###
#plot_sig(tbdata, sigtb4, binedge,aper=4,format = '+',
#                 colour='C1',label='2"')# label='Ensemble Variable Sources')
#plot_sig(chandata, sigtb4, binedge,aper=4,format = '+',
#                 colour='C1',label='2"')#, label='Ensemble Variable X-ray Sources')
#plot_sig(fullxray, sigtb4, binedge,aper=4, format = '+',
#                 colour='k')#, label='Ensemble Non Variable X-ray Sources')

### plot for 1.5" ###
#plot_sig(tbdata, sigtb3, binedge,aper=3,format = 'd',
#                 colour='C2', label='1.5"')#label='Ensemble Variable Sources')
#plot_sig(chandata, sigtb3, binedge,aper=3,format = 'd',
#                 colour='C2', label='1.5"')#, label='Ensemble Variable X-ray Sources')
#plot_sig(fullxray, sigtb3, binedge,aper=3, format = 'd',
#                 colour='k')#, label='Ensemble Non Variable X-ray Sources')

### plot for 1" ###
#plot_sig(tbdata, sigtb2, binedge,aper=2,format = 's',
#                 colour='C3', label='1"')#label='Ensemble Variable Sources')
#plot_sig(chandata, sigtb2, binedge,aper=2,format = 's',
#                 colour='C3', label='1"')#, label='Ensemble Variable X-ray Sources')
#plot_sig(fullxray, sigtb2, binedge,aper=2, format = 's',
#                 colour='k')#, label='Ensemble Non Variable X-ray Sources')

### plot for 0.67" ###
#plot_sig(tbdata, sigtb1, binedge,aper=1,format = 'x',
#                 colour='C4', label='0.67"')#label='Ensemble Variable Sources')
plot_sig(chandata, sigtb1, binedge,aper=1,format = 'x',
                 colour='C4', label='0.67"')#, label='Ensemble Variable X-ray Sources')
#plot_sig(fullxray, sigtb1, binedge,aper=1, format = 'x',
#                 colour='k')#, label='Ensemble Non Variable X-ray Sources')

plt.xscale('log')
plt.xlabel('X-ray Luminosity (ergs/s)')
plt.ylabel(r'$\sigma$')
#plt.ylim(ymin=-0.1,ymax=0.5)
plt.xlim(xmin=1e43)
plt.legend()





















end = time.time()
print(end-start)

