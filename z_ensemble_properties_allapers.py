#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 29 12:13:52 2018

Code to look at ensemble variability of objects

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

### Get bins from chandata ###
chanz = chandata['z_spec']#[chanmask]
chanz[chanz==-1] = chandata['z_p'][chanz==-1]
sortedchanz = np.sort(chanz)
binsize = 10
binedge = np.empty(int(len(chanz)/binsize))
size = len(binedge)
for n in range(size):
    ### Define bin edges ###
    enz = sortedchanz[n*binsize:(n*binsize)+binsize]
    enmin = np.nanmin(enz)
    binedge[n] = enmin

def get_ensemble_sig(tbdata, sigtb, binedge, posvar, aper=5):
    # Extract magnitude table and error table
    flux = vari_funcs.flux_stacks(tbdata,aper)
    flux, tbdata = vari_funcs.noneg(flux, tbdata)
    flux, fluxerr, tbdata = vari_funcs.create_quad_error_array(sigtb, tbdata,aper)
    
    ### Normalise ###
    fluxnorm, fluxerrnorm = vari_funcs.normalise_flux_and_errors(flux, fluxerr)
    
    ### Find z ###
    z = tbdata['z_spec']#[mask]
    z[z==-1] = tbdata['z_p'][z==-1]
    
    ### Create dicts to save data into ###
    enflux = {}
    enfluxerr = {}
    enz = {}
    sig = np.empty(size)
    sigerr = np.empty(size)
    meanz = np.empty(size)
    for m, enmin in enumerate(binedge):
        ### Isolate data needed ###
        mask1 = z >= enmin
        if m != size-1:
            mask2 = z < binedge[m+1]
        else:
            mask2 = z < 4.5#np.ones(len(mask1))
        
        enmask = mask1*mask2.astype(bool)
    
        enflux[m] = fluxnorm[enmask]
        enfluxerr[m] = fluxerrnorm[enmask]
        enz[m] = z[enmask]
        
        ### Combine into one flux curve per bin ###
        enfluxcurve = np.ravel(enflux[m])
        enfluxcurveerr = np.ravel(enfluxerr[m])
        
        
        ### Find max likelihood sig of curve ###
        [sig[m],sigerr[m]] = vari_funcs.maximum_likelihood(enfluxcurve, enfluxcurveerr, 1, posvar)
        
        ### find mean z ###
        meanz[m] = np.nanmean(enz[m])
        
    return fluxnorm, fluxerrnorm, sig, sigerr, z, meanz
 

#posvar = np.linspace(0,2,5000)
#fluxnorm, fluxerrnorm, sig, sigerr, z, meanz = get_ensemble_sig(tbdata, sigtb5, binedge,posvar, aper=5)
#chanfluxnorm, chanerrnorm, csig, csigerr, chanz, cmeanz = get_ensemble_sig(chandata, sigtb5, binedge,posvar, aper=5)
#fullfluxnorm, fullerrnorm, fsig, fsigerr, fullz, fmeanz = get_ensemble_sig(fullxray, sigtb5, binedge,posvar, aper=5)

#### Get non ensemble results ###
#numobs = np.shape(fluxnorm)[0]
#meanflux = np.nanmean(fluxnorm, axis=1)
#out = np.array([vari_funcs.maximum_likelihood(fluxnorm[n,:], fluxerrnorm[n,:], meanflux[n], posvar) for n in range(numobs)])
#
#numobs = np.shape(chanfluxnorm)[0]
#meanchan = np.nanmean(chanfluxnorm, axis=1)
#chanout = np.array([vari_funcs.maximum_likelihood(chanfluxnorm[n,:], chanerrnorm[n,:], meanchan[n], posvar) for n in range(numobs)])
#
#numobs = np.shape(fullfluxnorm)[0]
#meanfull = np.nanmean(fullfluxnorm, axis=1)
#fullout = np.array([vari_funcs.maximum_likelihood(fullfluxnorm[n,:], fullerrnorm[n,:], meanfull[n], posvar) for n in range(numobs)])

#%% Plot results ###

plt.figure(figsize=[8,7])
def plot_sig(tbdata, sigtb, binedge, aper=5, format='o',colour='b', label='Ensemble'):
    posvar = np.linspace(0,2,5000)
    fluxnorm, fluxerrnorm, sig, sigerr, z, meanz = get_ensemble_sig(tbdata, sigtb, binedge,posvar, aper)
    binupper = np.append(binedge[1:],np.max(z))
    xlow = meanz-binedge
    xhigh = binupper - meanz
    plt.errorbar(meanz, sig, xerr=[xlow, xhigh], yerr=sigerr, fmt=format,
                 color=colour, label=label)

nondata = tbdata[tbdata['X-ray']==False] # to get non X-ary variables
### plot for 3" ###
plot_sig(tbdata, sigtb5, binedge,aper=5,format = 'o',
                 colour='C0', label='3"')#label='Ensemble Variable Sources')
#plot_sig(chandata, sigtb5, binedge,aper=5,format = 'o',
#                 colour='r')#, label='Ensemble Variable X-ray Sources')
#plot_sig(fullxray, sigtb5, binedge,aper=5,format = 'o',
#                 colour='k', label='Ensemble Non Variable X-ray Sources')

### plot for 2" ###
#plot_sig(tbdata, sigtb4, binedge,aper=4,format = '+',
#                 colour='C1',label='2"')# label='Ensemble Variable Sources')
#plot_sig(chandata, sigtb4, binedge,aper=4,format = '+',
#                 colour='r')#, label='Ensemble Variable X-ray Sources')
#plot_sig(fullxray, sigtb4, binedge,aper=4, format = '+',
#                 colour='k')#, label='Ensemble Non Variable X-ray Sources')

### plot for 1.5" ###
#plot_sig(tbdata, sigtb3, binedge,aper=3,format = 'd',
#                 colour='C2', label='1.5"')#label='Ensemble Variable Sources')
#plot_sig(chandata, sigtb3, binedge,aper=3,format = 'd',
#                 colour='r')#, label='Ensemble Variable X-ray Sources')
#plot_sig(fullxray, sigtb3, binedge,aper=3, format = 'd',
#                 colour='k')#, label='Ensemble Non Variable X-ray Sources')

### plot for 1" ###
#plot_sig(tbdata, sigtb2, binedge,aper=2,format = 's',
#                 colour='C3', label='1"')#label='Ensemble Variable Sources')
#plot_sig(chandata, sigtb2, binedge,aper=2,format = 's',
#                 colour='r')#, label='Ensemble Variable X-ray Sources')
#plot_sig(fullxray, sigtb2, binedge,aper=2, format = 's',
#                 colour='k')#, label='Ensemble Non Variable X-ray Sources')

### plot for 0.67" ###
plot_sig(tbdata, sigtb1, binedge,aper=1,format = 'x',
                 colour='C4', label='0.67"')#label='Ensemble Variable Sources')
#plot_sig(chandata, sigtb1, binedge,aper=1,format = 'x',
#                 colour='r')#, label='Ensemble Variable X-ray Sources')
#plot_sig(fullxray, sigtb1, binedge,aper=1, format = 'x',
#                 colour='k')#, label='Ensemble Non Variable X-ray Sources')
plt.xlabel('z')
plt.ylabel(r'$\sigma$')
#plt.ylim(ymin=-0.25,ymax=1.755)
plt.xlim(xmin=-0.25, xmax=4.5)
plt.legend()
plt.tight_layout()




















end = time.time()
print(end-start)

