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
sigtb = Table.read('sigma_tables/quad_epoch_sigma_table_extra_clean_no06_2arcsec.fits')

posvar = np.linspace(0,2,5000)

#%% create function to do things 
def get_luminosity_and_flux(tbdata, xmm=False):
#    
    ### Extract magnitude table and error table ###
    flux = vari_funcs.flux5_stacks(tbdata)
    flux, tbdata = vari_funcs.noneg(flux, tbdata)
#    tbdata = tbdata[np.nanmean(flux,axis=1)>1e4]
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

def make_ensemble(tbdata, xrayL, binedge):
    
    ### Extract magnitude table and error table ###
    flux = vari_funcs.flux5_stacks(tbdata)
    flux, tbdata = vari_funcs.noneg(flux, tbdata)
    flux, fluxerr, tbdata = vari_funcs.create_quad_error_array(sigtb, tbdata)
    
    ### Normalise ###
    fluxnorm, fluxerrnorm = vari_funcs.normalise_flux_and_errors(flux, fluxerr)
    
    ### Create dicts to save data into ###
    enflux = {}
    enfluxerr = {}
    enXrayL = {}
    
    ### loop over bins ###
    for m, enmin in enumerate(binedge):
        ### Isolate data needed ###
        mask1 = xrayL >= enmin
        if m != len(binedge)-1:
            mask2 = xrayL < binedge[m+1]
        else:
            mask2 = np.ones(len(mask1))
        enmask = mask1*mask2.astype(bool)
        enflux[m] = fluxnorm[enmask]
        enfluxerr[m] = fluxerrnorm[enmask]
        enXrayL[m] = xrayL[enmask]
        
    return enflux, enfluxerr, enXrayL

def get_ensemble_sig(enflux, enfluxerr, enXrayL, posvar):
    ### Create dicts to save data into ###
    size = len(enflux)
    sig = np.empty(size)
    sigerr = np.empty(size)
    meanxrayL = np.empty(size)
    for m in enflux:
        ### Combine into one flux curve per bin ###
        enfluxcurve = np.ravel(enflux[m])
        enfluxcurveerr = np.ravel(enfluxerr[m])
        
        ### Find max likelihood sig of curve ###
        [sig[m],sigerr[m]] = vari_funcs.maximum_likelihood(enfluxcurve, enfluxcurveerr, 1, posvar)

        ### find mean xrayL ###
        meanxrayL[m] = np.nanmean(enXrayL[m])
    
    return sig, sigerr, meanxrayL
        
def z_split(tbdata):
    z = tbdata['z_spec']#[mask]
    z[z==-1] = tbdata['z_p'][z==-1]
    tbhigh = tbdata[z > 1.6]
    tblow = tbdata[z <= 1.6]
    return tbhigh, tblow

#%% remove edges ###
chandata = vari_funcs.remove_edges(chandata)
xmmdata = vari_funcs.remove_edges(xmmdata)
fullxray = vari_funcs.remove_edges(fullxray)
#%% Split by z ###
xmmhigh, xmmlow = z_split(xmmdata)
chanhigh, chanlow = z_split(chandata)
fullhigh, fulllow = z_split(fullxray)
#%%
chanhigh, highchanxrayL, highchanfluxnorm, highchanerrnorm= get_luminosity_and_flux(chanhigh)
chanlow, lowchanxrayL, lowchanfluxnorm, lowchanerrnorm= get_luminosity_and_flux(chanlow)
xmmhigh, highxmmxrayL, highxmmfluxnorm, highxmmerrnorm = get_luminosity_and_flux(xmmhigh, xmm=True)
xmmlow, lowxmmxrayL, lowxmmfluxnorm, lowxmmerrnorm = get_luminosity_and_flux(xmmlow, xmm=True)
fullhigh, highfullxrayL, highfullfluxnorm, highfullerrnorm= get_luminosity_and_flux(fullhigh)
fulllow, lowfullxrayL, lowfullfluxnorm, lowfullerrnorm= get_luminosity_and_flux(fulllow)

#%% Split into bins according to X-ray luminosity ###
hxrayL = np.append(highxmmxrayL, highchanxrayL)
lxrayL = np.append(lowxmmxrayL, lowchanxrayL)
sortedhighxrayL = np.sort(hxrayL)
sortedlowxrayL = np.sort(lxrayL)
#sortedxrayL = np.sort(fullxrayL)

binsize = 8
highbinedge = np.empty(int(len(hxrayL)/binsize))
lowbinedge = np.empty(int(len(lxrayL)/binsize))
size = len(highbinedge)
for n in range(size):
    ### Define bin edges ###
    enXrayL = sortedhighxrayL[n*binsize:(n*binsize)+binsize]
    enmin = np.nanmin(enXrayL)
    highbinedge[n] = enmin
size = len(lowbinedge)
for n in range(size):
    ### Define bin edges ###
    enXrayL = sortedlowxrayL[n*binsize:(n*binsize)+binsize]
    enmin = np.nanmin(enXrayL)
    lowbinedge[n] = enmin
    
sortedhighfullxrayL = np.sort(highfullxrayL)
sortedlowfullxrayL = np.sort(lowfullxrayL)
highfullbinedge = np.empty(int(len(highfullxrayL)/binsize))
size = len(highfullbinedge)
for n in range(size):
    ### Define bin edges ###
    enXrayL = sortedhighfullxrayL[n*binsize:(n*binsize)+binsize]
    enmin = np.nanmin(enXrayL)
    highfullbinedge[n] = enmin
lowfullbinedge = np.empty(int(len(lowfullxrayL)/binsize))
size = len(lowfullbinedge)
for n in range(size):
    ### Define bin edges ###
    enXrayL = sortedlowfullxrayL[n*binsize:(n*binsize)+binsize]
    enmin = np.nanmin(enXrayL)
    lowfullbinedge[n] = enmin

highchanenflux, highchanenfluxerr, highchanenXrayL = make_ensemble(chanhigh, highchanxrayL, highbinedge)
lowchanenflux, lowchanenfluxerr, lowchanenXrayL = make_ensemble(chanlow, lowchanxrayL, lowbinedge)
highxmmenflux, highxmmenfluxerr, highxmmenXrayL = make_ensemble(xmmhigh, highxmmxrayL, highbinedge)
lowxmmenflux, lowxmmenfluxerr, lowxmmenXrayL = make_ensemble(xmmlow, lowxmmxrayL, lowbinedge)
highfullenflux, highfullenfluxerr, highfullenXrayL = make_ensemble(fullhigh, highfullxrayL, highfullbinedge)
lowfullenflux, lowfullenfluxerr, lowfullenXrayL = make_ensemble(fulllow, lowfullxrayL, lowfullbinedge)


### Combine chan and xmm ##
henflux = {}
henfluxerr = {}
henXrayL = {}
for m, enmin in enumerate(highbinedge):
    henflux[m] = np.vstack([highxmmenflux[m],highchanenflux[m]])
    henfluxerr[m] = np.vstack([highxmmenfluxerr[m],highchanenfluxerr[m]])
    henXrayL[m] = np.append(highxmmenXrayL[m],highchanenXrayL[m])
lenflux = {}
lenfluxerr = {}
lenXrayL = {}
for m, enmin in enumerate(lowbinedge):
    lenflux[m] = np.vstack([lowxmmenflux[m],lowchanenflux[m]])
    lenfluxerr[m] = np.vstack([lowxmmenfluxerr[m],lowchanenfluxerr[m]])
    lenXrayL[m] = np.append(lowxmmenXrayL[m],lowchanenXrayL[m])

### Get ensemble maximum likelihood ###
highsig, highsigerr, highmeanxrayL = get_ensemble_sig(henflux, henfluxerr, henXrayL, posvar)
lowsig, lowsigerr, lowmeanxrayL = get_ensemble_sig(lenflux, lenfluxerr, lenXrayL, posvar)
highfullsig, highfullsigerr, highfullmeanxrayL = get_ensemble_sig(highfullenflux, highfullenfluxerr, highfullenXrayL, posvar)
lowfullsig, lowfullsigerr, lowfullmeanxrayL = get_ensemble_sig(lowfullenflux, lowfullenfluxerr, lowfullenXrayL, posvar)

#%% Plot results ###
highbinupper = np.append(highbinedge[1:],np.max(sortedhighxrayL))
lowbinupper = np.append(lowbinedge[1:],np.max(sortedlowxrayL))
highfullbinupper = np.append(highfullbinedge[1:],np.max(sortedhighfullxrayL))
lowfullbinupper = np.append(lowfullbinedge[1:],np.max(sortedlowfullxrayL))
highxlow = highmeanxrayL-highbinedge
highxhigh = highbinupper - highmeanxrayL
lowxlow = lowmeanxrayL-lowbinedge
lowxhigh = lowbinupper - lowmeanxrayL
highfullxlow = highfullmeanxrayL-highfullbinedge
highfullxhigh = highfullbinupper - highfullmeanxrayL
lowfullxlow = lowfullmeanxrayL-lowfullbinedge
lowfullxhigh = lowfullbinupper - lowfullmeanxrayL

plt.figure(figsize=[12,7])
#plt.errorbar(fmeanxrayL, fsig, xerr=[fxlow, fxhigh], yerr=fsigerr, fmt='o',
#             color='k', label='Ensemble Non Variable X-ray Sources')
plt.subplot(121)
plt.errorbar(lowmeanxrayL, lowsig, xerr=[lowxlow, lowxhigh], yerr=lowsigerr, fmt='o',
             color='b', label='z <= 1.6')
plt.errorbar(lowfullmeanxrayL, lowfullsig, xerr=[lowfullxlow, lowfullxhigh], 
             yerr=lowfullsigerr, fmt='o',color='tab:grey', label='Non Variable X-ray')
plt.xscale('log')
plt.xlabel('X-ray Luminosity (ergs/s)')
plt.ylabel(r'$\sigma$')
plt.ylim(ymin=-0.05,ymax=0.4)
plt.xlim(xmin=7e40,xmax=7e45)
plt.legend()
plt.subplot(122)
plt.errorbar(highmeanxrayL, highsig, xerr=[highxlow, highxhigh], yerr=highsigerr, fmt='o',
             color='g', label='z > 1.6')
plt.errorbar(highfullmeanxrayL, highfullsig, xerr=[highfullxlow, highfullxhigh], 
             yerr=highfullsigerr, fmt='o',color='tab:grey', label='Non Variable X-ray')
plt.xscale('log')
plt.xlabel('X-ray Luminosity (ergs/s)')
plt.ylabel(r'$\sigma$')
plt.ylim(ymin=-0.05,ymax=0.4)
plt.xlim(xmin=7e40,xmax=7e45)
plt.legend()
plt.tight_layout()





















end = time.time()
print(end-start)

