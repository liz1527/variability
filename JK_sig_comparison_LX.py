#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  1 13:30:42 2018

Code to compute maximum likelihood variability

@author: ppxee
"""
import time
start = time.time()
print(start)

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
import matplotlib.colors as colors
plt.close('all') #close any open plots
#from numpy.lib.recfunctions import append_fields

#%% Open the fits files and get data ###
### Import data for variables selected in K ###
KKdata = fits.open('variable_tables/no06_variables_chi30_2arcsec_chandata_DR11data_restframe.fits')[1].data
KJdata = fits.open('variable_tables/no06_variables_chi30_2arcsec_DR11data_J_best.fits')[1].data

### Import data for variables selected in J ###
JKdata = fits.open('variable_tables/J_variables_chi40_noneg_DR11data_Kdata.fits')[1].data
JJdata = fits.open('variable_tables/J_variables_chi40_noneg_chanXray_DR11data.fits')[1].data

### Import sig data ###
Jsigtb = Table.read('sigma_tables/quad_epoch_sigma_table_extra_clean_2arcsec_J_noneg.fits')
Ksigtb = Table.read('sigma_tables/quad_epoch_sigma_table_extra_clean_no06_2arcsec.fits')

#Jxraydata = Jdata[Jdata['X-ray']==True]

#### Limit to Chandra region for simplicity ###
#KKdata = vari_funcs.field_funcs.chandra_only(KKdata)
#KJdata = vari_funcs.field_funcs.chandra_only(KJdata)
#JKdata = vari_funcs.field_funcs.chandra_only(JKdata)
#JJdata = vari_funcs.field_funcs.chandra_only(JJdata)

### Extract magnitude table and error tables ###
KKflux, KKfluxerr, KKdata = vari_funcs.k_mag_flux.create_quad_error_array(Ksigtb, KKdata, aper=4)
KJflux, KJfluxerr, KJdata = vari_funcs.j_mag_flux.create_quad_error_array_J(Jsigtb, KJdata, aper=4)
JKflux, JKfluxerr, JKdata = vari_funcs.k_mag_flux.create_quad_error_array(Ksigtb, JKdata, aper=4)
JJflux, JJfluxerr, JJdata = vari_funcs.j_mag_flux.create_quad_error_array_J(Jsigtb, JJdata, aper=4)

### Normalise ###
KKfluxnorm, KKfluxerrnorm = vari_funcs.flux_funcs.normalise_flux_and_errors(KKflux, KKfluxerr)
KJfluxnorm, KJfluxerrnorm = vari_funcs.flux_funcs.normalise_flux_and_errors(KJflux, KJfluxerr)
JKfluxnorm, JKfluxerrnorm = vari_funcs.flux_funcs.normalise_flux_and_errors(JKflux, JKfluxerr)
JJfluxnorm, JJfluxerrnorm = vari_funcs.flux_funcs.normalise_flux_and_errors(JJflux, JJfluxerr)

#%% All points
posvar = np.linspace(0,2,5000)
KKmeanflux = np.nanmean(KKfluxnorm, axis=1)
KJmeanflux = np.nanmean(KJfluxnorm, axis=1)
JKmeanflux = np.nanmean(JKfluxnorm, axis=1)
JJmeanflux = np.nanmean(JJfluxnorm, axis=1)

#start = time.time()

### Set up arrays for K selected ###
xKKout = np.array([])
xKKouterr = np.array([])
xKJout = np.array([])
xKJouterr = np.array([])
x_K_J_K = np.array([])
zeroindsKK = np.array([])
zeroindsKJ = np.array([])
x_K_z = np.array([])
x_K_LX = np.array([])
KJinds = np.arange(len(KJdata))

### Find z and x-ray lum for full tbdata ###
K_FX, K_LX, Kz = vari_funcs.xray_funcs.get_xray_L(KKdata)

for n in range(len(KKdata)): #loop over the selection band
    obnum = KKdata['NUMBER_1'][n] #get DR11 number
    KJmask = np.isin(KJdata['NUMBER_1'], obnum) #find equivilant J
    if ~np.any(KJmask):
        continue
    
    ### Get maximum likelihoods in J and K for that object ###
    KKout = vari_funcs.vary_stats.maximum_likelihood(KKfluxnorm[n,:], 
                                                     KKfluxerrnorm[n,:], 
                                                     KKmeanflux[n], posvar)
    
    KJout = vari_funcs.vary_stats.maximum_likelihood(KJfluxnorm[KJmask,:].reshape(8), 
                                                     KJfluxerrnorm[KJmask,:].reshape(8), 
                                                     KJmeanflux[KJmask], posvar)
    
    
    ### save index if sigma==0 in J ###
    if KJout[0] == 0:
        zeroindsKJ = np.append(zeroindsKJ, KJinds[KJmask])
        zeroindsKK = np.append(zeroindsKK, n)
        
    ### save output into x-ray and band specfic arrays ###
    xKKout = np.append(xKKout, KKout[0])
    xKJout = np.append(xKJout, KJout[0])
    xKKouterr = np.append(xKKouterr, KKout[1])
    xKJouterr = np.append(xKJouterr, KJout[1])
    
    ### also save the J-K colour for that object ###
    J_K = KKdata['JMAG_20'][n] - KKdata['KMAG_20'][n]
    x_K_J_K = np.append(x_K_J_K, J_K)
    
    ### save z of object ###
    x_K_z = np.append(x_K_z, Kz[n])
    
    ### Save X-ray lum ###
    x_K_LX = np.append(x_K_LX, K_LX[n])
        
### Set up arrays for J selected ###
xJKout = np.array([])
xJKouterr = np.array([])
xJJout = np.array([])
xJJouterr = np.array([])
x_J_J_K = np.array([])
zeroindsJK = np.array([])
zeroindsJJ = np.array([])
x_J_z = np.array([])
x_J_LX = np.array([])
JKinds = np.arange(len(JKdata))

### Find z and x-ray lum for full tbdata ###
J_FX, J_LX, Jz = vari_funcs.xray_funcs.get_xray_L(JJdata)

for n in range(len(JJdata)): #loop over the selection band
    obnum = JJdata['ID'][n] #find DR11 number
    JKmask = np.isin(JKdata['NUMBER_1'], obnum) #find that object in K array
    if ~np.any(JKmask):
        continue
    
    ### Get maximum likelihoods in J and K for that object ###
    JJout = vari_funcs.vary_stats.maximum_likelihood(JJfluxnorm[n,:], 
                                                     JJfluxerrnorm[n,:], 
                                                     JJmeanflux[n], posvar)
    
    JKout = vari_funcs.vary_stats.maximum_likelihood(JKfluxnorm[JKmask,:].reshape(7), 
                                                     JKfluxerrnorm[JKmask,:].reshape(7), 
                                                     JKmeanflux[JKmask], posvar)


    ### save index if sigma==0 in K ###
    if JKout[0] == 0:
        zeroindsJK = np.append(zeroindsJJ, JKinds[JKmask])
        zeroindsJJ = np.append(zeroindsJJ, n)
#        JKou
        
    ### save output into x-ray and band specfic arrays ###
    xJKout = np.append(xJKout, JKout[0])
    xJJout = np.append(xJJout, JJout[0])
    xJKouterr = np.append(xJKouterr, JKout[1])
    xJJouterr = np.append(xJJouterr, JJout[1])
    
    ### also save the J-K colour for that object ###
    J_K = JJdata['JMAG_20'][n] - JJdata['KMAG_20'][n]
    x_J_J_K = np.append(x_J_J_K, J_K)
    
    ### save z of object ###
    x_J_z = np.append(x_J_z, Jz[n])
    
    ### Save X-ray lum ###
    x_J_LX = np.append(x_J_LX, J_LX[n])
    
### edit 0 sigs so they shou on plot as empty circle upper limits ###
zerosigval = 1.5e-3
xKJout[xKJout==0] = zerosigval
xJKout[xJKout==0] = zerosigval
#%% Plot just the X-ray points with LX colours on ###
x = np.linspace(0,2.5,10)
y = x

plt.figure()
plt.errorbar(xKKout, xKJout, yerr=xKJouterr, xerr=xKKouterr, color='tab:grey', 
             fmt='.', zorder=0, alpha=0.5)
plt.scatter(xKKout, xKJout, c=x_K_LX.value, marker='s', norm=colors.LogNorm())
plt.errorbar(xJKout, xJJout, yerr=xJJouterr, xerr=xJKouterr, color='tab:grey', 
             fmt='.', zorder=0, alpha=0.5)
plt.scatter(xJKout, xJJout, c=x_J_LX.value, marker='s', norm=colors.LogNorm())

plt.errorbar(xKKout[xKJout==zerosigval], xKJout[xKJout==zerosigval], 
             yerr=0.25e-3, fmt='r.', zorder=0, uplims=True)
plt.errorbar(xJKout[xJKout==zerosigval], xJJout[xJKout==zerosigval], 
             xerr=0.25e-3,fmt= 'r.', zorder=0, xuplims=True)

plt.xlabel('$\sigma_{K}$')
plt.ylabel('$\sigma_{J}$')
plt.title('X-ray Detected')
plt.xscale('log')
plt.yscale('log')
cbar=plt.colorbar()
cbar.set_label(r'$L_{X}$')
plt.plot(x,y,'k')
plt.xlim(xmin=1e-3,xmax=2.3)
plt.ylim(ymin=1e-3,ymax=2.3)
plt.tight_layout()
#plt.savefig('plots/JK_sig_comp/JK_sig_comp_Xray.png')
#plt.savefig('plots/JK_sig_comp/JK_sig_comp_Xray_K_variables.png')
#plt.savefig('plots/JK_sig_comp/JK_sig_comp_Xray_J_variables.png')


end = time.time()
print(end-start)







