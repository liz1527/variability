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
plt.close('all') #close any open plots
#from numpy.lib.recfunctions import append_fields

#%% Open the fits files and get data ###
### Import data for variables selected in K ###
KKdata = fits.open('variable_tables/no06_variables_chi30_2arcsec_DR11data_restframe.fits')[1].data
KJdata = fits.open('variable_tables/no06_variables_chi30_2arcsec_DR11data_J_best.fits')[1].data

### Import data for variables selected in J ###
JKdata = fits.open('variable_tables/J_variables_chi40_noneg_DR11data_Kdata.fits')[1].data
JJdata = fits.open('variable_tables/J_variables_chi40_noneg_DR11data.fits')[1].data

### Import sig data ###
Jsigtb = Table.read('sigma_tables/quad_epoch_sigma_table_extra_clean_2arcsec_J.fits')
Ksigtb = Table.read('sigma_tables/quad_epoch_sigma_table_extra_clean_no06_2arcsec_neg.fits')

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
noxKKout = np.array([])
xKKouterr = np.array([])
noxKKouterr = np.array([])
xKJout = np.array([])
noxKJout = np.array([])
xKJouterr = np.array([])
noxKJouterr = np.array([])
x_K_J_K = np.array([])
nox_K_J_K = np.array([])
zeroindsKK = np.array([])
zeroindsKJ = np.array([])
x_K_z = np.array([])
nox_K_z = np.array([])
KJinds = np.arange(len(KJdata))

### Find z for full tbdata ###
Kz = vari_funcs.get_z(KKdata)

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

    ### Plot separate points for X-ray and non X-ray detected ###
    if KJdata['X-ray'][KJmask] == True:
        
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
        
#        ### save x-ray flux for chan objects ###
#        xflux = KKdata[]
    else:
        
        ### save output into non-x-ray and band specfic arrays ###
        noxKKout = np.append(noxKKout, KKout[0])
        noxKJout = np.append(noxKJout, KJout[0])
        noxKKouterr = np.append(noxKKouterr, KKout[1])
        noxKJouterr = np.append(noxKJouterr, KJout[1])
        
        ### also save the J-K colour for that object ###
        J_K = KKdata['JMAG_20'][n] - KKdata['KMAG_20'][n]
        nox_K_J_K = np.append(nox_K_J_K, J_K)
        
        ### save z of object ###
        nox_K_z = np.append(nox_K_z, Kz[n])
        
### Set up arrays for J selected ###
xJKout = np.array([])
noxJKout = np.array([])
xJKouterr = np.array([])
noxJKouterr = np.array([])
xJJout = np.array([])
noxJJout = np.array([])
xJJouterr = np.array([])
noxJJouterr = np.array([])
x_J_J_K = np.array([])
nox_J_J_K = np.array([])
zeroindsJK = np.array([])
zeroindsJJ = np.array([])
x_J_z = np.array([])
nox_J_z = np.array([])
JKinds = np.arange(len(JKdata))

### Find z for full tbdata ###
Jz = vari_funcs.get_z(JJdata)

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
        
    ### Plot separate points for X-ray and non X-ray detected ###
    if JKdata['X-ray'][JKmask] == True:
        
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
        
    else:

        
        ### save output into non-x-ray and band specfic arrays ###
        noxJKout = np.append(noxJKout, JKout[0])
        noxJJout = np.append(noxJJout, JJout[0])
        noxJKouterr = np.append(noxJKouterr, JKout[1])
        noxJJouterr = np.append(noxJJouterr, JJout[1])
        
        ### also save the J-K colour for that object ###
        J_K = JJdata['JMAG_20'][n] - JJdata['KMAG_20'][n]
        nox_J_J_K = np.append(nox_J_J_K, J_K)
        
        ### save z of object ###
        nox_J_z = np.append(nox_J_z, Jz[n])
        
### edit 0 sigs so they shou on plot as empty circle upper limits ###
zerosigval = 1.5e-3
xKJout[xKJout==0] = zerosigval
xJKout[xJKout==0] = zerosigval
noxKJout[noxKJout==0] = zerosigval
noxJKout[noxJKout==0] = zerosigval

x = np.linspace(0,2.5,10)
y = x
#%% Create figure with X-ray/non X-ray split ###

plt.figure()
plt.errorbar(xKKout, xKJout, xerr=xKKouterr, yerr=xKJouterr, fmt='.', 
             color='tab:grey', zorder=0, alpha=0.5)
plt.errorbar(noxKKout, noxKJout, xerr=noxKKouterr, yerr=noxKJouterr, fmt='.', 
             color='tab:grey', zorder=0, alpha=0.5)
plt.plot(xKKout, xKJout, 'rs', zorder=2)
plt.plot(noxKKout, noxKJout, 'bo', zorder=1)
plt.errorbar(xJKout, xJJout, xerr=xJKouterr, yerr=xJJouterr, fmt='.', 
             color='tab:grey', zorder=0, alpha=0.5)
plt.errorbar(noxJKout, noxJJout, xerr=noxJKouterr, yerr=noxJJouterr, fmt='.', 
             color='tab:grey', zorder=0, alpha=0.5)
plt.plot(xJKout, xJJout, 'rs', zorder=2)
plt.plot(noxJKout, noxJJout, 'bo', zorder=1)

plt.errorbar(xKKout[xKJout==zerosigval], xKJout[xKJout==zerosigval], yerr=0.25e-3, fmt='r.', zorder=0, uplims=True)
plt.errorbar(noxKKout[noxKJout==zerosigval], noxKJout[noxKJout==zerosigval], yerr=0.25e-3, fmt='b.', zorder=0, uplims=True)
plt.errorbar(noxJKout[noxJKout==zerosigval], noxJJout[noxJKout==zerosigval], xerr=0.25e-3, fmt='b.', zorder=0, xuplims=True)
plt.errorbar(xJKout[xJKout==zerosigval], xJJout[xJKout==zerosigval], xerr=0.25e-3,fmt= 'r.', zorder=0, xuplims=True)

plt.xlabel('$\sigma_{K}$')
plt.ylabel('$\sigma_{J}$')
plt.xscale('log')
plt.yscale('log')
plt.plot(x,y,'k')
plt.xlim(xmin=1e-3,xmax=2.3)
plt.ylim(ymin=1e-3,ymax=2.3)
plt.tight_layout()
#plt.savefig('plots/JK_sig_comp/JK_sig_comp.png')
#plt.savefig('plots/JK_sig_comp/JK_sig_comp_K_variables.png')
#plt.savefig('plots/JK_sig_comp/JK_sig_comp_J_variables.png')

#%% Create figure with J/K selection split ###

plt.figure()
#plt.errorbar(xKKout, xKJout, xerr=xKKouterr, yerr=xKJouterr, fmt='.', 
#             color='tab:grey', zorder=0, alpha=0.5)
#plt.errorbar(noxKKout, noxKJout, xerr=noxKKouterr, yerr=noxKJouterr, fmt='.', 
#             color='tab:grey', zorder=0, alpha=0.5)
#plt.errorbar(xJKout, xJJout, xerr=xJKouterr, yerr=xJJouterr, fmt='.', 
#             color='tab:grey', zorder=0, alpha=0.5)
#plt.errorbar(noxJKout, noxJJout, xerr=noxJKouterr, yerr=noxJJouterr, fmt='.', 
#             color='tab:grey', zorder=0, alpha=0.5)
plt.plot(xKKout, xKJout, 'go')
plt.plot(xJKout, xJJout, 'ms', mfc='None', markersize=10)
plt.plot(noxKKout, noxKJout, 'go')
plt.plot(noxJKout, noxJJout, 'ms', mfc='None', markersize=10)

plt.errorbar(xKKout[xKJout==zerosigval], xKJout[xKJout==zerosigval], yerr=0.25e-3, fmt='r.', zorder=0, uplims=True)
plt.errorbar(noxKKout[noxKJout==zerosigval], noxKJout[noxKJout==zerosigval], yerr=0.25e-3, fmt='b.', zorder=0, uplims=True)
plt.errorbar(xJKout[xJKout==zerosigval], xJJout[xJKout==zerosigval], xerr=0.25e-3,fmt= 'r.', zorder=0, xuplims=True)
plt.errorbar(noxJKout[noxJKout==zerosigval], noxJJout[noxJKout==zerosigval], xerr=0.25e-3, fmt='b.', zorder=0, xuplims=True)

plt.xlabel('$\sigma_{K}$')
plt.ylabel('$\sigma_{J}$')
plt.xscale('log')
plt.yscale('log')
x = np.linspace(0,2,10)
y = x
plt.plot(x,y,'k')
plt.xlim(xmin=1e-3,xmax=2.3)
plt.ylim(ymin=1e-3,ymax=2.3)
plt.tight_layout()
plt.savefig('plots/JK_sig_comp/JK_sig_comp_selection_split.png')
#%% Plot just the X-ray points with J-K colours on ###

plt.figure()
plt.errorbar(xKKout, xKJout, yerr=xKJouterr, xerr=xKKouterr, color='tab:grey', 
             fmt='.', zorder=0, alpha=0.5)
plt.scatter(xKKout, xKJout, c=x_K_J_K, vmin=-1, vmax=2.4, marker='s')
plt.errorbar(xJKout, xJJout, yerr=xJJouterr, xerr=xJKouterr, color='tab:grey', 
             fmt='.', zorder=0, alpha=0.5)
plt.scatter(xJKout, xJJout, c=x_J_J_K, vmin=-1, vmax=2.4, marker='s')

plt.errorbar(xKKout[xKJout==zerosigval], xKJout[xKJout==zerosigval], yerr=0.25e-3, fmt='r.', zorder=0, uplims=True)
plt.errorbar(xJKout[xJKout==zerosigval], xJJout[xJKout==zerosigval], xerr=0.25e-3,fmt= 'r.', zorder=0, xuplims=True)

plt.xlabel('$\sigma_{K}$')
plt.ylabel('$\sigma_{J}$')
plt.title('X-ray Detected')
plt.xscale('log')
plt.yscale('log')
cbar=plt.colorbar()
cbar.set_label('J-K')
plt.plot(x,y,'k')
plt.xlim(xmin=1e-3,xmax=2.3)
plt.ylim(ymin=1e-3,ymax=2.3)
plt.tight_layout()
#plt.savefig('plots/JK_sig_comp/JK_sig_comp_Xray.png')
#plt.savefig('plots/JK_sig_comp/JK_sig_comp_Xray_K_variables.png')
#plt.savefig('plots/JK_sig_comp/JK_sig_comp_Xray_J_variables.png')

#%% Plot just the non-X-ray points with J-K colours on ###

plt.figure()
plt.errorbar(noxKKout[noxKJout!=0], noxKJout[noxKJout!=0], 
             yerr=noxKJouterr[noxKJout!=0], xerr=noxKKouterr[noxKJout!=0], 
             color='tab:grey', fmt='.', zorder=0, alpha=0.5)
plt.scatter(noxKKout, noxKJout, c=nox_K_J_K, vmin=-1, vmax=2.4)
plt.errorbar(noxJKout[noxJKout!=0], noxJJout[noxJKout!=0], 
             yerr=noxJJouterr[noxJKout!=0], xerr=noxJKouterr[noxJKout!=0], 
             color='tab:grey', fmt='.', zorder=0, alpha=0.5)
plt.scatter(noxJKout, noxJJout, c=nox_J_J_K, vmin=-1, vmax=2.4)

plt.errorbar(noxKKout[noxKJout==zerosigval], noxKJout[noxKJout==zerosigval], yerr=0.25e-3, fmt='b.', zorder=0, uplims=True)
plt.errorbar(noxJKout[noxJKout==zerosigval], noxJJout[noxJKout==zerosigval], xerr=0.25e-3, fmt='b.', zorder=0, xuplims=True)

plt.xlabel('$\sigma_{K}$')
plt.ylabel('$\sigma_{J}$')
plt.xscale('log')
plt.yscale('log')
plt.title('Not X-ray Detected')
cbar=plt.colorbar()
cbar.set_label('J-K')
plt.plot(x,y,'k')
plt.xlim(xmin=1e-3,xmax=2.3)
plt.ylim(ymin=1e-3,ymax=2.3)
plt.tight_layout()
#plt.savefig('plots/JK_sig_comp/JK_sig_comp_not_Xray.png')
#plt.savefig('plots/JK_sig_comp/JK_sig_comp_not_Xray_K_variables.png')
#plt.savefig('plots/JK_sig_comp/JK_sig_comp_not_Xray_J_variables.png')

#%% Plot just the X-ray points with z colours on ###

plt.figure()
plt.errorbar(xKKout, xKJout, yerr=xKJouterr, xerr=xKKouterr, color='tab:grey', 
             fmt='.', zorder=0, alpha=0.5)
plt.scatter(xKKout, xKJout, c=x_K_z, vmin=0, vmax=6, marker='s')
plt.errorbar(xJKout, xJJout, yerr=xJJouterr, xerr=xJKouterr, color='tab:grey', 
             fmt='.', zorder=0, alpha=0.5)
plt.scatter(xJKout, xJJout, c=x_J_z, vmin=0, vmax=6, marker='s')

plt.errorbar(xKKout[xKJout==zerosigval], xKJout[xKJout==zerosigval], yerr=0.25e-3, fmt='r.', zorder=0, uplims=True)
plt.errorbar(xJKout[xJKout==zerosigval], xJJout[xJKout==zerosigval], xerr=0.25e-3,fmt= 'r.', zorder=0, xuplims=True)

plt.xlabel('$\sigma_{K}$')
plt.ylabel('$\sigma_{J}$')
plt.title('X-ray Detected')
plt.xscale('log')
plt.yscale('log')
cbar=plt.colorbar()
cbar.set_label('z')
plt.plot(x,y,'k')
plt.xlim(xmin=1e-3,xmax=2.3)
plt.ylim(ymin=1e-3,ymax=2.3)
plt.tight_layout()
#plt.savefig('plots/JK_sig_comp/JK_sig_comp_Xray_zcolours.png')
#plt.savefig('plots/JK_sig_comp/JK_sig_comp_Xray_zcolours_K_variables.png')
#plt.savefig('plots/JK_sig_comp/JK_sig_comp_Xray_zcolours_J_variables.png')

#%% Plot just the non-X-ray points with z colours on ###

plt.figure()
plt.errorbar(noxKKout[noxKJout!=0], noxKJout[noxKJout!=0], 
             yerr=noxKJouterr[noxKJout!=0], xerr=noxKKouterr[noxKJout!=0], 
             color='tab:grey', fmt='.', zorder=0, alpha=0.5)
plt.scatter(noxKKout, noxKJout, c=nox_K_z, vmin=0, vmax=6)
plt.errorbar(noxJKout[noxJKout!=0], noxJJout[noxJKout!=0], 
             yerr=noxJJouterr[noxJKout!=0], xerr=noxJKouterr[noxJKout!=0], 
             color='tab:grey', fmt='.', zorder=0, alpha=0.5)
plt.scatter(noxJKout, noxJJout, c=nox_J_z, vmin=0, vmax=6)

plt.errorbar(noxKKout[noxKJout==zerosigval], noxKJout[noxKJout==zerosigval], yerr=0.25e-3, fmt='b.', zorder=0, uplims=True)
plt.errorbar(noxJKout[noxJKout==zerosigval], noxJJout[noxJKout==zerosigval], xerr=0.25e-3, fmt='b.', zorder=0, xuplims=True)

plt.xlabel('$\sigma_{K}$')
plt.ylabel('$\sigma_{J}$')
plt.xscale('log')
plt.yscale('log')
plt.title('Not X-ray Detected')
cbar=plt.colorbar()
cbar.set_label('z')
plt.plot(x,y,'k')
plt.xlim(xmin=1e-3,xmax=2.3)
plt.ylim(ymin=1e-3,ymax=2.3)
plt.tight_layout()
#plt.savefig('plots/JK_sig_comp/JK_sig_comp_not_Xray_zcolours.png')
#plt.savefig('plots/JK_sig_comp/JK_sig_comp_not_Xray_zcolours_K_variables.png')
#plt.savefig('plots/JK_sig_comp/JK_sig_comp_not_Xray_zcolours_J_variables.png')

#%% Save table of objects that have sigma_J == 0 ###
#zerostbdata = tbdata[zeroindsK.astype(int)] # creates full data table with K band lc data
#zerosJdata = Jdata[zeroindsJ.astype(int)] # creates table of J lcs
#
##tK = Table(zerostbdata)
##tK.write('variable_tables/zero_sigma_J_variables_K_data.fits')
##tJ = Table(zerosJdata)
##tJ.write('variable_tables/zero_sigma_J_variables_J_data.fits')

end = time.time()
print(end-start)












