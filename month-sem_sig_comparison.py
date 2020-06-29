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
#varydata = Table.read('variable_tables/J_and_K_variables_varystats_DR11data.fits')
#xvarydata = varydata[varydata['X-ray']==True]
#noxvarydata = varydata[varydata['X-ray']==False]

### Import matched sample ###
varydata = Table.read('variable_tables/J_and_K_variables_month_varystats_DR11data.fits')


#### restrict to high flux ###
#monthflux = varydata['Month_Flux_K']
#meanflux = np.nanmean(monthflux, axis=1)
#varydata = varydata[meanflux>1e4]

xvarydata = varydata[varydata['X-ray']==True]
noxvarydata = varydata[varydata['X-ray']==False]

#### Find stellarity for full tbdata ##
#phot = Table.read('UDS_catalogues/DR11-2arcsec-Jun-30-2019_best_photometry.fits')
#small = phot['MAG_APER'][:,0] # 0.7 arcsec
#big = phot['MAG_APER'][:,3] # 2 arcsec
#stell = big - small
#
#### set values where mag == 99 to nan ### 
#small[small==99] = np.nan
#big[big==99] = np.nan
#
#mask = np.isin(phot['ID'], varydata['ID'])
#
#varystell = stell[mask]
#xstell = varystell[varydata['X-ray']==True]
#noxstell = varystell[varydata['X-ray']==False]


### Get data ###
#ksemout = varydata['sig_K']
#ksemouterr = varydata['sig_K_err']
kxsemout = xvarydata['sig_K']
kxsemouterr = xvarydata['sig_K_err']
knoxsemout = noxvarydata['sig_K']
knoxsemouterr = noxvarydata['sig_K_err']
#kmonthout = varydata['Month_sig_K']
#kmonthouterr = varydata['Month_sig_K_err']
kxmonthout = xvarydata['Month_sig_K']
kxmonthouterr = xvarydata['Month_sig_K_err']
knoxmonthout = noxvarydata['Month_sig_K']
knoxmonthouterr = noxvarydata['Month_sig_K_err']
#jsemout = varydata['sig_J']
#jsemouterr = varydata['sig_J_err']
jxsemout = xvarydata['sig_J']
jxsemouterr = xvarydata['sig_J_err']
jnoxsemout = noxvarydata['sig_J']
jnoxsemouterr = noxvarydata['sig_J_err']
#jmonthout = varydata['Month_sig_J']
#jmonthouterr = varydata['Month_sig_J_err']
jxmonthout = xvarydata['Month_sig_J']
jxmonthouterr = xvarydata['Month_sig_J_err']
jnoxmonthout = noxvarydata['Month_sig_J']
jnoxmonthouterr = noxvarydata['Month_sig_J_err']

### Indentify those with 0 sig in all, xray and nonxray ###
#inds = np.arange(0,len(varydata))
#zeroindssem = inds[varydata['sig_K']==0] #indicies where the sig_sem is zero
#zeroindsmonth = inds[varydata['Month_sig_K']==0] #indicies where the Month_sig_K is zero
#zeroinds = np.append(zeroindssem, zeroindsmonth)
#notzeroinds = inds[~np.isin(inds, zeroinds)]

kxinds = np.arange(0,len(xvarydata))
kxzeroindssem = kxinds[xvarydata['sig_K']==0] #indicies where the sig_sem is zero
kxzeroindsmonth = kxinds[xvarydata['Month_sig_K']==0] #indicies where the Month_sig_K is zero
kxzeroinds = np.append(kxzeroindssem, kxzeroindsmonth)
kxnotzeroinds = kxinds[~np.isin(kxinds, kxzeroinds)]

knoxinds = np.arange(0,len(noxvarydata))
knoxzeroindssem = knoxinds[noxvarydata['sig_K']==0] #indicies where the sig_sem is zero
knoxzeroindsmonth = knoxinds[noxvarydata['Month_sig_K']==0] #indicies where the Month_sig_K is zero
knoxzeroinds = np.append(knoxzeroindssem, knoxzeroindsmonth)
knoxnotzeroinds = knoxinds[~np.isin(knoxinds, knoxzeroinds)]

jxinds = np.arange(0,len(xvarydata))
jxzeroindssem = jxinds[xvarydata['sig_K']==0] #indicies where the sig_sem is zero
jxzeroindsmonth = jxinds[xvarydata['Month_sig_K']==0] #indicies where the Month_sig_K is zero
jxzeroinds = np.append(jxzeroindssem, jxzeroindsmonth)
jxnotzeroinds = jxinds[~np.isin(jxinds, jxzeroinds)]

jnoxinds = np.arange(0,len(noxvarydata))
jnoxzeroindssem = jnoxinds[noxvarydata['sig_J']==0] #indicies where the sig_sem is zero
jnoxzeroindsmonth = jnoxinds[noxvarydata['Month_sig_J']==0] #indicies where the Month_sig_K is zero
jnoxzeroinds = np.append(jnoxzeroindssem, jnoxzeroindsmonth)
jnoxnotzeroinds = jnoxinds[~np.isin(jnoxinds, jnoxzeroinds)]

### Find z for full tbdata ###
#z = varydata['z']
xz = xvarydata['z']
noxz = noxvarydata['z']

### edit 0 sigs so they shou on plot as empty circle upper limits ###
zerosigval = 1.5e-3
kxsemout[kxsemout==0] = zerosigval
kxmonthout[kxmonthout==0] = zerosigval
knoxsemout[knoxsemout==0] = zerosigval
knoxmonthout[knoxmonthout==0] = zerosigval
jxsemout[jxsemout==0] = zerosigval
jxmonthout[jxmonthout==0] = zerosigval
jnoxsemout[jnoxsemout==0] = zerosigval
jnoxmonthout[jnoxmonthout==0] = zerosigval


#### get flux data ###
#monthflux = varydata['Month_Flux_K']
#monthfluxerr = varydata['Month_Fluxerr_K']
#semflux = varydata['Flux_K']
#semfluxerr = varydata['Fluxerr_K']
#xmonthflux = xvarydata['Month_Flux_K']
#xmonthfluxerr = xvarydata['Month_Fluxerr_K']
#noxmonthflux = noxvarydata['Month_Flux_K']
#noxmonthfluxerr = noxvarydata['Month_Fluxerr_K']
#
#### calculate mean flux ###
#monthmeanflux = np.nanmean(monthflux, axis=1)
#semmeanflux = np.nanmean(semflux, axis=1)
#xmonthmeanflux = np.nanmean(xmonthflux, axis=1)
#noxmonthmeanflux = np.nanmean(noxmonthflux, axis=1)
#xsemmeanflux = np.nanmean(xmonthflux, axis=1)
#noxsemmeanflux = np.nanmean(noxmonthflux, axis=1)


x = np.linspace(0,2.5,10)
y = x
#%% Create K figure with X-ray/non X-ray split ###

plt.figure()
plt.errorbar(kxsemout[kxnotzeroinds], kxmonthout[kxnotzeroinds], 
             xerr=kxsemouterr[kxnotzeroinds], yerr=kxmonthouterr[kxnotzeroinds], fmt='.', 
             color='tab:grey', zorder=0, alpha=0.5)
plt.errorbar(knoxsemout[knoxnotzeroinds], knoxmonthout[knoxnotzeroinds], 
             xerr=knoxsemouterr[knoxnotzeroinds], yerr=knoxmonthouterr[knoxnotzeroinds], 
             fmt='.', color='tab:grey', zorder=0, alpha=0.5)
plt.plot(kxsemout, kxmonthout, 'rs', zorder=2)
plt.plot(knoxsemout, knoxmonthout, 'bo', zorder=1)

plt.errorbar(kxsemout[kxzeroindsmonth], kxmonthout[kxzeroindsmonth], 
             yerr=0.25e-3, fmt='r.', zorder=0, uplims=True)
plt.errorbar(knoxsemout[knoxzeroindsmonth], knoxmonthout[knoxzeroindsmonth], 
             yerr=0.25e-3, fmt='b.', zorder=0, uplims=True)
plt.errorbar(kxsemout[kxzeroindssem], kxmonthout[kxzeroindssem], xerr=0.25e-3, 
             fmt='r.', zorder=0, xuplims=True)
plt.errorbar(knoxsemout[knoxzeroindssem], knoxmonthout[knoxzeroindssem], 
             xerr=0.25e-3, fmt='b.', zorder=0, xuplims=True)

plt.xlabel('$\sigma_{sem}$')
plt.ylabel('$\sigma_{month}$')
plt.xscale('log')
plt.yscale('log')
plt.plot(x,y,'k')
plt.xlim(xmin=1e-3,xmax=2.3)
plt.ylim(ymin=1e-3,ymax=2.3)
plt.title('K-band')
plt.tight_layout()
#plt.savefig('plots/new_catalogue/JK_sig_comp/JK_sig_comp_matched.png')
#plt.savefig('plots/new_catalogue/JK_sig_comp/JK_sig_comp_K_variables.png')
#plt.savefig('plots/new_catalogue/JK_sig_comp/JK_sig_comp_J_variables.png')
#
#%% Create J figure with X-ray/non X-ray split ###

plt.figure()
plt.errorbar(jxsemout[jxnotzeroinds], jxmonthout[jxnotzeroinds], 
             xerr=jxsemouterr[jxnotzeroinds], yerr=jxmonthouterr[jxnotzeroinds], fmt='.', 
             color='tab:grey', zorder=0, alpha=0.5)
plt.errorbar(jnoxsemout[jnoxnotzeroinds], jnoxmonthout[jnoxnotzeroinds], 
             xerr=jnoxsemouterr[jnoxnotzeroinds], yerr=jnoxmonthouterr[jnoxnotzeroinds], 
             fmt='.', color='tab:grey', zorder=0, alpha=0.5)
plt.plot(jxsemout, jxmonthout, 'rs', zorder=2)
plt.plot(jnoxsemout, jnoxmonthout, 'bo', zorder=1)

plt.errorbar(jxsemout[jxzeroindsmonth], jxmonthout[jxzeroindsmonth], 
             yerr=0.25e-3, fmt='r.', zorder=0, uplims=True)
plt.errorbar(jnoxsemout[jnoxzeroindsmonth], jnoxmonthout[jnoxzeroindsmonth], 
             yerr=0.25e-3, fmt='b.', zorder=0, uplims=True)
plt.errorbar(jxsemout[jxzeroindssem], jxmonthout[jxzeroindssem], xerr=0.25e-3, 
             fmt='r.', zorder=0, xuplims=True)
plt.errorbar(jnoxsemout[jnoxzeroindssem], jnoxmonthout[jnoxzeroindssem], 
             xerr=0.25e-3, fmt='b.', zorder=0, xuplims=True)

plt.xlabel('$\sigma_{sem}$')
plt.ylabel('$\sigma_{month}$')
plt.xscale('log')
plt.yscale('log')
plt.plot(x,y,'k')
plt.xlim(xmin=1e-3,xmax=2.3)
plt.ylim(ymin=1e-3,ymax=2.3)
plt.title('J-band')
plt.tight_layout()
#plt.savefig('plots/new_catalogue/JK_sig_comp/JK_sig_comp_matched.png')
#plt.savefig('plots/new_catalogue/JK_sig_comp/JK_sig_comp_K_variables.png')
#plt.savefig('plots/new_catalogue/JK_sig_comp/JK_sig_comp_J_variables.png')
#
##%% Plot just the X-ray points with J-K colours on ###
#
#plt.figure()
#plt.errorbar(xKout[xnotzeroinds], xJout[xnotzeroinds], 
#             yerr=xJouterr[xnotzeroinds], xerr=xKouterr[xnotzeroinds], 
#             color='tab:grey', fmt='.', zorder=0, alpha=0.5)
#plt.scatter(xKout, xJout, c=x_J_K, vmin=-1, vmax=2.4, marker='s')
#
#
#plt.errorbar(xKout[xzeroindsJ], xJout[xzeroindsJ], yerr=0.25e-3, fmt='r.', zorder=0, 
#             uplims=True)
#plt.errorbar(xKout[xzeroindsK], xJout[xzeroindsK], xerr=0.25e-3, fmt='r.', zorder=0, 
#             xuplims=True)
#
#plt.xlabel('$\sigma_{K}$')
#plt.ylabel('$\sigma_{J}$')
#plt.title('X-ray Detected')
#plt.xscale('log')
#plt.yscale('log')
#cbar=plt.colorbar()
#cbar.set_label('J-K')
#plt.plot(x,y,'k')
#plt.xlim(xmin=1e-3,xmax=2.3)
#plt.ylim(ymin=1e-3,ymax=2.3)
#plt.tight_layout()
##plt.savefig('plots/new_catalogue/JK_sig_comp/JK_sig_comp_Xray_matched.png')
##plt.savefig('plots/new_catalogue/JK_sig_comp/JK_sig_comp_Xray_K_variables.png')
##plt.savefig('plots/new_catalogue/JK_sig_comp/JK_sig_comp_Xray_J_variables.png')
#
##%% Plot just the non-X-ray points with J-K colours on ###
#
#plt.figure()
#plt.errorbar(noxKout[noxnotzeroinds], noxJout[noxnotzeroinds], 
#             yerr=noxJouterr[noxnotzeroinds], xerr=noxKouterr[noxnotzeroinds], 
#             color='tab:grey', fmt='.', zorder=0, alpha=0.5)
#plt.scatter(noxKout, noxJout, c=nox_J_K, vmin=-1, vmax=2.4)
#
#plt.errorbar(noxKout[noxzeroindsJ], noxJout[noxzeroindsJ], yerr=0.25e-3, fmt='b.', 
#             zorder=0, uplims=True)
#plt.errorbar(noxKout[noxzeroindsK], noxJout[noxzeroindsK], xerr=0.25e-3, fmt='b.', 
#             zorder=0, xuplims=True)
#
#plt.xlabel('$\sigma_{K}$')
#plt.ylabel('$\sigma_{J}$')
#plt.xscale('log')
#plt.yscale('log')
#plt.title('Not X-ray Detected')
#cbar=plt.colorbar()
#cbar.set_label('J-K')
#plt.plot(x,y,'k')
#plt.xlim(xmin=1e-3,xmax=2.3)
#plt.ylim(ymin=1e-3,ymax=2.3)
#plt.tight_layout()
##plt.savefig('plots/new_catalogue/JK_sig_comp/JK_sig_comp_not_Xray_matched.png')
##plt.savefig('plots/new_catalogue/JK_sig_comp/JK_sig_comp_not_Xray_K_variables.png')
##plt.savefig('plots/new_catalogue/JK_sig_comp/JK_sig_comp_not_Xray_J_variables.png')
#
##%% Plot just the X-ray points with z colours on ###
#
#plt.figure()
#plt.errorbar(xKout[xnotzeroinds], xJout[xnotzeroinds], 
#             yerr=xJouterr[xnotzeroinds], xerr=xKouterr[xnotzeroinds], 
#             color='tab:grey', fmt='.', zorder=0, alpha=0.5)
#plt.scatter(xKout, xJout, c=xz, vmin=0, vmax=5, marker='s')
#
#
#plt.errorbar(xKout[xzeroindsJ], xJout[xzeroindsJ], yerr=0.25e-3, fmt='r.', zorder=0, 
#             uplims=True)
#plt.errorbar(xKout[xzeroindsK], xJout[xzeroindsK], xerr=0.25e-3, fmt='r.', zorder=0, 
#             xuplims=True)
#
#plt.xlabel('$\sigma_{K}$')
#plt.ylabel('$\sigma_{J}$')
#plt.title('X-ray Detected')
#plt.xscale('log')
#plt.yscale('log')
#cbar=plt.colorbar()
#cbar.set_label('z')
#plt.plot(x,y,'k')
#plt.xlim(xmin=1e-3,xmax=2.3)
#plt.ylim(ymin=1e-3,ymax=2.3)
#plt.tight_layout()
##plt.savefig('plots/new_catalogue/JK_sig_comp/JK_sig_comp_Xray_zcolours_matched.png')
##plt.savefig('plots/new_catalogue/JK_sig_comp/JK_sig_comp_Xray_zcolours_K_variables.png')
##plt.savefig('plots/new_catalogue/JK_sig_comp/JK_sig_comp_Xray_zcolours_J_variables.png')
#
##%% Plot just the non-X-ray points with z colours on ###
#
#plt.figure()
#plt.errorbar(noxKout[noxnotzeroinds], noxJout[noxnotzeroinds], 
#             yerr=noxJouterr[noxnotzeroinds], xerr=noxKouterr[noxnotzeroinds], 
#             color='tab:grey', fmt='.', zorder=0, alpha=0.5)
#plt.scatter(noxKout, noxJout, c=noxz, vmin=0, vmax=5)
#
#plt.errorbar(noxKout[noxzeroindsJ], noxJout[noxzeroindsJ], yerr=0.25e-3, fmt='b.', 
#             zorder=0, uplims=True)
#plt.errorbar(noxKout[noxzeroindsK], noxJout[noxzeroindsK], xerr=0.25e-3, fmt='b.', 
#             zorder=0, xuplims=True)
#
#plt.xlabel('$\sigma_{K}$')
#plt.ylabel('$\sigma_{J}$')
#plt.xscale('log')
#plt.yscale('log')
#plt.title('Not X-ray Detected')
#cbar=plt.colorbar()
#cbar.set_label('z')
#plt.plot(x,y,'k')
#plt.xlim(xmin=1e-3,xmax=2.3)
#plt.ylim(ymin=1e-3,ymax=2.3)
#plt.tight_layout()
##plt.savefig('plots/new_catalogue/JK_sig_comp/JK_sig_comp_not_Xray_zcolours_matched.png')
##plt.savefig('plots/new_catalogue/JK_sig_comp/JK_sig_comp_not_Xray_zcolours_K_variables.png')
##plt.savefig('plots/new_catalogue/JK_sig_comp/JK_sig_comp_not_Xray_zcolours_J_variables.png')

##%% Plot just the X-ray points with stellarity colours on ###
#
#### Change stell arrays so colour will be bi-modal instead of a spectrum ###
#def bimodal_colours(stell):
#    stell[stell < -1.2] = 0
#    stell[stell > -0.9] = 0
#    stell[stell != 0] = 1
#    return stell
#xstell = bimodal_colours(xstell)
#noxstell = bimodal_colours(noxstell)
#
#plt.figure()
#plt.errorbar(xKout[xnotzeroinds], xJout[xnotzeroinds], 
#             yerr=xJouterr[xnotzeroinds], xerr=xKouterr[xnotzeroinds], 
#             color='tab:grey', fmt='.', zorder=0, alpha=0.5)
#plt.scatter(xKout, xJout, c=xstell)#, vmin=-2, vmax=1.5)
#
#
#plt.errorbar(xKout[xzeroindsJ], xJout[xzeroindsJ], yerr=0.25e-3, fmt='r.', zorder=0, 
#             uplims=True)
#plt.errorbar(xKout[xzeroindsK], xJout[xzeroindsK], xerr=0.25e-3, fmt='r.', zorder=0, 
#             xuplims=True)
#
#plt.xlabel('$\sigma_{K}$')
#plt.ylabel('$\sigma_{J}$')
#plt.xscale('log')
#plt.yscale('log')
#plt.title('X-ray Detected')
#cbar=plt.colorbar()
##plt.clim(-2,1.5)
#cbar.set_label(r'$K_{2^{\prime\prime}} - K_{0.7^{\prime\prime}}$')
#plt.plot(x,y,'k')
#plt.xlim(xmin=1e-3,xmax=2.3)
#plt.ylim(ymin=1e-3,ymax=2.3)
#plt.tight_layout()
#plt.savefig('plots/new_catalogue/JK_sig_comp/JK_sig_comp_Xray_stellcolours_noxtalkcontam.png')
#
##%% Plot just the non-X-ray points with stellarity colours on ###
#
#plt.figure()
#plt.errorbar(noxKout[noxnotzeroinds], noxJout[noxnotzeroinds], 
#             yerr=noxJouterr[noxnotzeroinds], xerr=noxKouterr[noxnotzeroinds], 
#             color='tab:grey', fmt='.', zorder=0, alpha=0.5)
#plt.scatter(noxKout, noxJout, c=noxstell)#, vmin=-2, vmax=1.5)
#
#plt.errorbar(noxKout[noxzeroindsJ], noxJout[noxzeroindsJ], yerr=0.25e-3, fmt='b.', 
#             zorder=0, uplims=True)
#plt.errorbar(noxKout[noxzeroindsK], noxJout[noxzeroindsK], xerr=0.25e-3, fmt='b.', 
#             zorder=0, xuplims=True)
#
#plt.xlabel('$\sigma_{K}$')
#plt.ylabel('$\sigma_{J}$')
#plt.xscale('log')
#plt.yscale('log')
#plt.title('Not X-ray Detected')
#cbar=plt.colorbar()
##plt.clim(-2,0)
#cbar.set_label(r'$K_{2^{\prime\prime}} - K_{0.7^{\prime\prime}}$')
#plt.plot(x,y,'k')
#plt.xlim(xmin=1e-3,xmax=2.3)
#plt.ylim(ymin=1e-3,ymax=2.3)
#plt.tight_layout()
#plt.savefig('plots/new_catalogue/JK_sig_comp/JK_sig_comp_not_Xray_stellcolours_noxtalkcontam.png')

end = time.time()
print(end-start)












