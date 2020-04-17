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
varydata = Table.read('variable_tables/J_and_K_variables_varystats_DR11data.fits')
#varydata = varydata[varydata['KMAG_20']<=23]
xvarydata = varydata[varydata['X-ray']==True]
noxvarydata = varydata[varydata['X-ray']==False]

xraydata = Table.read('UDS_catalogues/xray_varystats_noneg.fits')

### remove variables from x-ray table ###
mask = np.isin(xraydata['ID'], varydata['ID'])
xraydata = xraydata[~mask]

### Get data ###
Kout = varydata['sig_K']
Kouterr = varydata['sig_K_err']
xKout = xvarydata['sig_K']
xKouterr = xvarydata['sig_K_err']
noxKout = noxvarydata['sig_K']
noxKouterr = noxvarydata['sig_K_err']
allxKout = xraydata['sig_K']
allxKouterr = xraydata['sig_K_err']
Jout = varydata['sig_J']
Jouterr = varydata['sig_J_err']
xJout = xvarydata['sig_J']
xJouterr = xvarydata['sig_J_err']
noxJout = noxvarydata['sig_J']
noxJouterr = noxvarydata['sig_J_err']
allxJout = xraydata['sig_J']
allxJouterr = xraydata['sig_J_err']
J_K = varydata['JMAG_20'] - varydata['KMAG_20']
x_J_K = xvarydata['JMAG_20'] - xvarydata['KMAG_20']
nox_J_K = noxvarydata['JMAG_20'] - noxvarydata['KMAG_20']

### Indentify those with 0 sig in all, xray and nonxray ###
inds = np.arange(0,len(varydata))
zeroindsK = inds[varydata['sig_K']==0] #indicies where the sig_K is zero
zeroindsJ = inds[varydata['sig_J']==0] #indicies where the sig_J is zero
zeroinds = np.append(zeroindsK, zeroindsJ)
notzeroinds = inds[~np.isin(inds, zeroinds)]

xinds = np.arange(0,len(xvarydata))
xzeroindsK = xinds[xvarydata['sig_K']==0] #indicies where the sig_K is zero
xzeroindsJ = xinds[xvarydata['sig_J']==0] #indicies where the sig_J is zero
xzeroinds = np.append(xzeroindsK, xzeroindsJ)
xnotzeroinds = xinds[~np.isin(xinds, xzeroinds)]

noxinds = np.arange(0,len(noxvarydata))
noxzeroindsK = noxinds[noxvarydata['sig_K']==0] #indicies where the sig_K is zero
noxzeroindsJ = noxinds[noxvarydata['sig_J']==0] #indicies where the sig_J is zero
noxzeroinds = np.append(noxzeroindsK, noxzeroindsJ)
noxnotzeroinds = noxinds[~np.isin(noxinds, noxzeroinds)]

allxinds = np.arange(0,len(xraydata))
allxzeroindsK = allxinds[xraydata['sig_K']==0] #indicies where the sig_K is zero
allxzeroindsJ = allxinds[xraydata['sig_J']==0] #indicies where the sig_J is zero
allxzeroinds = np.append(allxzeroindsK, allxzeroindsJ)
allxnotzeroinds = allxinds[~np.isin(allxinds, allxzeroinds)]

### Find z for full tbdata ###
z = varydata['z']
xz = xvarydata['z']
noxz = noxvarydata['z']
allxz = xraydata['z']


#%% Plot distance from 1-1 line in y vs z ###
### Calc distance from 1-1 line in sig_y ###
diff = Jout - Kout
xdiff = xJout - xKout
noxdiff = noxJout - noxKout
allxdiff = allxJout - allxKout

### fit lines to populations  ###
xfit = np.polyfit(xz, xdiff,1)
noxfit = np.polyfit(noxz, noxdiff,1)
allxfit = np.polyfit(allxz, allxdiff,1)

x = np.linspace(0,4.5,10)
yxfit = xfit[0] * x + xfit[1]
ynoxfit = noxfit[0] * x + noxfit[1]
yallxfit = allxfit[0] * x + allxfit[1]

### set values fro those off limits ###
z[z > 4.5] = 4.35
xz[xz > 4.5] = 4.35
noxz[noxz > 4.5] = 4.35
allxz[allxz > 4.5] = 4.35

### Plot ###
plt.figure()
#plt.plot(allxz, allxdiff, 'k+')
plt.plot(noxz, noxdiff, 'bo')
plt.plot(xz, xdiff, 'ro')
plt.errorbar(noxz[noxz==4.35], noxdiff[noxz==4.35], xerr=0.09, fmt='b.', 
             zorder=0, xlolims=True)
plt.errorbar(xz[xz==4.35], xdiff[xz==4.35], xerr=0.09, fmt='r.', 
             zorder=0, xlolims=True)
#plt.errorbar(allxz[allxz==4.35], allxdiff[allxz==4.35], xerr=0.09, fmt='k+', 
#             zorder=0, xlolims=True)
plt.xlabel('Redshift')
plt.ylabel('Difference between $\sigma_{J}$ and $\sigma_{K}$')
plt.xlim(-0.1, 4.5)
plt.hlines(0,-0.1, 4.5)
plt.tight_layout()

### Plot with fits ###
plt.figure()
#plt.plot(allxz, allxdiff, 'k+')
plt.plot(noxz, noxdiff, 'bo', alpha=0.25)
plt.plot(xz, xdiff, 'ro', alpha=0.25)
plt.errorbar(noxz[noxz==4.35], noxdiff[noxz==4.35], xerr=0.09, fmt='b.', 
             zorder=0, xlolims=True, alpha=0.5)
plt.errorbar(xz[xz==4.35], xdiff[xz==4.35], xerr=0.09, fmt='r.', 
             zorder=0, xlolims=True, alpha=0.2)
#plt.errorbar(allxz[allxz==4.35], allxdiff[allxz==4.35], xerr=0.09, fmt='k+', 
#             zorder=0, xlolims=True)
plt.plot(x, yxfit,'r-', label='Fit to vary X-ray')
plt.plot(x, ynoxfit,'b-', label='Fit to vary non-X-ray')
#plt.plot(x, yallxfit,'k--', label='Fit to X-ray')
plt.xlabel('Redshift')
plt.ylabel('Difference between $\sigma_{J}$ and $\sigma_{K}$')
plt.xlim(-0.1, 4.5)
plt.hlines(0,-0.1, 4.5)
plt.tight_layout()
### edit 0 sigs so they shou on plot as empty circle upper limits ###
zerosigval = 1.5e-3
xKout[xKout==0] = zerosigval
xJout[xJout==0] = zerosigval
noxKout[noxKout==0] = zerosigval
noxJout[noxJout==0] = zerosigval


#x = np.linspace(0,2.5,10)
#y = x
##%% Plot just the X-ray points with z colours on ###
#
#plt.figure()
#plt.errorbar(xKout[xnotzeroinds], xJout[xnotzeroinds], 
#             yerr=xJouterr[xnotzeroinds], xerr=xKouterr[xnotzeroinds], 
#             color='tab:grey', fmt='.', zorder=0, alpha=0.5)
#plt.scatter(xKout, xJout, c=xz, vmin=0, vmax=6, marker='s')
#
#plt.errorbar(xKout[xzeroindsK], xJout[xzeroindsK], yerr=0.25e-3, fmt='r.', 
#             zorder=0, uplims=True)
#plt.errorbar(xKout[xzeroindsJ], xJout[xzeroindsJ], xerr=0.25e-3, fmt='r.', 
#             zorder=0, xuplims=True)
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
##plt.savefig('plots/new_catalogue/JK_sig_comp/JK_sig_comp_Xray_zcolours.png')
##plt.savefig('plots/new_catalogue/JK_sig_comp/JK_sig_comp_Xray_zcolours_K_variables.png')
##plt.savefig('plots/new_catalogue/JK_sig_comp/JK_sig_comp_Xray_zcolours_J_variables.png')
#
##%% Plot just the non-X-ray points with z colours on ###
#
#plt.figure()
#plt.errorbar(noxKout[noxnotzeroinds], noxJout[noxnotzeroinds], 
#             yerr=noxJouterr[noxnotzeroinds], xerr=noxKouterr[noxnotzeroinds], 
#             color='tab:grey', fmt='.', zorder=0, alpha=0.5)
#plt.scatter(noxKout, noxJout, c=noxz, vmin=0, vmax=6)
#
#plt.errorbar(noxKout[noxzeroinds], noxJout[noxzeroinds], yerr=0.25e-3, fmt='b.', zorder=0, uplims=True)
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
##plt.savefig('plots/new_catalogue/JK_sig_comp/JK_sig_comp_not_Xray_zcolours.png')
##plt.savefig('plots/new_catalogue/JK_sig_comp/JK_sig_comp_not_Xray_zcolours_K_variables.png')
##plt.savefig('plots/new_catalogue/JK_sig_comp/JK_sig_comp_not_Xray_zcolours_J_variables.png')

end = time.time()
print(end-start)












