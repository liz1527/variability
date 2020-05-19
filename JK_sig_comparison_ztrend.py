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
#xvarydata = varydata[varydata['X-ray']==True]
#noxvarydata = varydata[varydata['X-ray']==False]

xvarydata = Table.read('variable_tables/matched_Xray_J_and_K_variables_varystats_DR11data_0.25_0.5_0.2.fits')
noxvarydata = Table.read('variable_tables/matched_notXray_J_and_K_variables_varystats_DR11data_0.25_0.5_0.2.fits')

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

### calculate errors ###
differr = np.sqrt(np.square(Jouterr) + np.square(Kouterr))
xdifferr = np.sqrt(np.square(xJouterr) + np.square(xKouterr))
noxdifferr = np.sqrt(np.square(noxJouterr) + np.square(noxKouterr))
allxdifferr = np.sqrt(np.square(allxJouterr) + np.square(allxKouterr))

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
fig = plt.figure()
ax1 = fig.add_subplot(111) # Subplot covering whole plot
ax2 = ax1.twiny() # Twin of the first subplot
#plt.plot(allxz, allxdiff, 'k+')
plt.errorbar(noxz, noxdiff, yerr=noxdifferr, color='tab:grey', fmt='.',
             alpha=0.5, zorder=0)
plt.errorbar(xz, xdiff, yerr=xdifferr, color='tab:grey', fmt='.', alpha=0.5, 
             zorder=0)
plt.plot(noxz, noxdiff, 'bo')
plt.plot(xz, xdiff, 'ro')
plt.errorbar(noxz[noxz==4.35], noxdiff[noxz==4.35], xerr=0.09, fmt='b.', 
             zorder=0, xlolims=True)
plt.errorbar(xz[xz==4.35], xdiff[xz==4.35], xerr=0.09, fmt='r.', 
             zorder=0, xlolims=True)
#plt.errorbar(allxz[allxz==4.35], allxdiff[allxz==4.35], xerr=0.09, fmt='k+', 
#             zorder=0, xlolims=True)
ax1.set_xlabel('Redshift')
ax1.set_ylabel('Difference between $\sigma_{J}$ and $\sigma_{K}$')
ax1.set_xlim(-0.1, 4.5)
plt.hlines(0,-0.1, 4.5)

J_wave = 1.2
K_wave = 2.2
tick_z = np.array([0,1,2,3,4])
tick_J = np.round(J_wave/(1+tick_z),1) 
ax2.set_xlim(ax1.get_xlim()) #set twin to same limits
ax1.set_xticks(tick_z) #set z ticks
ax2.set_xticks(tick_z) #set z ticks
ax2.set_xticklabels(tick_J) #set J wave ticks
ax2.set_xlabel('J Restframe Wavelength ($\mu m$)')

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


end = time.time()
print(end-start)












