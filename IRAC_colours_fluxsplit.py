#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 16 10:47:04 2019

@author: ppxee
"""

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

    
### Get data ###
tbdata = fits.open('UDS_catalogues/DR11-2arcsec-June24-2018+plusXY_best_spec_z.fits')[1].data
varydata = fits.open('variable_tables/no06_variables_chi30_2arcsec_noXray_DR11data_restframe.fits')[1].data
#devdata = fits.open('variable_tables/no06_variables_chi30_2arcsec_deviant_DR11data_restframe.fits')[1].data
devdata = fits.open('variable_tables/no06_variables_chi30_2arcsec_DR11data_SpUDSdata.fits')[1].data
xvarydata = fits.open('variable_tables/no06_variables_chi30_2arcsec_Xray_DR11data_restframe.fits')[1].data

#def get_data(tbdata):
#    ### get IRAC colours ### 
#    I1 = tbdata['I1MAG_20']
#    I2 = tbdata['I2MAG_20']
#    
#    ### remove 99s ###
#    I1[I1 == 99] = np.nan # remove those with null values
#    I1[I1 == -99] = np.nan # remove those with null values
#    mask1 = ~np.isnan(I1)
#    I2[I2 == 99] = np.nan # remove those with null values
#    I2[I2 == -99] = np.nan # remove those with null values
#    mask2 = ~np.isnan(I2)
#    mask = mask1*mask2.astype('bool')
#    I1 = I1[mask]
#    I2 = I2[mask]
#    
#    ### Get z ###
#    z = tbdata['z_spec']
#    z[z==-1] = tbdata['z_p'][z==-1]
#    z = z[mask]
#    
#    return I1, I2, z

def get_data(tbdata):
    ### get IRAC colours ### 
    I1 = tbdata['I1MAG_20']
    I2 = tbdata['I2MAG_20']
    
    ### remove 99s ###
    I1[I1 == 99] = np.nan # remove those with null values
    I1[I1 == -99] = np.nan # remove those with null values
    mask1 = ~np.isnan(I1)
    I2[I2 == 99] = np.nan # remove those with null values
    I2[I2 == -99] = np.nan # remove those with null values
    mask2 = ~np.isnan(I2)
    mask = mask1*mask2.astype('bool')
    I1 = I1[mask]
    I2 = I2[mask]
    
    I1minI2 = I1 - I2
    
    ### Get z ###
    z = tbdata['z_spec']
    z[z==-1] = tbdata['z_p'][z==-1]
    z = z[mask]
    
    ### find ones out of range of the plot ###
    zlimsmask = z>4
    Iuplimsmask = I1minI2 > 1
    Ilolimsmask = I1minI2 < -0.7
#    Ilims = I1minI2[mask]
    
    return I1minI2, z, zlimsmask, Iuplimsmask, Ilolimsmask

### remove those with very high phot z ###
#xvarydata = xvarydata[xvarydata['z_p']<6]
#varydata = varydata[varydata['z_p']<6]

xvarydataup = vari_funcs.flux_split(xvarydata, 'upper')
varydataup = vari_funcs.flux_split(varydata, 'upper')

xvarydatalow = vari_funcs.flux_split(xvarydata, 'lower')
varydatalow = vari_funcs.flux_split(varydata, 'lower')

I1minI2, z, zlimsmask, Iuplimsmask, Ilolimsmask= get_data(tbdata)
varyI1minI2, varyz, varyzlimsmask, varyIuplimsmask, varyIlolimsmask = get_data(varydata)
xvaryI1minI2, xvaryz, xvaryzlimsmask, xvaryIuplimsmask, xvaryIlolimsmask = get_data(xvarydata)
devI1minI2, devz, devzlimsmask, devIuplimsmask, devIlolimsmask = get_data(devdata)

varyI1minI2up, varyzup, varyzlimsmaskup, varyIuplimsmaskup, varyIlolimsmaskup = get_data(varydataup)
xvaryI1minI2up, xvaryzup, xvaryzlimsmaskup, xvaryIuplimsmaskup, xvaryIlolimsmaskup = get_data(xvarydataup)

varyI1minI2low, varyzlow, varyzlimsmasklow, varyIuplimsmasklow, varyIlolimsmasklow = get_data(varydatalow)
xvaryI1minI2low, xvaryzlow, xvaryzlimsmasklow, xvaryIuplimsmasklow, xvaryIlolimsmasklow = get_data(xvarydatalow)


### Plot result ###
plt.figure(figsize=[8,8])
#plt.figure(figsize=[9,8])

### plot those within range
plt.plot(z, I1minI2,'k.', markersize=1, label='Galaxy')
plt.plot(varyz, varyI1minI2,'bo', label='Non X-ray Variable')
plt.plot(xvaryz, xvaryI1minI2,'ro', label='X-ray Variable')

#plt.plot(varyzup, varyI1minI2up,'go', label='High flux, not X-ray')
#plt.plot(xvaryzup, xvaryI1minI2up,'gd', label='High flux, X-ray')
#
#plt.plot(varyzlow, varyI1minI2low,'mo', label='Low flux, not X-ray')
#plt.plot(xvaryzlow, xvaryI1minI2low,'md', label='Low flux, X-ray')
#plt.plot(devz, devI1minI2,'ms',markersize=10,mfc='None', label='SpUDS Variable')
plt.xlim(xmin=0, xmax=4.1)         
plt.ylim(ymin=-0.7, ymax=1)
#plt.xlim(xmin=-0.1,xmax=10)
#plt.ylim(ymin=-1, ymax=1)
plt.xlabel('redshift')
plt.ylabel('3.6 - 4.5')
plt.legend(loc='upper right', frameon=False)
plt.tight_layout()

### plot out of range z ###
z[zlimsmask]=4
varyz[varyzlimsmask]=4
xvaryz[xvaryzlimsmask]=4
varyzup[varyzlimsmaskup]=4
xvaryzlow[xvaryzlimsmasklow]=4
varyzlow[varyzlimsmasklow]=4
xvaryz[xvaryzlimsmask]=4
devz[devzlimsmask]=4
plt.errorbar(z[zlimsmask], I1minI2[zlimsmask],fmt='k.',xerr=0.05, markersize=1, xlolims=True)
plt.errorbar(varyz[varyzlimsmask], varyI1minI2[varyzlimsmask],fmt='bo',xerr=0.05, xlolims=True)
plt.errorbar(xvaryz[xvaryzlimsmask], xvaryI1minI2[xvaryzlimsmask],fmt='ro',xerr=0.05, xlolims=True)

#plt.errorbar(varyzup[varyzlimsmaskup], varyI1minI2up[varyzlimsmaskup],fmt='go',xerr=0.05, xlolims=True)
#plt.errorbar(xvaryzup[xvaryzlimsmaskup], xvaryI1minI2up[xvaryzlimsmaskup],fmt='gd',xerr=0.05, xlolims=True)
#
#plt.errorbar(varyzlow[varyzlimsmasklow], varyI1minI2low[varyzlimsmasklow],fmt='mo',xerr=0.05, xlolims=True)
#plt.errorbar(xvaryzlow[xvaryzlimsmasklow], xvaryI1minI2low[xvaryzlimsmasklow],fmt='md',xerr=0.05, xlolims=True)

### plot out of range low I1minI2 ###
I1minI2[Ilolimsmask]=-0.65
#varyI1minI2[varyIlolimsmask]=-0.65
#xvaryI1minI2[xvaryIlolimsmask]=-0.65
varyI1minI2up[varyIlolimsmaskup]=-0.65
xvaryI1minI2up[xvaryIlolimsmaskup]=-0.65
varyI1minI2low[varyIlolimsmasklow]=-0.65
xvaryI1minI2low[xvaryIlolimsmasklow]=-0.65
devI1minI2[devIlolimsmask]=-0.65
plt.errorbar(z[Ilolimsmask], I1minI2[Ilolimsmask],fmt='k.',yerr=0.03, markersize=1, uplims=True)
plt.errorbar(varyz[varyIlolimsmask], varyI1minI2[varyIlolimsmask],fmt='bo',yerr=0.03, uplims=True)
plt.errorbar(xvaryz[xvaryIlolimsmask], xvaryI1minI2[xvaryIlolimsmask],fmt='ro',yerr=0.03, uplims=True)

#plt.errorbar(varyzup[varyIlolimsmaskup], varyI1minI2up[varyIlolimsmaskup],fmt='go',yerr=0.03, uplims=True)
#plt.errorbar(xvaryzup[xvaryIlolimsmaskup], xvaryI1minI2up[xvaryIlolimsmaskup],fmt='gd',yerr=0.03, uplims=True)
#
#plt.errorbar(varyzlow[varyIlolimsmasklow], varyI1minI2low[varyIlolimsmasklow],fmt='mo',yerr=0.03, uplims=True)
#plt.errorbar(xvaryzlow[xvaryIlolimsmasklow], xvaryI1minI2low[xvaryIlolimsmasklow],fmt='md',yerr=0.03, uplims=True)

### plot out of range high I1minI2 ###
I1minI2[Iuplimsmask]=0.95
#varyI1minI2[varyIuplimsmask]=0.95
#xvaryI1minI2[xvaryIuplimsmask]=0.95
varyI1minI2up[varyIuplimsmaskup]=0.95
xvaryI1minI2up[xvaryIuplimsmaskup]=0.95
varyI1minI2low[varyIuplimsmasklow]=0.95
xvaryI1minI2low[xvaryIuplimsmasklow]=0.95
devI1minI2[devIuplimsmask]=0.95
plt.errorbar(z[Iuplimsmask], I1minI2[Iuplimsmask],fmt='k.',yerr=0.03, markersize=1, lolims=True)
plt.errorbar(varyz[varyIuplimsmask], varyI1minI2[varyIuplimsmask],fmt='bo',yerr=0.03, lolims=True)
plt.errorbar(xvaryz[xvaryIuplimsmask], xvaryI1minI2[xvaryIuplimsmask],fmt='ro',yerr=0.03, lolims=True)

#plt.errorbar(varyzup[varyIuplimsmaskup], varyI1minI2up[varyIuplimsmaskup],fmt='go',yerr=0.03, lolims=True)
#plt.errorbar(xvaryzup[xvaryIuplimsmaskup], xvaryI1minI2up[xvaryIuplimsmaskup],fmt='gd',yerr=0.03, lolims=True)
#
#plt.errorbar(varyzlow[varyIuplimsmasklow], varyI1minI2low[varyIuplimsmasklow],fmt='mo',yerr=0.03, lolims=True)
#plt.errorbar(xvaryzlow[xvaryIuplimsmasklow], xvaryI1minI2low[xvaryIuplimsmasklow],fmt='md',yerr=0.03, lolims=True)
