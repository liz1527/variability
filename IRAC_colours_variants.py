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
varydata = fits.open('variable_tables/no06_variables_chi30_2arcsec_DR11data_restframe.fits')[1].data
devdata = fits.open('variable_tables/no06_variables_chi30_2arcsec_deviant_DR11data_restframe.fits')[1].data
noxvarydata = fits.open('variable_tables/no06_variables_chi30_2arcsec_noXray_DR11data_restframe.fits')[1].data
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
    
    return I1minI2, z, zlimsmask, Iuplimsmask, Ilolimsmask, tbdata[mask]

I1minI2, z, zlimsmask, Iuplimsmask, Ilolimsmask, tbdata= get_data(tbdata)
varyI1minI2, varyz, varyzlimsmask, varyIuplimsmask, varyIlolimsmask, varydata = get_data(varydata)
noxvaryI1minI2, noxvaryz, noxvaryzlimsmask, noxvaryIuplimsmask, noxvaryIlolimsmask, noxvarydata = get_data(noxvarydata)
xvaryI1minI2, xvaryz, xvaryzlimsmask, xvaryIuplimsmask, xvaryIlolimsmask, xvarydata = get_data(xvarydata)
devI1minI2, devz, devzlimsmask, devIuplimsmask, devIlolimsmask, devdata = get_data(devdata)

varyL = varydata['M_K_z_p']
xvaryL = xvarydata['M_K_z_p']
noxvaryL = noxvarydata['M_K_z_p']
clim = -17
### Plot result ###
plt.figure(figsize=[8,8])

### plot those within range
plt.plot(z, I1minI2,'k.', markersize=1, label='Galaxy', zorder=0)
plt.scatter(varyz, varyI1minI2, c=varyL, vmax=clim, marker='o', label='Variable')
plt.plot(noxvaryz, noxvaryI1minI2,'bo', markersize=10, mfc='None', label='Non X-ray Variable')
plt.plot(xvaryz, xvaryI1minI2,'ro', markersize=10, mfc='None', label='X-ray Variable')
#plt.plot(devz, devI1minI2,'mx', markersize=10, mfc='None', label='deviant')
plt.xlim(xmin=0, xmax=4.1)         
plt.ylim(ymin=-0.7, ymax=1)
#plt.xlim(xmin=-0.1,xmax=10)
#plt.ylim(ymin=-1, ymax=1)
plt.xlabel('redshift')
plt.ylabel('3.6 - 4.5')
plt.legend(loc='upper right', frameon=False)
cbar = plt.colorbar()
cbar.set_label('K Band Absolute Magnitude')
plt.tight_layout()

### plot out of range z ###
z[zlimsmask]=4
varyz[varyzlimsmask]=4
xvaryz[xvaryzlimsmask]=4
devz[devzlimsmask]=4
plt.errorbar(z[zlimsmask], I1minI2[zlimsmask],fmt='k.',xerr=0.05, markersize=1, 
             xlolims=True, zorder=0)
plt.errorbar(varyz[varyzlimsmask], varyI1minI2[varyzlimsmask],fmt='bo',
             xerr=0.05, markersize=10, mfc='None', xlolims=True, zorder=0)
plt.errorbar(xvaryz[xvaryzlimsmask], xvaryI1minI2[xvaryzlimsmask],fmt='ro',
             xerr=0.05, markersize=10, mfc='None', xlolims=True, zorder=0)
plt.errorbar(noxvaryz[noxvaryzlimsmask], noxvaryI1minI2[noxvaryzlimsmask],
             fmt='bo',xerr=0.05, markersize=10, mfc='None', xlolims=True, zorder=0)
plt.scatter(varyz[varyzlimsmask], varyI1minI2[varyzlimsmask], 
            c=varyL[varyzlimsmask], vmax=clim, marker='o')
plt.scatter(xvaryz[xvaryzlimsmask], xvaryI1minI2[xvaryzlimsmask], 
            c=xvaryL[xvaryzlimsmask], vmax=clim, marker='o')
plt.scatter(noxvaryz[noxvaryzlimsmask], noxvaryI1minI2[noxvaryzlimsmask], 
            c=noxvaryL[noxvaryzlimsmask], vmax=clim, marker='o')


### plot out of range low I1minI2 ###
I1minI2[Ilolimsmask]=-0.65
varyI1minI2[varyIlolimsmask]=-0.65
xvaryI1minI2[xvaryIlolimsmask]=-0.65
devI1minI2[devIlolimsmask]=-0.65
plt.errorbar(z[Ilolimsmask], I1minI2[Ilolimsmask],fmt='k.',yerr=0.03, 
             markersize=1, uplims=True, zorder=0)
plt.errorbar(varyz[varyIlolimsmask], varyI1minI2[varyIlolimsmask], fmt='bo',
             yerr=0.03, markersize=10, mfc='None', uplims=True, zorder=0)
plt.errorbar(noxvaryz[noxvaryIlolimsmask], noxvaryI1minI2[noxvaryIlolimsmask],
             fmt='bo', markersize=10, mfc='None', yerr=0.03, uplims=True, zorder=0)
plt.scatter(varyz[varyIlolimsmask], varyI1minI2[varyIlolimsmask], 
            c=varyL[varyIlolimsmask], vmax=clim, marker='o')
plt.scatter(noxvaryz[noxvaryIlolimsmask], noxvaryI1minI2[noxvaryIlolimsmask], 
            c=noxvaryL[noxvaryIlolimsmask], vmax=clim, marker='o')

### plot out of range high I1minI2 ###
I1minI2[Iuplimsmask]=0.95
varyI1minI2[varyIuplimsmask]=0.95
xvaryI1minI2[xvaryIuplimsmask]=0.95
devI1minI2[devIuplimsmask]=0.95
plt.errorbar(z[Iuplimsmask], I1minI2[Iuplimsmask],fmt='k.', yerr=0.03, 
             markersize=1, lolims=True, zorder=0)
plt.errorbar(varyz[varyIuplimsmask], varyI1minI2[varyIuplimsmask], fmt='bo',
             yerr=0.03, markersize=10, mfc='None', lolims=True, zorder=0)
plt.errorbar(noxvaryz[noxvaryIuplimsmask], noxvaryI1minI2[noxvaryIuplimsmask],
             fmt='bo', markersize=10, mfc='None', yerr=0.03, lolims=True)
plt.scatter(varyz[varyIuplimsmask], varyI1minI2[varyIuplimsmask], 
            c=varyL[varyIuplimsmask], vmax=clim, marker='o')
plt.scatter(noxvaryz[noxvaryIuplimsmask], noxvaryI1minI2[noxvaryIuplimsmask], 
            c=noxvaryL[noxvaryIuplimsmask], vmax=clim, marker='o')