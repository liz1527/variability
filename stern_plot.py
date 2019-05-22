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
#tbdata = fits.open('UDS_catalogues/SpUDS_IRAC_catalogue_DR11data.fits')[1].data
tbdata = fits.open('UDS_catalogues/SpUDS_IRAC_catalogue.fits')[1].data
varydata = fits.open('variable_tables/no06_variables_chi30_2arcsec_DR11data_SpUDSdata.fits')[1].data
xvarydata = fits.open('variable_tables/no06_variables_chi30_2arcsec_Xray_DR11data_SpUDSdata.fits')[1].data

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
    I1 = tbdata['mag_ap_36']
    I2 = tbdata['mag_ap_45']
    
    
    ### get IRAC colours ### 
    I3 = tbdata['mag_ap_56']
    I4 = tbdata['mag_ap_80']
    
    ### remove 99s ###
    I1[I1 == 99] = np.nan # remove those with null values
    I1[I1 == -99] = np.nan # remove those with null values
    mask1 = ~np.isnan(I1)
    I2[I2 == 99] = np.nan # remove those with null values
    I2[I2 == -99] = np.nan # remove those with null values
    mask2 = ~np.isnan(I2)
#    mask = mask1*mask2.astype('bool')
    I3[I3 == 99] = np.nan # remove those with null values
    I3[I3 == -99] = np.nan # remove those with null values
    mask3 = ~np.isnan(I3)
    I4[I4 == 99] = np.nan # remove those with null values
    I4[I4 == -99] = np.nan # remove those with null values
    mask4 = ~np.isnan(I4)
    mask = mask1*mask2*mask3*mask4.astype('bool')
    I1 = I1[mask]
    I2 = I2[mask]
    I3 = I3[mask]
    I4 = I4[mask]
    
    ### Convert to vega ###
    I1 = I1 - 2.788
    I2 = I2 - 3.255
    I3 = I3 - 3.743
    I4 = I4 - 4.372
    
    ### Calculate differences ###
    I1minI2 = I1 - I2
    I3minI4 = I3 - I4
    
    ### find ones out of range of the plot ###
    I12uplimsmask = I1minI2 > 1.5
    I12lolimsmask = I1minI2 < -0.3
    I34uplimsmask = I3minI4 > 3.5
    I34lolimsmask = I3minI4 < -0.3
#    Ilims = I1minI2[mask]
    
    return I1minI2, I12uplimsmask, I12lolimsmask, I3minI4, I34uplimsmask, I34lolimsmask

I1minI2, I12uplimsmask, I12lolimsmask, I3minI4, I34uplimsmask, I34lolimsmask = get_data(tbdata)
varyI1minI2, varyI12uplimsmask, varyI12lolimsmask, varyI3minI4, varyI34uplimsmask, varyI34lolimsmask = get_data(varydata)
xvaryI1minI2, xvaryI12uplimsmask, xvaryI12lolimsmask, xvaryI3minI4, xvaryI34uplimsmask, xvaryI34lolimsmask = get_data(xvarydata)


### Plot result ###
plt.figure(figsize=[10,7])

### plot those within range
plt.plot(I3minI4, I1minI2,'k.', markersize=1, label='SpUDS source')
plt.plot(varyI3minI4, varyI1minI2,'bo', label='Variable SpUDS source')
plt.plot(xvaryI3minI4, xvaryI1minI2,'ro', label='Variable SpUDS source')
#plt.xlim(xmin=-1.2, xmax=2.8)         
#plt.ylim(ymin=-0.7, ymax=1)
plt.xlim(xmin=-0.3,xmax=3.5)
plt.ylim(ymin=-0.3, ymax=1.5)
plt.xlabel('5.6 - 8.0')
plt.ylabel('3.6 - 4.5')
plt.legend(loc='upper right', frameon=False)
plt.tight_layout()

### plot out of range low I1minI2 ###
I1minI2[I12lolimsmask]=-0.25
varyI1minI2[varyI12lolimsmask]=-0.25
plt.errorbar(I3minI4[I12lolimsmask], I1minI2[I12lolimsmask],fmt='k.',yerr=0.03, markersize=1, uplims=True)
plt.errorbar(varyI3minI4[varyI12lolimsmask], varyI1minI2[varyI12lolimsmask],fmt='bo',yerr=0.03, uplims=True)
#plt.errorbar(xvaryI3minI4[xvaryI12lolimsmask], xvaryI1minI2[xvaryI12lolimsmask],fmt='bo',yerr=0.03, uplims=True)


### plot out of range high I1minI2 ###
I1minI2[I12uplimsmask]= 1.45
varyI1minI2[varyI12uplimsmask]= 1.45
plt.errorbar(I3minI4[I12uplimsmask], I1minI2[I12uplimsmask],fmt='k.',yerr=0.03, markersize=1, lolims=True)
plt.errorbar(varyI3minI4[varyI12uplimsmask], varyI1minI2[varyI12uplimsmask],fmt='bo',yerr=0.03, lolims=True)


### plot out of range low I3minI4 ###
I1minI2[I12lolimsmask]=-0.25
varyI1minI2[varyI12lolimsmask]=-0.25
plt.errorbar(I3minI4[I12lolimsmask], I1minI2[I12lolimsmask],fmt='k.',yerr=0.03, markersize=1, uplims=True)
plt.errorbar(varyI3minI4[varyI12lolimsmask], varyI1minI2[varyI12lolimsmask],fmt='bo',yerr=0.03, uplims=True)


### plot out of range high I3minI4 ###
I1minI2[I34uplimsmask]= 3.45
varyI1minI2[varyI34uplimsmask]= 3.45
plt.errorbar(I3minI4[I34uplimsmask], I1minI2[I34uplimsmask],fmt='k.',yerr=0.03, markersize=1, lolims=True)
plt.errorbar(varyI3minI4[varyI34uplimsmask], varyI1minI2[varyI34uplimsmask],fmt='bo',yerr=0.03, lolims=True)

### Plot wedge ###
#plt.vlines(-0.029, -0.7, 1,linestyles='dashed')
#x = np.linspace(-0.029, 1.5)
#y1 = 0.2*(x-0.629) + 0.18
#y2 = 2.5*(x-0.629) - 3.5
plt.vlines(0.6, 0.3, 1.5,linestyles='dashed')
x1 = np.linspace(0.6, 1.6)
x2 = np.linspace(1.6, 2)
y1 = 0.2*(x1) + 0.18
y2 = 2.5*(x2) - 3.5
plt.plot(x1,y1,'k--')
plt.plot(x2,y2,'k--')
