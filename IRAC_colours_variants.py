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
import matplotlib.colors as colors
plt.close('all') #close any open plots

    
### Get data ###
tbdata = fits.open('UDS_catalogues/DR11-2arcsec-June24-2018+plusXY_best_spec_z.fits')[1].data
varydata = fits.open('variable_tables/no06_variables_chi30_2arcsec_DR11data_restframe.fits')[1].data
devdata = fits.open('variable_tables/no06_variables_chi30_2arcsec_deviant_DR11data_restframe.fits')[1].data
noxvarydata = fits.open('variable_tables/no06_variables_chi30_2arcsec_noXray_DR11data_restframe.fits')[1].data
xvarydata = fits.open('variable_tables/no06_variables_chi30_2arcsec_Xray_DR11data_restframe.fits')[1].data
xdata = fits.open('UDS_catalogues/SpUDS_IRAC_catalogue_DR11data_xray.fits')[1].data

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
xI1minI2, xz, xzlimsmask, xIuplimsmask, xIlolimsmask, xdata= get_data(xdata)

### flux ###
varyL = vari_funcs.flux4_stacks(varydata)
varyL = np.nanmean(varyL, axis=1)
xvaryL = vari_funcs.flux4_stacks(xvarydata)
xvaryL = np.nanmean(xvaryL, axis=1)
noxvaryL = vari_funcs.flux4_stacks(noxvarydata)
noxvaryL = np.nanmean(noxvaryL, axis=1)
#colmap = 'cool' 
#
#### M_star ###
#varyL = varydata['Mstar_z_p']
#xvaryL = xvarydata['Mstar_z_p']
#noxvaryL = noxvarydata['Mstar_z_p']
#def ignore_zeros(m, z, I1minI2, zlimsmask, Iuplimsmask, Ilolimsmask, tbdata):
#    m[m==0] = np.nan
#    mask = ~np.isnan(m)
#    z = z[mask]
#    I1minI2 = I1minI2[mask]
#    m = m[mask]
#    zlimsmask = zlimsmask[mask]
#    Iuplimsmask = Iuplimsmask[mask]
#    Ilolimsmask = Ilolimsmask[mask]
#    tbdata = tbdata[mask]
#    return m, z, I1minI2, zlimsmask, Iuplimsmask, Ilolimsmask, tbdata
#
#varyL, varyz, varyI1minI2, varyzlimsmask, \
#varyIuplimsmask, varyIlolimsmask, varydata  = ignore_zeros(varyL, varyz, 
#                                                           varyI1minI2, 
#                                                           varyzlimsmask,
#                                                           varyIuplimsmask, 
#                                                           varyIlolimsmask, 
#                                                           varydata)
#xvaryL, xvaryz, xvaryI1minI2, xvaryzlimsmask, \
#xvaryIuplimsmask, xvaryIlolimsmask, xvarydata  = ignore_zeros(xvaryL, xvaryz, 
#                                                           xvaryI1minI2, 
#                                                           xvaryzlimsmask,
#                                                           xvaryIuplimsmask, 
#                                                           xvaryIlolimsmask, 
#                                                           xvarydata)
#noxvaryL, noxvaryz, noxvaryI1minI2, noxvaryzlimsmask, \
#noxvaryIuplimsmask, noxvaryIlolimsmask, noxvarydata  = ignore_zeros(noxvaryL, noxvaryz, 
#                                                           noxvaryI1minI2, 
#                                                           noxvaryzlimsmask,
#                                                           noxvaryIuplimsmask, 
#                                                           noxvaryIlolimsmask, 
#                                                           noxvarydata)
### Absolute magnitude ###
varyL = varydata['M_K_z_p']
xvaryL = xvarydata['M_K_z_p']
noxvaryL = noxvarydata['M_K_z_p']
def ignore_99s(m, z, I1minI2, zlimsmask, Iuplimsmask, Ilolimsmask, tbdata):
    m[m==99] = np.nan
    mask = ~np.isnan(m)
    z = z[mask]
    I1minI2 = I1minI2[mask]
    m = m[mask]
    zlimsmask = zlimsmask[mask]
    Iuplimsmask = Iuplimsmask[mask]
    Ilolimsmask = Ilolimsmask[mask]
    tbdata = tbdata[mask]
    return m, z, I1minI2, zlimsmask, Iuplimsmask, Ilolimsmask, tbdata

varyL, varyz, varyI1minI2, varyzlimsmask, \
varyIuplimsmask, varyIlolimsmask, varydata  = ignore_99s(varyL, varyz, 
                                                           varyI1minI2, 
                                                           varyzlimsmask,
                                                           varyIuplimsmask, 
                                                           varyIlolimsmask, 
                                                           varydata)
xvaryL, xvaryz, xvaryI1minI2, xvaryzlimsmask, \
xvaryIuplimsmask, xvaryIlolimsmask, xvarydata  = ignore_99s(xvaryL, xvaryz, 
                                                           xvaryI1minI2, 
                                                           xvaryzlimsmask,
                                                           xvaryIuplimsmask, 
                                                           xvaryIlolimsmask, 
                                                           xvarydata)
noxvaryL, noxvaryz, noxvaryI1minI2, noxvaryzlimsmask, \
noxvaryIuplimsmask, noxvaryIlolimsmask, noxvarydata  = ignore_99s(noxvaryL, noxvaryz, 
                                                           noxvaryI1minI2, 
                                                           noxvaryzlimsmask,
                                                           noxvaryIuplimsmask, 
                                                           noxvaryIlolimsmask, 
                                                           noxvarydata)

#clim = np.nanmax(varyL)
#
### Plot result ###
plt.figure(figsize=[8,8])
#
#### find max and min ###
cmax = -17#np.nanmax(varyL)
#cmin = np.nanmin(varyL)#1e8

### Impose axes limits ###
z[zlimsmask]=4
xz[xzlimsmask]=4
varyz[varyzlimsmask]=4
xvaryz[xvaryzlimsmask]=4
noxvaryz[noxvaryzlimsmask]=4
devz[devzlimsmask]=4

I1minI2[Ilolimsmask]=-0.65
xI1minI2[xIlolimsmask]=-0.65
varyI1minI2[varyIlolimsmask]=-0.65
xvaryI1minI2[xvaryIlolimsmask]=-0.65
noxvaryI1minI2[noxvaryIlolimsmask]=-0.65
devI1minI2[devIlolimsmask]=-0.65

I1minI2[Iuplimsmask]=0.95
xI1minI2[xIuplimsmask]=0.95
varyI1minI2[varyIuplimsmask]=0.95
xvaryI1minI2[xvaryIuplimsmask]=0.95
noxvaryI1minI2[noxvaryIuplimsmask]=0.95
devI1minI2[devIuplimsmask]=0.95

### plot background galaxies and X-ray sources ###
plt.plot(z, I1minI2,'.', markersize=1, color='tab:gray', label='Galaxy',
         alpha=0.5, zorder=0)
#plt.plot(xz, xI1minI2,'k+',label='X-ray source', zorder=0)

### plot x and non x as full circles ###
#plt.plot(noxvaryz, noxvaryI1minI2,'bo', label='Non X-ray Variable')
#plt.plot(xvaryz, xvaryI1minI2,'ro', label='X-ray Variable')

### Plot x and non x as boundary circles with colours from scatter ###
#plt.scatter(varyz, varyI1minI2, c=varyL, 
#            norm=colors.LogNorm(vmin=cmin, vmax=cmax), marker='o', label='Variable')
#plt.plot(noxvaryz, noxvaryI1minI2,'bo', mfc='None', label='Non X-ray Variable')
#plt.plot(xvaryz, xvaryI1minI2,'ro', mfc='None', label='X-ray Variable')

### Plot deviants ###
#plt.plot(devz, devI1minI2,'mx', markersize=10, mfc='None', label='deviant')

#### Plot x and non x as different shapes with colours - log scale ###
#plt.scatter(noxvaryz, noxvaryI1minI2, c=noxvaryL,
#            norm=colors.LogNorm(vmin=cmin, vmax=cmax), marker='x',
#            label='Non X-ray Variable')
#plt.scatter(xvaryz, xvaryI1minI2, c=xvaryL,
#            norm=colors.LogNorm(vmin=cmin, vmax=cmax), marker='o', #s=15,
#            label='X-ray Variable')

### Plot x and non x as different shapes with colours - linear scale ###
plt.scatter(noxvaryz, noxvaryI1minI2, c=noxvaryL, vmax=cmax, marker='x',
            label='Non X-ray Variable')
plt.scatter(xvaryz, xvaryI1minI2, c=xvaryL, vmax=cmax, marker='o', #s=15,
            label='X-ray Variable')

plt.xlim(xmin=0, xmax=4.1)         
plt.ylim(ymin=-0.7, ymax=1)
plt.xlabel('z')
plt.ylabel('3.6 - 4.5')
plt.legend(loc='upper right', frameon=False)
cbar = plt.colorbar()
cbar.set_label('K Band Absolute Magnitude')
#cbar.set_label('2" K Flux')
#cbar.set_label('$M_{star}$')
plt.tight_layout()

### plot out of range z ###
plt.errorbar(z[zlimsmask], I1minI2[zlimsmask],fmt='.',xerr=0.05, markersize=1, 
             xlolims=True, zorder=0, color='tab:gray', alpha=0.5)
#plt.errorbar(xz[xzlimsmask], xI1minI2[xzlimsmask],fmt='k+',xerr=0.05, 
#             xlolims=True, zorder=0)
plt.errorbar(varyz[varyzlimsmask], varyI1minI2[varyzlimsmask],fmt='C0.',
             xerr=0.05, xlolims=True, zorder=0)
#plt.errorbar(xvaryz[xvaryzlimsmask], xvaryI1minI2[xvaryzlimsmask],fmt='ro',
#             xerr=0.05, mfc='None', xlolims=True, zorder=0)
#plt.errorbar(noxvaryz[noxvaryzlimsmask], noxvaryI1minI2[noxvaryzlimsmask],
#             fmt='bo',xerr=0.05, mfc='None', xlolims=True, zorder=0)

### plot out of range low I1minI2 ###
plt.errorbar(z[Ilolimsmask], I1minI2[Ilolimsmask],fmt='.',yerr=0.03, 
             markersize=1, uplims=True, zorder=0, color='tab:gray', alpha=0.5)
#plt.errorbar(xz[xIlolimsmask], xI1minI2[xIlolimsmask],fmt='k+',yerr=0.03, 
#             uplims=True, zorder=0)
plt.errorbar(varyz[varyIlolimsmask], varyI1minI2[varyIlolimsmask], fmt='C0.',
             yerr=0.03, uplims=True, zorder=0)
#plt.errorbar(xvaryz[xvaryIlolimsmask], xvaryI1minI2[xvaryIlolimsmask], fmt='ro',
#             yerr=0.03, mfc='None', uplims=True, zorder=0)
#plt.errorbar(noxvaryz[noxvaryIlolimsmask], noxvaryI1minI2[noxvaryIlolimsmask],
#             fmt='bo', mfc='None', yerr=0.03, uplims=True, zorder=0)

#### plot out of range high I1minI2 ###
plt.errorbar(z[Iuplimsmask], I1minI2[Iuplimsmask],fmt='.', yerr=0.03, 
             markersize=1, lolims=True, zorder=0, color='tab:gray', alpha=0.5)
#plt.errorbar(xz[xIuplimsmask], xI1minI2[xIuplimsmask],fmt='.', yerr=0.03, 
#             markersize=1, lolims=True, zorder=0, color='tab:gray', alpha=0.5)
plt.errorbar(varyz[varyIuplimsmask], varyI1minI2[varyIuplimsmask], fmt='C0.',
             yerr=0.03, lolims=True, zorder=0)
#plt.errorbar(xvaryz[xvaryIuplimsmask], xvaryI1minI2[xvaryIuplimsmask], fmt='ro',
#             yerr=0.03, mfc='None', lolims=True, zorder=0)
#plt.errorbar(noxvaryz[noxvaryIuplimsmask], noxvaryI1minI2[noxvaryIuplimsmask],
#             fmt='bo', mfc='None', yerr=0.03, lolims=True)