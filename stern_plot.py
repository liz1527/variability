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
tbdata = fits.open('UDS_catalogues/SpUDS_IRAC_catalogue_DR11data.fits')[1].data
xdata = fits.open('UDS_catalogues/SpUDS_IRAC_catalogue_DR11data_chandra.fits')[1].data
#tbdata = fits.open('UDS_catalogues/SpUDS_IRAC_catalogue.fits')[1].data
varydata = fits.open('variable_tables/no06_variables_chi30_2arcsec_DR11data_SpUDSdata_IRAC.fits')[1].data
xvarydata = fits.open('variable_tables/no06_variables_chi30_2arcsec_Xray_DR11data_SpUDSdata_IRAC.fits')[1].data
noxvarydata = fits.open('variable_tables/no06_variables_chi30_2arcsec_noXray_DR11data_SpUDSdata_IRAC.fits')[1].data
#%%
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
#    I1[I1 == -99] = np.nan # remove those with null values
    mask1 = ~np.isnan(I1)
    I2[I2 == 99] = np.nan # remove those with null values
#    I2[I2 == -99] = np.nan # remove those with null values
    mask2 = ~np.isnan(I2)
#    mask = mask1*mask2.astype('bool')
    I3[I3 == 99] = np.nan # remove those with null values
#    I3[I3 == -99] = np.nan # remove those with null values
    mask3 = ~np.isnan(I3)
    I4[I4 == 99] = np.nan # remove those with null values
#    I4[I4 == -99] = np.nan # remove those with null values
    mask4 = ~np.isnan(I4)
    mask = mask1*mask2*mask3*mask4.astype('bool')
    I1 = I1[mask]
    I2 = I2[mask]
    I3 = I3[mask]
    I4 = I4[mask]
    tbdata = tbdata[mask]
    
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
    I34lolimsmask = I3minI4 < -0.5
#    Ilims = I1minI2[mask]
    
    return I1minI2, I12uplimsmask, I12lolimsmask, I3minI4, I34uplimsmask, \
            I34lolimsmask, tbdata
### Limit to Chandra region ###
#noxvarydata = vari_funcs.chandra_only(noxvarydata)
#xvarydata = vari_funcs.chandra_only(xvarydata)

### remove those with very high phot z ###
#xvarydata = xvarydata[xvarydata['z_p']<6]
#noxvarydata = noxvarydata[noxvarydata['z_p']<6]

#xvarydataup = vari_funcs.flux_split(xvarydata, 'upper')
#noxvarydataup = vari_funcs.flux_split(noxvarydata, 'upper')
#
#xvarydatalow = vari_funcs.flux_split(xvarydata, 'lower')
#noxvarydatalow = vari_funcs.flux_split(noxvarydata, 'lower')

I1minI2, I12uplimsmask, I12lolimsmask, I3minI4, I34uplimsmask, I34lolimsmask, \
tbdata = get_data(tbdata)

xI1minI2, xI12uplimsmask, xI12lolimsmask, xI3minI4, xI34uplimsmask, \
xI34lolimsmask, xdata = get_data(xdata)

varyI1minI2, varyI12uplimsmask, varyI12lolimsmask, varyI3minI4, \
varyI34uplimsmask, varyI34lolimsmask, varydata = get_data(varydata)


xvaryI1minI2, xvaryI12uplimsmask, xvaryI12lolimsmask, xvaryI3minI4, \
xvaryI34uplimsmask, xvaryI34lolimsmask, xvarydata = get_data(xvarydata)

noxvaryI1minI2, noxvaryI12uplimsmask, noxvaryI12lolimsmask, noxvaryI3minI4, \
noxvaryI34uplimsmask, noxvaryI34lolimsmask, noxvarydata = get_data(noxvarydata)

#xvaryI1minI2up, xvaryI12uplimsmaskup, xvaryI12lolimsmaskup, xvaryI3minI4up, \
#xvaryI34uplimsmaskup, xvaryI34lolimsmaskup, xvarydataup = get_data(xvarydataup)
#
#noxvaryI1minI2up, noxvaryI12uplimsmaskup, noxvaryI12lolimsmaskup, \
#noxvaryI3minI4up, noxvaryI34uplimsmaskup, noxvaryI34lolimsmaskup, \
#noxvarydataup = get_data(noxvarydataup)
#
#xvaryI1minI2low, xvaryI12uplimsmasklow, xvaryI12lolimsmasklow, \
#xvaryI3minI4low, xvaryI34uplimsmasklow, xvaryI34lolimsmasklow, \
#xvarydatalow = get_data(xvarydatalow)
#
#noxvaryI1minI2low, noxvaryI12uplimsmasklow, noxvaryI12lolimsmasklow, \
#noxvaryI3minI4low, noxvaryI34uplimsmasklow, noxvaryI34lolimsmasklow, \
#noxvarydatalow = get_data(noxvarydatalow)


### Plot result ###
#plt.figure(figsize=[11,7])
plt.figure(figsize=[10,7])

#z = xvarydata['z_spec']
#z[z==-1] = xvarydata['z_p'][z==-1]

#varyL = xvarydata['M_K_z_p']

### flux ###
#varyL = vari_funcs.flux4_stacks(varydata)
#varyL = np.nanmean(varyL, axis=1)
#xvaryL = vari_funcs.flux4_stacks(xvarydata)
#xvaryL = np.nanmean(xvaryL, axis=1)
#noxvaryL = vari_funcs.flux4_stacks(noxvarydata)
#noxvaryL = np.nanmean(noxvaryL, axis=1)
#colmap='cool'

### M_star ###
#varyL = varydata['Mstar_z_p']
#xvaryL = xvarydata['Mstar_z_p']
#noxvaryL = noxvarydata['Mstar_z_p']
#
#def ignore_zeros(m, z, I1minI2):
#    m[m==0] = np.nan
#    mask = ~np.isnan(m)
#    z = z[mask]
#    I1minI2 = I1minI2[mask]
#    m = m[mask]
#    return m, z, I1minI2
#
#varyL, varyI3minI4, varyI1minI2,  = ignore_zeros(varyL, varyI3minI4, 
#                                                           varyI1minI2)
#xvaryL, xvaryI3minI4, xvaryI1minI2,  = ignore_zeros(xvaryL, xvaryI3minI4, 
#                                                           xvaryI1minI2)
#noxvaryL, noxvaryI3minI4, noxvaryI1minI2,  = ignore_zeros(noxvaryL, noxvaryI3minI4, 
#                                                           noxvaryI1minI2)

##colmap='summer'
#### find max and min ###
#cmax = np.nanmax(varyL)
#cmin = np.nanmin(varyL)
#
##zcolour = vari_funcs.get_z(tbdata)
#zcolour = tbdata['z_spec']
#zcolour[zcolour==-1] = np.nan
#clim = np.nanmax(zcolour)

### set values out of range to their plottable values ###
#I1minI2[I12lolimsmask]=-0.25
#I1minI2[I12uplimsmask]= 1.47
#I3minI4[I34lolimsmask]=-0.25
#I3minI4[I34uplimsmask]= 3.45

#### plot with lum as colour ####
plt.plot(I3minI4, I1minI2,'.', color='tab:grey', markersize=1, label='Galaxy', 
         zorder=0, alpha=0.5)
plt.plot(xI3minI4, xI1minI2,'ks', markersize=5,  label='X-ray AGN',zorder=0)
#plt.scatter(I3minI4, I1minI2, c=zcolour, marker='.',vmax=clim, label='Galaxy')#, zorder=0, alpha=0.5)
#plt.scatter(varyI3minI4, varyI1minI2, c=varyL,
#            norm=colors.LogNorm(vmin=cmin, vmax=cmax), marker='o')
#plt.plot(noxvaryI3minI4, noxvaryI1minI2,'bo', mfc='None', label='Non X-ray Variable')
#plt.plot(xvaryI3minI4, xvaryI1minI2,'ro', mfc='None', label='X-ray Variable')

#plt.scatter(xvaryI3minI4, xvaryI1minI2, c=xvaryL,
#            norm=colors.LogNorm(vmin=cmin, vmax=cmax), marker='o')
#plt.scatter(noxvaryI3minI4, noxvaryI1minI2, c=noxvaryL,
#            norm=colors.LogNorm(vmin=cmin, vmax=cmax), marker='x')


### plot those within range
#plt.plot(I3minI4, I1minI2,'k.', markersize=1, label='SpUDS source')
#plt.plot(varyI3minI4, varyI1minI2,'kd', label='Variable SpUDS source')
plt.plot(noxvaryI3minI4, noxvaryI1minI2,'bo', label='Variable Non-X-ray AGN')
plt.plot(xvaryI3minI4, xvaryI1minI2,'ro', label='Variable X-ray AGN')

#plt.plot(noxvaryI3minI4up, noxvaryI1minI2up,'go', label='High flux, not X-ray')
#plt.plot(xvaryI3minI4up, xvaryI1minI2up,'gd',label='High flux, X-ray')
#
#plt.plot(noxvaryI3minI4low, noxvaryI1minI2low,'mo',label='Low flux, not X-ray')
#plt.plot(xvaryI3minI4low, xvaryI1minI2low,'md',label='Low flux, X-ray')
#plt.xlim(xmin=-1.2, xmax=2.8)         
#plt.ylim(ymin=-0.7, ymax=1)
plt.xlim(xmin=-0.5,xmax=3.5)
plt.ylim(ymin=-0.3, ymax=1.5)
plt.xlabel('5.6 - 8.0')
plt.ylabel('3.6 - 4.5')
plt.legend(loc='upper right')
#cbar = plt.colorbar()
#cbar.set_label('2" K Flux')
#cbar.set_label('$M_{star}$')
plt.tight_layout()


### find variables that satisfy ###
def stern_selection(I3minI4, I1minI2, tbdata):
    mask1 = I3minI4 > 0.6
    mask2 = I1minI2 > 0.2*(I3minI4) + 0.18
    mask3 = I1minI2 > 2.5*(I3minI4) - 3.5
    mask = mask1*mask2*mask3.astype(bool)
    data = np.copy(tbdata)
    sterndata = data[mask]
    return sterndata

stern = stern_selection(I3minI4, I1minI2, tbdata)
sternvary = stern_selection(varyI3minI4, varyI1minI2, varydata)
xsternvary = stern_selection(xvaryI3minI4, xvaryI1minI2, xvarydata)
xstern = stern_selection(xI3minI4, xI1minI2, xdata)
noxsternvary = stern_selection(noxvaryI3minI4, noxvaryI1minI2, noxvarydata)

print('Total stern: '+str(len(stern)))
print('Total stern x-ray: '+str(len(xstern)))
print('Total stern variables: '+str(len(sternvary)))
print('Total stern x-ray variables: '+str(len(xsternvary)))
print('Total stern not x-ray variables: '+str(len(noxsternvary)))

#
#### plot out of range low I1minI2 ###
#I1minI2[I12lolimsmask]=-0.25
#xI1minI2[xI12lolimsmask]=-0.25
##varyI1minI2[varyI12lolimsmask]=-0.25
#plt.errorbar(I3minI4[I12lolimsmask], I1minI2[I12lolimsmask],fmt='.', 
#             color='tab:grey',yerr=0.03, markersize=1, uplims=True, alpha=0.5)
#plt.errorbar(xI3minI4[xI12lolimsmask], xI1minI2[xI12lolimsmask],fmt='k+', 
#             yerr=0.03, uplims=True)
##plt.errorbar(varyI3minI4[varyI12lolimsmask], varyI1minI2[varyI12lolimsmask],
##             fmt='bo',yerr=0.03, uplims=True)
##plt.errorbar(xvaryI3minI4[xvaryI12lolimsmask], xvaryI1minI2[xvaryI12lolimsmask],
##             fmt='bo',yerr=0.03, uplims=True)
#
#
#### plot out of range high I1minI2 ###
#I1minI2[I12uplimsmask]= 1.45
#xI1minI2[xI12uplimsmask]= 1.45
##varyI1minI2[varyI12uplimsmask]= 1.45
#plt.errorbar(I3minI4[I12uplimsmask], I1minI2[I12uplimsmask],fmt='.', 
#             color='tab:grey',yerr=0.03, markersize=1, lolims=True, alpha=0.5)
#plt.errorbar(xI3minI4[xI12uplimsmask], xI1minI2[xI12uplimsmask],fmt='k+', 
#             yerr=0.03, lolims=True)
#plt.errorbar(varyI3minI4[varyI12uplimsmask], varyI1minI2[varyI12uplimsmask],
#             fmt='bo',yerr=0.03, lolims=True)
#
#
#### plot out of range low I3minI4 ###
#I3minI4[I34lolimsmask]=-0.45
#xI3minI4[xI34lolimsmask]=-0.45
##varyI1minI2[varyI12lolimsmask]=-0.25
#plt.errorbar(I3minI4[I34lolimsmask], I1minI2[I34lolimsmask],fmt='.', 
#             color='tab:grey',xerr=0.03, markersize=1, xuplims=True, alpha=0.5)
#plt.errorbar(xI3minI4[xI34lolimsmask], xI1minI2[xI34lolimsmask],fmt='k+', 
#             xerr=0.03, xuplims=True)
#plt.errorbar(varyI3minI4[varyI12lolimsmask], varyI1minI2[varyI12lolimsmask],
#             fmt='bo',yerr=0.03, uplims=True)
#
#
#### plot out of range high I3minI4 ###
#I3minI4[I34uplimsmask]= 3.45
#xI3minI4[xI34uplimsmask]= 3.45
##varyI1minI2[varyI34uplimsmask]= 3.45
#plt.errorbar(I3minI4[I34uplimsmask], I1minI2[I34uplimsmask],fmt='.', 
#             color='tab:grey',xerr=0.03, markersize=1, xlolims=True, alpha=0.5)
#plt.errorbar(xI3minI4[xI34uplimsmask], xI1minI2[xI34uplimsmask],fmt='k+', 
#             xerr=0.03, xlolims=True)
#plt.errorbar(varyI3minI4[varyI34uplimsmask], varyI1minI2[varyI34uplimsmask],
#            fmt='bo',yerr=0.03, lolims=True)

### Plot wedge ###
#plt.vlines(-0.029, -0.7, 1,linestyles='dashed')
#x = np.linspace(-0.029, 1.5)
#y1 = 0.2*(x-0.629) + 0.18
#y2 = 2.5*(x-0.629) - 3.5
plt.vlines(0.6, 0.3, 1.5, 'g',linestyles='dashed', zorder=3)
x1 = np.linspace(0.6, 1.6)
x2 = np.linspace(1.6, 2)
y1 = 0.2*(x1) + 0.18
y2 = 2.5*(x2) - 3.5
plt.plot(x1,y1,'g--')
plt.plot(x2,y2,'g--')
plt.tight_layout()

#sterntb = Table(noxsternvary)
#sterntb.write('variable_tables/no06_variables_chi30_2arcsec_DR11data_SpUDSdata_IRAC_noXray_stern.fits')

#sterntb = Table(stern)
#sterntb.write('UDS_catalogues/SpUDS_IRAC_catalogue_DR11data_stern.fits')
#varyI1minI2, varyI12uplimsmask, varyI12lolimsmask, varyI3minI4, \
#varyI34uplimsmask, varyI34lolimsmask, varydata = get_data(sternvary)
##
#xvaryI1minI2, xvaryI12uplimsmask, xvaryI12lolimsmask, xvaryI3minI4, \
#xvaryI34uplimsmask, xvaryI34lolimsmask, xvarydata = get_data(xsternvary)
##
#noxvaryI1minI2, noxvaryI12uplimsmask, noxvaryI12lolimsmask, noxvaryI3minI4, \
#noxvaryI34uplimsmask, noxvaryI34lolimsmask, noxvarydata = get_data(noxsternvary)
##    
#plt.plot(varyI3minI4, varyI1minI2,'mo',mfc='None', label='Variable SpUDS source')
#plt.plot(noxvaryI3minI4, noxvaryI1minI2,'mo',mfc='None', label='Variable source')
#plt.plot(xvaryI3minI4, xvaryI1minI2,'mo',mfc='None', label='Variable X-ray source')