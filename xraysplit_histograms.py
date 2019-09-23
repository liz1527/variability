#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 28 16:25:12 2019

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
plt.close('all')


### Get fits ###
tbdata = fits.open('mag_flux_tables/mag_flux_table_best_extra_clean_no06.fits')[1].data
dr11 = fits.open('UDS_catalogues/DR11-2arcsec-June24-2018+plusXY_best.fits')[1].data
varydata = fits.open('variable_tables/no06_variables_chi30_2arcsec_DR11data_restframe.fits')[1].data
fullxray = Table.read('UDS_catalogues/DR11-2arcsec-June24-2018+plusXY_best_chandra.fits')
deviant = fits.open('variable_tables/no06_variables_chi30_2arcsec_deviant_DR11data_restframe.fits')[1].data

### with x/nox split ###
noxvarydata = fits.open('variable_tables/no06_variables_chi30_2arcsec_noXray_DR11data_restframe.fits')[1].data
xvarydata = fits.open('variable_tables/no06_variables_chi30_2arcsec_Xray_DR11data_restframe.fits')[1].data
#noxvarydata = fits.open('variable_tables/no06_variables_chi30_2arcsec_nochanXray_DR11data_restframe.fits')[1].data
#xvarydata = fits.open('variable_tables/no06_variables_chi30_2arcsec_chandata_DR11data_restframe.fits')[1].data

#%% restrict to just chandra ###
tbdata = vari_funcs.chandra_only(tbdata)
varydata = vari_funcs.chandra_only(varydata)
noxvarydata = vari_funcs.chandra_only(noxvarydata)
xvarydata = vari_funcs.chandra_only(xvarydata)
deviant = vari_funcs.chandra_only(deviant)

#%% histograms of z ###
z = vari_funcs.get_z(dr11)
noxz = vari_funcs.get_z(noxvarydata)
xz = vari_funcs.get_z(xvarydata)
allxz = vari_funcs.get_z(fullxray)

plt.figure()
#plt.hist([noxz, xz], color=['b','r'],histtype='step', 
#         label=['Non X-ray',r'X-ray'])#,
##         normed=True)
plt.hist([xz, allxz], bins=50, color=['r','k'], 
         histtype='step', label=['X-ray Variable','X-ray Non-Variable'])#,
#         normed=True)
plt.legend()
plt.xlabel('z')
plt.ylabel('Number')

#%% histograms of flux ###


### get fluxes from mag-flux ###

def get_mean_flux(tbdata):
    flux = vari_funcs.flux4_stacks(tbdata)
    meanflux = np.nanmean(flux, axis=1)
    return meanflux

def get_jansky_flux(tbdata):
    meanmag = tbdata['KMAG_20']
    meanflux = 10**(23-((meanmag+48.6)/2.5))
    return meanflux
    
meanflux = get_mean_flux(tbdata)
meannoxvary = get_mean_flux(noxvarydata)
meanxvary = get_mean_flux(xvarydata)
meanvary = get_mean_flux(varydata)

### get fluxs from DR11data ###
#meanflux = get_jansky_flux(dr11)
#dr11vary = get_jansky_flux(varydata)

### plot ###
#bins = np.logspace(-8,-3,50)
#bins = np.logspace(0,7,50)    
bins = np.arange(13,25,0.2)
bins = np.append(bins, [25])
bins = 10**((30-bins)/2.5)
bins = np.flip(bins, axis=0)

plt.figure()
plt.hist(meanflux,  bins, histtype='step', normed='True')
plt.hist(meannoxvary,  bins, histtype='step', color='b', normed='True')
plt.xscale('log')
plt.xlabel('2" K band flux')
plt.ylabel('Number')
plt.tight_layout()

plt.figure()
#plt.hist(meanvary,  bins, histtype='step', color='k')#, normed='True')
plt.hist(meannoxvary,  bins, histtype='step', color='b', label='Non X-ray')#, normed='True')
plt.hist(meanxvary,  bins, histtype='step', color='r', label='X-ray')#, normed='True')
plt.xscale('log')
plt.xlabel('2" K band flux')
plt.ylabel('Number')
plt.legend()
plt.tight_layout()

#%% Histograms of luminosity ###

L = varydata['M_K_z_p']
xL = xvarydata['M_K_z_p']
noxL = noxvarydata['M_K_z_p']
allxL = fullxray['M_K_z_p']
devL = deviant['M_K_z_p']

def remove_99s(L):
    ### remove any that have val of 99 ###
    L[L == 99] = np.nan # remove those with null values
    L[L == 999] = np.nan # remove those with null values
    mask = ~np.isnan(L)
    return L[mask]

L = remove_99s(L)
xL = remove_99s(xL)
noxL = remove_99s(noxL)
allxL = remove_99s(allxL)
devL = remove_99s(devL)

#maxL = np.nanmax([np.nanmax(noxL), np.nanmax(xL), np.nanmax(allxL)])
#minL = np.nanmin([np.nanmin(noxL), np.nanmin(xL), np.nanmin(allxL)])
#bins = np.linspace(minL,maxL,50)

plt.figure(figsize=[7,5])
#plt.hist([noxL, xL], color=['b','r'], bins=50,histtype='step', 
#         label=['Non-X-ray Variable','X-ray Variable'],
#         normed=True)
#plt.hist([xL, allxL], bins=50, color=['r','k'], 
#         histtype='step', label=['X-ray Variable','X-ray Non-Variable'])#,
#         normed=True)
_, bins, _ = plt.hist([noxL, xL, allxL], bins=50, color=['b','r','k'], #linestyle=['dashed','dashed','dashed'],
         histtype='step', label=['Variable Non-X-ray AGN','Variable X-ray AGN','X-ray AGN'],
         normed=True)
plt.figure(figsize=[7,5])
plt.hist(noxL, bins, color='b', linestyle='dashdot',
         histtype='step', label='Variable Non-X-ray AGN',
         normed=True, linewidth=1.5)
plt.hist(xL, bins, color='r', linestyle='solid',
         histtype='step', label='Variable X-ray AGN',
         normed=True, linewidth=1.5)
plt.hist(allxL, bins, color='k', linestyle='dashed',
         histtype='step', label='X-ray AGN',
         normed=True, linewidth=1.5)
#plt.hist([noxL, xL, devL], bins=50, color=['b','r','y'], 
#         histtype='step', label=['Variable Non-X-ray','Variable X-ray','deviant'])#,
plt.legend(loc='upper left')
plt.xlabel('K-band Absolute Magnitude')
plt.xlim(xmax=-15)
plt.ylabel('Normalised Frequency')
plt.gca().invert_xaxis()
plt.tight_layout()


#%% Histograms of stellar mass ###
m = dr11['Mstar_z_p']
noxm = noxvarydata['Mstar_z_p']
xm = xvarydata['Mstar_z_p']
allxm = fullxray['Mstar_z_p']

### ignore 0s ###
def ignore_zeros(m):
    m[m==0] = np.nan
    mask = ~np.isnan(m)
    m = m[mask]
    return m, mask

m, mask = ignore_zeros(m)
noxm, noxmask = ignore_zeros(noxm)
xm, xmask = ignore_zeros(xm)
allxm, allxmask = ignore_zeros(allxm)

### plot hist for info ###
plt.figure()
#plt.hist([noxm, xm], bins=np.logspace(4.5,12), color=['b','r'], histtype='step', label=['Non X-ray','X-ray'])
#plt.hist([noxm, xm, allxm], bins=np.logspace(4.5,12), color=['b','r','k'], 
#         histtype='step', label=['Non X-ray Variable','X-ray Variable','X-ray Non-Variable'])
plt.hist([xm, allxm], bins=np.logspace(4.5,12), color=['r','k'], 
         histtype='step', label=['X-ray Variable','X-ray Non-Variable'])
#         , normed=True)
plt.xscale('log')
plt.xlabel(r'$M_{star}$')
plt.ylabel('Number')
plt.legend(loc='upper left')
plt.tight_layout()

#%% Stellar mass vs z ###
import matplotlib.colors as colors
spudsvary = fits.open('variable_tables/no06_variables_chi30_2arcsec_DR11data_SpUDSdata_IRAC_stern.fits')[1].data
spudsm = spudsvary['Mstar_z_p']
spudsm, spudsmask = ignore_zeros(spudsm)
spudsz = vari_funcs.get_z(spudsvary)

### mask z arrays ###
z = z[mask]
allxz = allxz[allxmask]
xz = xz[xmask]
noxz = noxz[noxmask]
spudsz = spudsz[spudsmask]

#plt.figure(figsize=[9,7])
plt.figure(figsize=[10,7])
plt.plot(z, m, '.',markersize=1, color='tab:gray', alpha=0.25, label='UDS Galaxy')
plt.plot(allxz, allxm, 'k+', label='X-ray Non-Variable')
plt.plot(noxz, noxm, 'bo', label='Non-X-ray Variable')
plt.plot(xz, xm, 'ro', label='X-ray Variable')
plt.plot(spudsz, spudsm, 'md', mfc='none', label='SpUDS Variable')

### mask flux arrays and make flux colour ###
meannoxvary = meannoxvary[noxmask]
meanxvary = meanxvary[xmask]

### find max and min ###
cmax = np.nanmax([np.nanmax(meannoxvary), np.nanmax(meanxvary)])
cmin = np.nanmin([np.nanmin(meannoxvary), np.nanmin(meanxvary)])

#plt.scatter(noxz, noxm, marker='o', c=meannoxvary, label='Non-X-ray Variable', 
#            norm=colors.LogNorm(vmin=cmin, vmax=cmax), zorder=3)
#plt.scatter(xz, xm, marker='o', c=meanxvary, label='X-ray Variable', 
#         norm=colors.LogNorm(vmin=cmin, vmax=cmax), zorder=3)
#cbar=plt.colorbar()
#cbar.set_label('Mean 2" Flux')
plt.yscale('log')
plt.xlim(xmin=-0.1, xmax=4.5)
plt.ylim(ymin=1e4, ymax=3e12)
plt.legend(loc='lower right')
plt.xlabel('z')
plt.ylabel('$M_{star}$')
plt.tight_layout()

#%% chi squared distribution ###
sigtb = Table.read('sigma_tables/quad_epoch_sigma_table_extra_clean_no06_2arcsec.fits')

xflux = vari_funcs.flux4_stacks(xvarydata)
xflux, xfluxerr, xvarydata = vari_funcs.create_quad_error_array(sigtb, xvarydata, aper=4)

noxflux = vari_funcs.flux4_stacks(noxvarydata)
noxflux, noxfluxerr, noxvarydata = vari_funcs.create_quad_error_array(sigtb, noxvarydata, aper=4)

xchi = vari_funcs.my_chisquare_err(xflux, xfluxerr)
noxchi = vari_funcs.my_chisquare_err(noxflux, noxfluxerr)

bins = np.logspace(np.log10(3e1),np.log10(5e3))
plt.figure()
plt.hist([noxchi, xchi], bins, color=['b','r'],histtype='step', 
         label=['Variable Non-X-ray','Variable X-ray'])#,
#         normed=True)
plt.xscale('log')
plt.ylabel('Number')
plt.xlabel('$\chi^{2}$')
plt.tight_layout()
