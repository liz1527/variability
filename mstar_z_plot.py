#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 12:23:28 2019

Code to plot M_star vs z

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
import matplotlib.colors as colors
plt.close('all')


### Get fits ###
tbdata = fits.open('mag_flux_tables/K/mag_flux_table_best_extra_clean_no06.fits')[1].data
dr11 = fits.open('UDS_catalogues/DR11-2arcsec-June24-2018+plusXY_best.fits')[1].data
varydata = fits.open('variable_tables/no06_variables_chi30_2arcsec_DR11data_restframe.fits')[1].data
fullxray = Table.read('UDS_catalogues/DR11-2arcsec-June24-2018+plusXY_best_chandra.fits')
sternvary = fits.open('variable_tables/no06_variables_chi30_2arcsec_DR11data_SpUDSdata_IRAC_stern.fits')[1].data
stern = fits.open('UDS_catalogues/SpUDS_IRAC_catalogue_DR11data_stern.fits')[1].data
nosternvarydata = fits.open('variable_tables/no06_variables_chi30_2arcsec_DR11data_SpUDSdata_IRAC.fits')[1].data
deviant = fits.open('variable_tables/no06_variables_chi30_2arcsec_deviant_DR11data_restframe.fits')[1].data
sn = fits.open('variable_tables/no06_variables_chi30_2arcsec_DR11data_restframe_SN.fits')[1].data
dev07 = fits.open('variable_tables/no06_variables_chi30_2arcsec_DR11data_restframe_07B.fits')[1].data
#deviant = fits.open('variable_tables/no06_variables_chi30_2arcsec_DR11data_restframe_allspecz.fits')[1].data

### with x/nox split ###
#noxvarydata = fits.open('variable_tables/no06_variables_chi30_2arcsec_noXray_DR11data_restframe.fits')[1].data
#xvarydata = fits.open('variable_tables/no06_variables_chi30_2arcsec_Xray_DR11data_restframe.fits')[1].data
noxvarydata = fits.open('variable_tables/no06_variables_chi30_2arcsec_nochanXray_DR11data_restframe.fits')[1].data
xvarydata = fits.open('variable_tables/no06_variables_chi30_2arcsec_chandata_DR11data_restframe.fits')[1].data

noxsternvarydata = fits.open('variable_tables/no06_variables_chi30_2arcsec_DR11data_SpUDSdata_IRAC_noXray_stern.fits')[1].data
xsternvarydata = fits.open('variable_tables/no06_variables_chi30_2arcsec_DR11data_SpUDSdata_IRAC_Xray_stern.fits')[1].data

def ignore_zeros(m):
    m[m==0] = np.nan
    mask = ~np.isnan(m)
    m = m[mask]
    return m, mask

def prep_variables(tbdata):
    z = vari_funcs.get_z(tbdata) # get redshift
    m = tbdata['Mstar_z_p'] #get mass
    m, mask = ignore_zeros(m) #mask those with null mass
    z = z[mask] #apply mask to z array
    return z, m


z, m = prep_variables(dr11)
noxz, noxm = prep_variables(noxvarydata)
xz, xm = prep_variables(xvarydata)
varyz, varym = prep_variables(varydata)
allxz, allxm = prep_variables(fullxray)
allsternz, allsternm = prep_variables(stern)
sternz, sternm = prep_variables(sternvary)
nosternz, nosternm = prep_variables(nosternvarydata)
noxsternz, noxsternm = prep_variables(noxsternvarydata)
xsternz, xsternm = prep_variables(xsternvarydata)
devz, devm = prep_variables(deviant)
snz, snm = prep_variables(sn)
dev07z, dev07m = prep_variables(dev07)

print(str(len(noxz[noxz!=-1]))+' non-x-ray variables')
print(str(len(xz[xz!=-1]))+' x-ray variables')

### KS test ###
from scipy import stats
[Dx_nox, px_nox] = stats.ks_2samp(noxm, xm)
[Dallx_vary, pallx_vary] = stats.ks_2samp(varym, allxm)
[Ddev_vary, pdev_vary] = stats.ks_2samp(varym, devm)
[Dx_allx, px_allx] = stats.ks_2samp(xm, allxm)
[Dnox_allx, pnox_allx] = stats.ks_2samp(noxm, allxm)

[Dstern_nostern, pstern_nostern] = stats.ks_2samp(nosternm, sternm)
[Dallstern_vary, pallstern_vary] = stats.ks_2samp(varym, allsternm)
[Dstern_allstern, pstern_allstern] = stats.ks_2samp(sternm, allsternm)
[Dnostern_allstern, pnostern_allstern] = stats.ks_2samp(nosternm, allsternm)

#plt.figure()
#plt.hist([noxm, xm, devm, varym], bins=np.logspace(4.5,12), color=['b','r','y','g'], 
#         histtype='step', label=['Non-X-ray Variable','X-ray Variable','Deviant','all variables'])
#plt.xlabel('$M_{star}$')
#plt.ylabel('Number')
#plt.xscale('log')
#
#plt.figure(figsize=[9,7])
plt.figure(figsize=[10,7])
plt.plot(z, m, '.',markersize=1, color='tab:gray', alpha=0.2, label='Galaxy')
#plt.plot(allsternz, allsternm, 'm.', label='Stern AGN', alpha=0.5)
plt.plot(allxz, allxm, 'ks',markersize=5, label='Non-Variable X-ray AGN')
plt.plot(xz, xm, 'ro', label='Variable X-ray AGN')
plt.plot(noxz, noxm, 'bo', label='Variable Non-X-ray AGN')
#plt.plot(noxz, noxm, 'rs',mfc='None', label='Variable AGN')
#plt.plot(xz, xm, 'rs', mfc='None')
#plt.plot(devz, devm, 'yd', label='Deviant in 07B')
#plt.plot(sternz, sternm, 'yd',markersize=10, mfc='none', label='Stern Variable')

#plt.hlines(2e9,-0.1,4.5,linestyle='dashed')

plt.yscale('log')
plt.xlim(xmin=-0.1, xmax=4.5)
plt.ylim(ymin=2e6, ymax=5e11)
plt.legend(loc='lower right')
plt.xlabel('Redshift')
plt.ylabel('Stellar Mass ($M_{\odot}$)')
plt.tight_layout()

### Plot with stern details
#plt.figure(figsize=[7,5])
plt.figure(figsize=[10,7])
font = {'family' : 'DejaVu Sans',
        'weight' : 'normal',
        'size'   : 18}
plt.rc('font', **font)
plt.plot(z, m, '.',markersize=1, color='tab:gray', alpha=0.2, label='Galaxy')
plt.plot(allsternz, allsternm, 'k.', label='Non-Variable Stern AGN', alpha=0.5)
#plt.plot(allxz, allxm, 'k+', label='X-ray Non-Variable')
#plt.plot(noxsternz, noxsternm, 'bo', label='Non-X-ray Variable')
#plt.plot(xsternz, xsternm, 'ro', label='X-ray Variable')
plt.plot(nosternz, nosternm, 'go', label='Variable Non-Stern AGN')
#plt.plot(xz, xm, 'go')#, label='X-ray Variable')
plt.plot(sternz, sternm, 'mo', label='Variable Stern AGN')

#plt.hlines(2e9,-0.1,4.5,linestyle='dashed')

plt.yscale('log')
plt.xlim(xmin=-0.1, xmax=4.5)
plt.ylim(ymin=2e6, ymax=5e11)
plt.legend(loc='lower right')
plt.xlabel('Redshift')
plt.ylabel('Stellar Mass ($M_{\odot}$)')
plt.tight_layout()

#%% Plot with deviant details ###
#### remove deviants from main sample ###
#mask = np.isin(noxvarydata['NUMBER_1'], deviant['NUMBER_1'])
#noxvarydata = noxvarydata[~mask]
#mask = np.isin(xvarydata['NUMBER_1'], deviant['NUMBER_1'])
#xvarydata = xvarydata[~mask]
#
#noxz, noxm = prep_variables(noxvarydata)
#xz, xm = prep_variables(xvarydata)
#
#### remove 07Bs for talk slide ###
#mask = np.isin(deviant['NUMBER_1'], dev07['NUMBER_1'])
#deviant = deviant[~mask]
#devz, devm = prep_variables(deviant)
#
##plt.figure(figsize=[9,7])
#plt.figure(figsize=[10,7])
#plt.plot(z, m, '.',markersize=1, color='tab:gray', alpha=0.25, label='Galaxy')
##plt.plot(allxz, allxm, 'ks', label='X-ray AGN')
#plt.plot(xz, xm, 'ro', label='Variable X-ray AGN')
#plt.plot(noxz, noxm, 'bo', label='Variable Non-X-ray AGN')
#plt.plot(devz, devm, 'gs', markersize=6, label='Deviant')
#plt.plot(snz, snm, 'y*', markersize=12, label='Potential SN')
##plt.plot(dev07z, dev07m, 'gd', markersize=8, label='Deviant in 07B')
#
##plt.hlines(2e9,-0.1,4.5,linestyle='dashed')
#
#plt.yscale('log')
#plt.xlim(xmin=-0.1, xmax=4.5)
#plt.ylim(ymin=2e5, ymax=5e11)
#plt.legend(loc='lower right')
#plt.xlabel('Redshift')
#plt.ylabel('Stellar Mass ($M_{\odot}$)')
##plt.xlabel('z')
##plt.ylabel('$M_{star}$')
#plt.tight_layout()
###%% Plot with flux colours ###
##### get fluxes from mag-flux ###
###
###def get_mean_flux(tbdata):
###    flux = vari_funcs.flux4_stacks(tbdata)
###    meanflux = np.nanmean(flux, axis=1)
###    return meanflux
###
###def get_jansky_flux(tbdata):
###    meanmag = tbdata['KMAG_20']
###    meanflux = 10**(23-((meanmag+48.6)/2.5))
###    return meanflux
###    
###meanflux = get_mean_flux(tbdata)
###meannoxvary = get_mean_flux(noxvarydata)
###meanxvary = get_mean_flux(xvarydata)
###meanvary = get_mean_flux(varydata)
###
###### mask flux arrays and make flux colour ###
###meannoxvary = meannoxvary[noxmask]
###meanxvary = meanxvary[xmask]
###
###### find max and min ###
###cmax = np.nanmax([np.nanmax(meannoxvary), np.nanmax(meanxvary)])
###cmin = np.nanmin([np.nanmin(meannoxvary), np.nanmin(meanxvary)])
###
###plt.figure(figsize=[10,7])
###plt.plot(z, m, '.',markersize=1, color='tab:gray', alpha=0.25, label='UDS Galaxy')
###plt.plot(allxz, allxm, 'k+', label='X-ray Non-Variable')
####plt.scatter(noxz, noxm, marker='o', c=meannoxvary, label='Non-X-ray Variable', 
####            norm=colors.LogNorm(vmin=cmin, vmax=cmax), zorder=3)
####plt.scatter(xz, xm, marker='o', c=meanxvary, label='X-ray Variable', 
####         norm=colors.LogNorm(vmin=cmin, vmax=cmax), zorder=3)
####cbar=plt.colorbar()
####cbar.set_label('Mean 2" Flux')
###plt.yscale('log')
###plt.xlim(xmin=-0.1, xmax=4.5)
###plt.ylim(ymin=1e4, ymax=3e12)
###plt.legend(loc='lower right')
###plt.xlabel('z')
###plt.ylabel('$M_{star}$')
###plt.tight_layout()
##
