#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 14:10:37 2019

Code to find alphs-kx distrutions of variable sources in chandra region

@author: ppxee
"""

### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
from astropy.table import Table #for handling tables
import numpy as np #for handling arrays
import vari_funcs #my module to help run code neatly
import matplotlib.colors as colors
# Imports for luminosity distance
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
# Set up cosmology
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
import scipy.integrate as integrate
from astropy import constants
plt.close('all')

### Get fits ###
#dr11 = fits.open('UDS_catalogues/DR11-2arcsec-June24-2018+plusXY_best.fits')[1].data
fullxray = fits.open('UDS_catalogues/DR11-2arcsec-June24-2018+plusXY_chandra_novarys.fits')[1].data

### with x/nox split ###
#noxvarydata = fits.open('variable_tables/no06_variables_chi30_2arcsec_noXray_DR11data_restframe.fits')[1].data
#xvarydata = fits.open('variable_tables/no06_variables_chi30_2arcsec_Xray_DR11data_restframe.fits')[1].data
noxvarydata = fits.open('variable_tables/no06_variables_chi30_2arcsec_nochanXray_DR11data_restframe.fits')[1].data
xvarydata = fits.open('variable_tables/no06_variables_chi30_2arcsec_chandata_DR11data_restframe.fits')[1].data


allx_L_k = vari_funcs.xray_funcs.get_L_k(fullxray)
nox_L_k = vari_funcs.xray_funcs.get_L_k(noxvarydata)
x_L_k = vari_funcs.xray_funcs.get_L_k(xvarydata)

##_, bins, _ = plt.hist([nox_L_k, x_L_k, allx_L_k], bins=50, color=['b','r','k'], #linestyle=['dashed','dashed','dashed'],
##         histtype='step', label=['Variable Non-X-ray AGN','Variable X-ray AGN','X-ray AGN'],
##         normed=True)
#bins=np.logspace(16,np.log10(1.2e25))
#plt.figure(figsize=[7,5])
#plt.hist(nox_L_k, bins, color='b', linestyle='dashdot',
#         histtype='step', label='Variable Non-X-ray AGN')#,
##         normed=True, linewidth=1.5)
#plt.hist(x_L_k, bins, color='r', linestyle='solid',
#         histtype='step', label='Variable X-ray AGN')#,
##         normed=True, linewidth=1.5)
#plt.hist(allx_L_k, bins, color='k', linestyle='dashed',
#         histtype='step', label='X-ray AGN')#,
##         normed=True, linewidth=1.5)
#plt.xscale('log')

### Need to calculate the monochromatic luminosity for 2keV ###
### First need to assume power law for flux density and find constant ###
### Then sub this equation for flux density into luminosity density eq ###


allx_L_2, allx_F_2, allx_flux, allx_z = vari_funcs.xray_funcs.get_xray_L_2(fullxray)
nox_L_2, nox_F_2, nox_flux, nox_z = vari_funcs.xray_funcs.get_xray_L_2(noxvarydata, Xray=False)
x_L_2, x_F_2, x_flux, x_z = vari_funcs.xray_funcs.get_xray_L_2(xvarydata)

#print('Start flux = '+str(allx_flux[10]))#+' erg/cm**2/s')
#print('z = '+str(allx_z[10]))
#print('Flux density at 2keV = '+str(allx_F_2[10]))#+' erg/cm**2/s/keV')
print('Luminosity density at 2keV = '+str(allx_L_2[10]))#+' erg/s/eV')
print('Luminosity density at K = '+str(allx_L_k[10]))#+' erg/s/eV')


allx_alpha_kx = vari_funcs.xray_funcs.calc_alpha_kx(allx_L_k, allx_L_2)
nox_alpha_kx = vari_funcs.xray_funcs.calc_alpha_kx(nox_L_k, nox_L_2)
x_alpha_kx = vari_funcs.xray_funcs.calc_alpha_kx(x_L_k, x_L_2)

plt.figure()
upplims = nox_L_2.value/3
#plt.plot(allx_L_k, allx_L_2,'ko')
#plt.plot(nox_L_k, nox_L_2,'bo')
plt.errorbar(nox_L_k.value, nox_L_2.value,yerr=upplims, fmt='bo', #mfc='None', #markersize=10,
             uplims=True, label='Variable Non-X-ray AGN')
plt.plot(x_L_k, x_L_2,'ro', label='Variable X-ray AGN')
plt.xlim(xmin=1e24, xmax = 1e32)
plt.ylim(ymin=1e20, ymax = 4e27)
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$L_{K}$')
plt.ylabel(r'$L_{2 keV}$')
plt.legend()
plt.tight_layout()

_, bins, _ = plt.hist([nox_alpha_kx, x_alpha_kx, allx_alpha_kx], bins=50)#, color=['b','r','k'], #linestyle=['dashed','dashed','dashed'],
#         histtype='step', label=['Variable Non-X-ray AGN','Variable X-ray AGN','X-ray AGN']),
#         normed=True)
f, (a0, a1, a2, a3) = plt.subplots(4, 1, gridspec_kw={'height_ratios': [10, 1, 1, 1]}, sharex=True,figsize=[7,7])
a0.hist(nox_alpha_kx, bins, color='b', linestyle='dashdot',
         histtype='step', label='Variable Non-X-ray AGN',
         density=True, linewidth=1.5)
a0.hist(x_alpha_kx, bins, color='r', linestyle='solid',
         histtype='step', label='Variable X-ray AGN',
         density=True, linewidth=1.5)
a0.hist(allx_alpha_kx, bins, color='k', linestyle='dashed',
         histtype='step', label='X-ray AGN',
         density=True, linewidth=1.5)
#a0.xlabel(r'$\alpha_{KX}$')
a0.legend(loc='upper left')
a0.set_ylim(ymax=5)

allx_meanalpha = np.nanmean(allx_alpha_kx)
nox_meanalpha = np.nanmean(nox_alpha_kx)
x_meanalpha = np.nanmean(x_alpha_kx)
allx_stdalpha = np.nanstd(allx_alpha_kx)
nox_stdalpha = np.nanstd(nox_alpha_kx)
x_stdalpha = np.nanstd(x_alpha_kx)
#plt.xscale('log')

a1.vlines(allx_meanalpha, 0, 3.5, 'k', linestyle='dashed', label=r'Mean $\alpha_{KX}$')
a2.vlines(nox_meanalpha, 0, 3.5, 'b', linestyle='dashdot', label=r'Mean $\alpha_{KX}$')
a3.vlines(x_meanalpha, 0, 3.5, 'r', label=r'Mean $\alpha_{KX}$')
a1.axvspan(allx_meanalpha-allx_stdalpha, allx_meanalpha+allx_stdalpha, alpha=0.5, color='k')
a2.axvspan(nox_meanalpha-nox_stdalpha, nox_meanalpha+nox_stdalpha, alpha=0.5, color='b')
a3.axvspan(x_meanalpha-x_stdalpha, x_meanalpha+x_stdalpha, alpha=0.5, color='r')

a1.set_yticks([])
a2.set_yticks([])
a3.set_yticks([])

a1.legend(loc='upper right', frameon=False)
a2.legend(loc='upper right', frameon=False)
a3.legend(loc='upper right', frameon=False)

plt.subplots_adjust(hspace=0)
a3.set_xlabel(r'$\alpha_{KX}$')