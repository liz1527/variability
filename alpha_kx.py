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

def func(x):
    ### power law to integrate ###
    return x**(-0.9)

def get_L_k(tbdata):
    ### get magnitude ###
    kmag = tbdata['M_K_z_p']
    
    ### convert to luminosity ###
    L_k = 10**((34.1-kmag)/2.5)
    L_k = L_k * u.W * (u.Hz)**-1 # Units of W/Hz
    
    L_k = L_k.to((u.erg) * (u.s)**-1 * (u.Hz)**-1, equivalencies=u.spectral())
    
    return L_k

allx_L_k = get_L_k(fullxray)
nox_L_k = get_L_k(noxvarydata)
x_L_k = get_L_k(xvarydata)

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

def get_xray_L_2(tbdata, Xray=True, band='Hard'):
    z = vari_funcs.get_z(tbdata)
#    z[z>4] = np.nan
    if band=='Hard': 
        upplim = 10 ## set band limits in keV
        lowlim = 2
        if Xray == True: # if it is an X-ray source, get flux from catalogue
            flux = tbdata['Hard_flux'] 
            flux# Units of erg cm**-2 s**-1
        else: # if it is non X-ray - use the upper limit
            flux = np.zeros(len(tbdata))
            flux += 6.5e-16 # Units of erg cm**-2 s**-1
    elif band=='Full': 
        upplim = 10
        lowlim = 0.5
        if Xray == True:
            flux = tbdata['Full_flux'] # Units of erg cm**-2 s**-1
        else:
            flux = 4.4e-16 # Units of erg cm**-2 s**-1
    elif band=='Soft': 
        upplim = 2
        lowlim = 0.5
        if Xray == True:
            flux = tbdata['Soft_flux'] # Units of erg cm**-2 s**-1
        else:
            flux = 1.4e-16 # Units of erg cm**-2 s**-1
    elif band=='Uhrd': 
        upplim = 10
        lowlim = 5
        if Xray == True:
            flux = tbdata['Uhrd_flux'] # Units of erg cm**-2 s**-1
        else:
            flux = 9.2e-15 # Units of erg cm**-2 s**-1
            
    ### Add units ###
    flux = flux* (u.erg) * (u.cm)**-2 * (u.s)**-1 
    upplim = upplim * u.keV
    lowlim = lowlim * u.keV
    
    ### redshift limits ###
#    upplim = upplim/(1+z)
#    lowlim = lowlim/(1+z)
    
    ### get integrated flux density ###
#    result = integrate.quad(func, lowlim.value, upplim.value)
#    denom = result[0] 
#    print(denom)
    denom = ((upplim**(0.1))/(0.1)) - ((lowlim**(0.1))/(0.1))
    print(denom)
    
    ### use this and flux value to get the power law constant ###
    const = flux / denom
    
    ### calculate flux density ###
    nu = 2 * u.keV # 2kev is value to evaluate at
    F_2 = const * (nu**(-0.9))
    
    ### get luminosity distance ###
    DL = cosmo.luminosity_distance(z) # need to sort out units
    DL = DL.to(u.cm)
    
    ### calculate luminosity density ###
    L_2 = 4 * np.pi * (DL**2) * F_2#const * (nu**-0.9)
    
    L_2 = L_2.to((u.erg) * (u.s)**-1 * (u.Hz)**-1, equivalencies=u.spectral())
    
    #.to(u.W * (u.Hz)**-1, equivalencies=u.spectral())
    
    L_2[L_2==0] = np.nan
    
    return L_2, F_2, flux, z #L_2_w_Hz

allx_L_2, allx_F_2, allx_flux, allx_z = get_xray_L_2(fullxray)
nox_L_2, nox_F_2, nox_flux, nox_z = get_xray_L_2(noxvarydata, Xray=False)
x_L_2, x_F_2, x_flux, x_z = get_xray_L_2(xvarydata)

#print('Start flux = '+str(allx_flux[10]))#+' erg/cm**2/s')
#print('z = '+str(allx_z[10]))
#print('Flux density at 2keV = '+str(allx_F_2[10]))#+' erg/cm**2/s/keV')
print('Luminosity density at 2keV = '+str(allx_L_2[10]))#+' erg/s/eV')
print('Luminosity density at K = '+str(allx_L_k[10]))#+' erg/s/eV')

def calc_alpha_kx(L_k, L_2, xband=2*u.keV, optband=1.6*u.um):
    
    ### convert units ###
    xband = xband.to(u.um, equivalencies=u.spectral())
    
    optband = optband.to(u.um, equivalencies=u.spectral())
    
    numer = np.log(optband.value * L_k.value) - np.log(xband.value * L_2.value)
    denom = np.log(optband.value) - np.log(xband.value)
    
    alpha = -(numer/denom) + 1
    
#    alpha2 = -0.3838 * (np.log(L_2.value)/np.log(L_k.value))
    
    mask = np.isinf(alpha)
    alpha[mask] = np.nan
    return alpha

allx_alpha_kx = calc_alpha_kx(allx_L_k, allx_L_2)
nox_alpha_kx = calc_alpha_kx(nox_L_k, nox_L_2)
x_alpha_kx = calc_alpha_kx(x_L_k, x_L_2)

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