#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 14:10:37 2019

Code to find alphs-ox distrutions of variable sources in chandra region

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

def get_L_O(tbdata):
    z = vari_funcs.get_z(tbdata)
    z[z<0.5] = np.nan
    z[z>4] = np.nan
    
    ### get luminosity distance ###
    DL = cosmo.luminosity_distance(z) # need to sort out units
    DL = DL.to(u.cm)
    
    ### get U band magnitude ###
    umag = tbdata['M_U_z_p']
    
    L_U = 10**((34.1-umag)/2.5)
    L_U = L_U * u.W * (u.Hz)**-1 # Units of W/Hz
    
    L_U = L_U.to((u.erg) * (u.s)**-1 * (u.Hz)**-1, equivalencies=u.spectral())
    
    ### Use frequency ratios to find Omag ###
    Ofreq = constants.c/(2500 * u.AA)
    Ufreq = constants.c/(3743 * u.AA)
    
    L_O = L_U * (Ufreq/Ofreq)     # assuming flux propotional to freq**-1
#    F_O = F_O.to((u.erg) * (u.s)**-1 * (u.cm)**2 * (u.Hz)**-1, equivalencies=u.spectral()) # units of ergs/s/cm**"/Hz
    
    
#    L_O = L_O.to((u.erg) * (u.s)**-1 * (u.Hz)**-1, equivalencies=u.spectral())
    
    return L_O

allx_L_O = get_L_O(fullxray)
nox_L_O = get_L_O(noxvarydata)
x_L_O = get_L_O(xvarydata)

##_, bins, _ = plt.hist([nox_L_O, x_L_O, allx_L_O], bins=50, color=['b','r','k'], #linestyle=['dashed','dashed','dashed'],
##         histtype='step', label=['Variable Non-X-ray AGN','Variable X-ray AGN','X-ray AGN'],
##         normed=True)
#bins=np.logspace(16,np.log10(1.2e25))
#plt.figure(figsize=[7,5])
#plt.hist(nox_L_O, bins, color='b', linestyle='dashdot',
#         histtype='step', label='Variable Non-X-ray AGN')#,
##         normed=True, linewidth=1.5)
#plt.hist(x_L_O, bins, color='r', linestyle='solid',
#         histtype='step', label='Variable X-ray AGN')#,
##         normed=True, linewidth=1.5)
#plt.hist(allx_L_O, bins, color='k', linestyle='dashed',
#         histtype='step', label='X-ray AGN')#,
##         normed=True, linewidth=1.5)
#plt.xscale('log')

### Need to calculate the monochromatic luminosity for 2keV ###
### First need to assume power law for flux density and find constant ###
### Then sub this equation for flux density into luminosity density eq ###

def get_xray_L_2(tbdata, Xray=True, band='Hard'):
    z = vari_funcs.get_z(tbdata)
    z[z<0.5] = np.nan
    z[z>4] = np.nan
    
    ### get luminosity distance ###
    DL = cosmo.luminosity_distance(z) # need to sort out units
    DL = DL.to(u.cm)
    
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
    
    
    ### calculate luminosity density ###
    L_2 = 4 * np.pi * (DL**2) * F_2#const * (nu**-0.9)
    
    L_2 = L_2.to((u.erg) * (u.s)**-1 * (u.Hz)**-1, equivalencies=u.spectral())
        
    L_2[L_2==0] = np.nan
    
    return L_2, F_2, flux, z #L_2_w_Hz

allx_L_2, allx_F_2, allx_flux, allx_z = get_xray_L_2(fullxray)
nox_L_2, nox_F_2, nox_flux, nox_z = get_xray_L_2(noxvarydata, Xray=False)
x_L_2, x_F_2, x_flux, x_z = get_xray_L_2(xvarydata)

#print('Start flux = '+str(allx_flux[10]))#+' erg/cm**2/s')
#print('z = '+str(allx_z[10]))
#print('Flux density at 2keV = '+str(allx_F_2[10]))#+' erg/cm**2/s/keV')
print('Luminosity density at 2keV = '+str(allx_L_2[10]))#+' erg/s/eV')
print('Luminosity density at O = '+str(allx_L_O[10]))#+' erg/s/eV')

def calc_alpha_Ox(L_O, L_2, xband=2*u.keV, optband=1.6*u.um):
    
    ### convert units ###
    xband = xband.to(u.um, equivalencies=u.spectral())
    
    optband = optband.to(u.um, equivalencies=u.spectral())
    
#    numer = np.log(optband.value * L_O.value) - np.log(xband.value * L_2.value)
#    denom = np.log(optband.value) - np.log(xband.value)
#    
#    alpha = -(numer/denom) + 1
    
    alpha = -0.3838 * (np.log10(L_2.value/L_O.value))
    
    mask = np.isinf(alpha)
    alpha[mask] = np.nan
    return alpha

allx_alpha_Ox = calc_alpha_Ox(allx_L_O, allx_L_2)
nox_alpha_Ox = calc_alpha_Ox(nox_L_O, nox_L_2)
x_alpha_Ox = calc_alpha_Ox(x_L_O, x_L_2)


### Calculate a-ox lines for plot ###
x = np.logspace(25,32,10)
y1 = 10 ** (np.log10(x)-1/0.3838)
y2 = 10 ** (np.log10(x) - 2/0.3838)
#y3 = 10 ** (np.log10(x) - 3/0.3838)

plt.figure()
upplims = nox_L_2.value/3
plt.plot(allx_L_O, allx_L_2,'k+', label='Non-Variable X-ray AGN',zorder=1)
dot = plt.plot(nox_L_O, nox_L_2,'bo', label='Variable Non-X-ray AGN')
plt.errorbar(nox_L_O.value, nox_L_2.value,yerr=upplims, fmt='bo', #mfc='None', #markersize=10,
             uplims=True,zorder=3)#, label='Variable Non-X-ray AGN'
plt.plot(x_L_O, x_L_2,'ro', label='Variable X-ray AGN',zorder=2)
plt.plot(x,y1, 'k', label=r'$\alpha_{OX} = 1$',zorder=0)
plt.plot(x,y2, 'k--', label=r'$\alpha_{OX} = 2$',zorder=0)
#plt.plot(x,y3, label=r'$\alpha_{OX} = 3$')

plt.xlim(xmin=1e26, xmax = 2e31)
plt.ylim(ymin=1e23, ymax = 2e28)
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$L_{O}\ (erg\ s^{-1}\ Hz^{-1})$')
plt.ylabel(r'$L_{2 keV}\ (erg\ s^{-1}\ Hz^{-1})$')
#arrow = plt.scatter(0 ,0, c='b',marker=r'${\downarrow}$',s=100, label='Variable Non-X-ray AGN' )
plt.legend(fontsize = 10)
plt.tight_layout()

plt.figure()
_, bins, _ = plt.hist([nox_alpha_Ox, x_alpha_Ox, allx_alpha_Ox], bins=50)#, color=['b','r','k'], #linestyle=['dashed','dashed','dashed'],
#         histtype='step', label=['Variable Non-X-ray AGN','Variable X-ray AGN','X-ray AGN']),
#         normed=True)
plt.close()
f, (a0, a1, a2, a3) = plt.subplots(4, 1, gridspec_kw={'height_ratios': [10, 1, 1, 1]}, sharex=True,figsize=[7,7])
a0.hist(nox_alpha_Ox, bins, color='b', linestyle='dashdot',
         histtype='step', label='Variable Non-X-ray AGN',
         density=True, linewidth=1.5)
a0.hist(x_alpha_Ox, bins, color='r', linestyle='solid',
         histtype='step', label='Variable X-ray AGN',
         density=True, linewidth=1.5)
a0.hist(allx_alpha_Ox, bins, color='k', linestyle='dashed',
         histtype='step', label='X-ray AGN',
         density=True, linewidth=1.5)
#a0.xlabel(r'$\alpha_{OX}$')
a0.legend(loc='upper left')
a0.set_ylim(ymax=6)

allx_meanalpha = np.nanmean(allx_alpha_Ox)
nox_meanalpha = np.nanmean(nox_alpha_Ox)
x_meanalpha = np.nanmean(x_alpha_Ox)
allx_stdalpha = np.nanstd(allx_alpha_Ox)
nox_stdalpha = np.nanstd(nox_alpha_Ox)
x_stdalpha = np.nanstd(x_alpha_Ox)
#plt.xscale('log')

a1.vlines(allx_meanalpha, 0, 3.5, 'k', linestyle='dashed', label=r'Mean $\alpha_{OX}$')
a2.vlines(nox_meanalpha, 0, 3.5, 'b', linestyle='dashdot', label=r'Mean $\alpha_{OX}$')
a3.vlines(x_meanalpha, 0, 3.5, 'r', label=r'Mean $\alpha_{OX}$')
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
a3.set_xlabel(r'$\alpha_{OX}$')