#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 11 10:10:23 2018

Investigate properties of the variable sources

@author: ppxee
"""

### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
from astropy.table import Table #for handling tables
import numpy as np #for handling arrays
#import math
from astropy.stats import median_absolute_deviation
import vari_funcs #my module to help run code neatly
# Imports for catalogue matching
from astropy.coordinates import match_coordinates_sky
from astropy.coordinates import SkyCoord
from astropy import units as u
# Imports for luminosity distance
from astropy.cosmology import FlatLambdaCDM
# Set up cosmology
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
plt.close('all') #close any open plots

## Import fits table
#varydata = fits.open('variable_tables/variable_month_3.5_alldetails.fits')[1].data
#
## Extract magnitude table and error table
#mag = vari_funcs.mag5_months(varydata)
#magerr = vari_funcs.magerr5_months(varydata)
#
## Calculate excess variance
#excess = vari_funcs.sigmasq(mag, magerr)

# Get X-ray data
chandata = fits.open('variable_tables/chan_variable_month_3.5_alldetails.fits')[1].data
xmmdata = fits.open('variable_tables/xmm_variable_month_3.5_alldetails.fits')[1].data

#xrayvaryidx = chanvaryidx
#xrayvaryidx = np.append(chanvaryidx, xmmvaryidx)

# Calculate luminosity distances with both photometric and spectroscopic data
chanDLspec = cosmo.luminosity_distance(chandata['z_spec'])
chanDLspec = chanDLspec.to(u.cm)
chanDLspec[chanDLspec==0.0] = np.nan
chanDLphot = cosmo.luminosity_distance(chandata['z_m2'])
chanDLphot = chanDLphot.to(u.cm)
xmmDLspec = cosmo.luminosity_distance(xmmdata['z_spec'])
xmmDLspec[xmmDLspec==0.0] = np.nan
xmmDLspec = xmmDLspec.to(u.cm)
xmmDLphot = cosmo.luminosity_distance(xmmdata['z_m2'])
xmmDLphot = xmmDLphot.to(u.cm)

# Calculate luminosity
chanF = chandata['Soft_flux']
chanF[chanF==0] = np.nan
xmmF = xmmdata['CR(S)']*0.167e-14

xmmLspec = xmmF*4*np.pi*xmmDLspec.value**2
chanLspec = chanF*4*np.pi*chanDLspec.value**2
xmmLphot = xmmF*4*np.pi*xmmDLphot.value**2
chanLphot = chanF*4*np.pi*chanDLphot.value**2

#get variability measures
chanmodz = chandata['mod_z_score']
xmmmodz = xmmdata['mod_z_score']

chanmag = vari_funcs.mag5_months(chandata)
chanmagerr = vari_funcs.magerr5_months(chandata)
xmmmag = vari_funcs.mag5_months(xmmdata)
xmmmagerr = vari_funcs.magerr5_months(xmmdata)
## Change 99s to nans so they are ignored ###
mask = chanmagerr >= 99
chanmag[mask] = np.nan
chanmagerr[mask] = np.nan
mask = xmmmagerr >= 99
xmmmag[mask] = np.nan
xmmmagerr[mask] = np.nan

chanexcess = vari_funcs.normsigmasq(chanmag, chanmagerr)
xmmexcess = vari_funcs.normsigmasq(xmmmag, xmmmagerr)

#plt.figure()
#plt.scatter(chanLspec, chanmodz, label='Chandra')
#plt.xlabel('X-ray Luminosity (ergs/s)')
#plt.ylabel('Modified z-score')
#plt.title('Soft x-ray luminosity using spectroscopic redshifts')
#plt.scatter(xmmLspec, xmmmodz, label='XMM')
#plt.yscale('log')
#plt.xscale('log')
#plt.legend()
#
#plt.figure()
#plt.scatter(chanLphot, chanmodz, label='Chandra')
#plt.xlabel('X-ray Luminosity (ergs/s)')
#plt.ylabel('Modified z-score')
#plt.title('Soft x-ray luminosity using photometric redshifts')
#plt.scatter(xmmLphot, xmmmodz, label='XMM')
#plt.yscale('log')
#plt.xscale('log')
#plt.legend()

plt.figure()
plt.scatter(chanLspec, chanexcess, label='Chandra')
plt.xlabel('X-ray Luminosity (ergs/s)')
plt.ylabel('Excess Variance')
plt.title('Soft x-ray luminosity using spectroscopic redshifts')
plt.scatter(xmmLspec, xmmexcess, label='XMM')
#plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.tight_layout()

plt.figure()
plt.scatter(chanLphot, chanexcess, label='Chandra')
plt.xlabel('X-ray Luminosity (ergs/s)')
plt.ylabel('Excess Variance')
plt.title('Soft x-ray luminosity using photometric redshifts')
plt.scatter(xmmLphot, xmmexcess, label='XMM')
#plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.tight_layout()
