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
chanfull = fits.open('mag_flux_tables/month_chan_mag_flux_table_DR11.fits')[1].data
xmmfull = fits.open('mag_flux_tables/month_xmm_mag_flux_table_DR11.fits')[1].data

# Calculate luminosity
chanF = chandata['Soft_flux']
chanF[chanF==0] = np.nan
fullchanF = chanfull['Soft_flux']
fullchanF[fullchanF==0] = np.nan
xmmF = xmmdata['CR(S)']*0.167e-14
fullxmmF = xmmfull['CR(S)']*0.167e-14

def spec_lum(tbdata, flux):
    DLspec = cosmo.luminosity_distance(tbdata['z_spec'])
    DLspec = DLspec.to(u.cm)
    DLspec[DLspec==0.0] = np.nan
    Lspec = flux*4*np.pi*DLspec.value**2
    return Lspec

def phot_lum(tbdata, flux):
    DLphot = cosmo.luminosity_distance(tbdata['z_m2'])
    DLphot = DLphot.to(u.cm)
    DLphot[DLphot==0.0] = np.nan
    Lphot = flux*4*np.pi*DLphot.value**2
    return Lphot


xmmLspec = spec_lum(xmmdata, xmmF)
chanLspec = spec_lum(chandata, chanF)
xmmLphot = phot_lum(xmmdata, xmmF)
chanLphot = phot_lum(chandata, chanF)
fullxmmLspec = spec_lum(xmmfull, fullxmmF)
fullchanLspec = spec_lum(chanfull, fullchanF)
fullxmmLphot = phot_lum(xmmfull, fullxmmF)
fullchanLphot = phot_lum(chanfull, fullchanF) 

#get variability measures
chanmag = vari_funcs.mag5_months(chandata)
chanmagerr = vari_funcs.magerr5_months(chandata)
xmmmag = vari_funcs.mag5_months(xmmdata)
xmmmagerr = vari_funcs.magerr5_months(xmmdata)
fullchanmag = vari_funcs.mag5_months(chanfull)
fullchanmagerr = vari_funcs.magerr5_months(chanfull)
fullxmmmag = vari_funcs.mag5_months(xmmfull)
fullxmmmagerr = vari_funcs.magerr5_months(xmmfull)

## Change 99s to nans so they are ignored ###
mask = chanmagerr >= 99
chanmag[mask] = np.nan
chanmagerr[mask] = np.nan
mask = xmmmagerr >= 99
xmmmag[mask] = np.nan
xmmmagerr[mask] = np.nan
mask = fullchanmagerr >= 99
fullchanmag[mask] = np.nan
fullchanmagerr[mask] = np.nan
mask = fullxmmmagerr >= 99
fullxmmmag[mask] = np.nan
fullxmmmagerr[mask] = np.nan

chanexcess = vari_funcs.normsigmasq(chanmag, chanmagerr)
xmmexcess = vari_funcs.normsigmasq(xmmmag, xmmmagerr)
fullchanexcess = vari_funcs.normsigmasq(fullchanmag, fullchanmagerr)
fullxmmexcess = vari_funcs.normsigmasq(fullxmmmag, fullxmmmagerr)
fullxmmexcess[fullxmmexcess<=0.0] = np.nan
fullchanexcess[fullchanexcess<=0.0] = np.nan

plt.figure()
plt.scatter(fullchanLspec, fullchanexcess, label='Chandra')
plt.scatter(fullxmmLspec, fullxmmexcess, label='XMM')
plt.plot(chanLspec, chanexcess, 'ro', markersize=10, label='Varying Chandra')
plt.plot(xmmLspec, xmmexcess, 'go', markersize=10, label='Varying XMM')
plt.xlabel('X-ray Luminosity (ergs/s)')
plt.ylabel('Excess Variance')
plt.title('Soft x-ray luminosity using spectroscopic redshifts')
plt.ylim(ymin=0.00002)
plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.tight_layout()

plt.figure()
plt.scatter(fullchanLphot, fullchanexcess, c=chanfull['KMAG_20'], label='Chandra')
plt.scatter(fullxmmLphot, fullxmmexcess, c=xmmfull['KMAG_20'], label='XMM')
plt.plot(chanLphot, chanexcess, 'ro', markersize=10, mfc='None', label='Varying Chandra')
plt.plot(xmmLphot, xmmexcess, 'go', markersize=10, mfc='None', label='Varying XMM')
plt.xlabel('X-ray Luminosity (ergs/s)')
plt.ylabel('Excess Variance')
plt.title('Soft x-ray luminosity using photometric redshifts')
plt.ylim(ymin=0.00002)
plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.tight_layout()
