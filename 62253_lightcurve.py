#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 10:38:09 2019

Code to plot photometrically correct lightcurves for the flare

Done by subtracting 3 arcsec fluxes

@author: ppxee
"""

### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
import matplotlib.animation as animation
#plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'
from astropy.io import fits #for handling fits
#from astropy.table import Table #for handling tables
from photutils import CircularAperture, aperture_photometry
import numpy as np #for handling arrays
#import math
#from astropy.stats import median_absolute_deviation
import vari_funcs #my module to help run code neatly
plt.close('all') #close any open plots

### get mag flux tables ###
ktbdata = fits.open('mag_flux_tables/mag_flux_table_best_extra_clean_no06.fits')[1].data
jtbdata = fits.open('mag_flux_tables/mag_flux_table_best_J.fits')[1].data
htbdata = fits.open('mag_flux_tables/mag_flux_table_best_H.fits')[1].data

### change to just 62253 ###
kobdata = ktbdata[ktbdata['NUMBER_05B'] == 62243]
jobdata = jtbdata[jtbdata['NUMBER_05B'] == 57824]
hobdata = htbdata[htbdata['NUMBER_06B'] == 48107]

def extract_flux(band, binlen, tbdata):
    ### first get the flux and error values for each stack
    binarr = range(binlen)
    shortflux = np.empty(len(binarr))
    shortfluxerr = np.empty(len(binarr))
    for n in binarr:
        seout = fits.open('SE_outputs_shortstacks/08B_'+str(band)+'_output_cam1_'+str(binlen)+'bin_'+str(n)+'.fits')[1].data
        if band == 'J':
            seout = seout[seout['NUMBER'] == 13363]
        elif band == 'K':
            seout = seout[seout['NUMBER'] == 15100]
        elif band == 'H':
            seout = seout[seout['NUMBER'] == 11900]
        shortflux[n] = seout['FLUX_APER'][:,4]
        shortfluxerr[n] = seout['FLUXERR_APER'][:,4]
    
    ### Average non-flare flux ###
    if band == 'J':
        flux = vari_funcs.jflux5_stacks(tbdata).reshape(4)
        fluxerr = vari_funcs.jfluxerr5_stacks(tbdata).reshape(4)
        flux[3] = np.nan
        fluxerr[3] = np.nan
    elif band == 'K':
        flux = vari_funcs.flux5_stacks(tbdata).reshape(7)
        fluxerr = vari_funcs.fluxerr5_stacks(tbdata).reshape(7)
        flux[2] = np.nan
        fluxerr[2] = np.nan
    elif band == 'H':
        flux = vari_funcs.hflux5_stacks(tbdata).reshape(7)
        fluxerr = vari_funcs.hfluxerr5_stacks(tbdata).reshape(7)
        flux[2] = np.nan
        fluxerr[2] = np.nan
        
    meanflux = np.nanmean(flux)
    
    ### get flare flux ###
    flareflux = shortflux - meanflux#flux[0]
    
    ### Propogate errors ###
    errmeanflux = np.sqrt(np.nansum(np.square(fluxerr)))/binlen
    flarefluxerr = np.sqrt(np.square(shortfluxerr)+np.square(errmeanflux))
    
    ### Convert to magnitudes ###
    if band == 'J':
        flaremag = 30 - 2.5*np.log10(flareflux) + 0.91
    elif band == 'H':
        flaremag = 30 - 2.5*np.log10(flareflux) + 1.39 
    elif band == 'K':
        flaremag = 30 - 2.5*np.log10(flareflux) + 1.85 
    
    flaremagerr = 1.086/(flareflux/flarefluxerr)
    
    return flareflux, flarefluxerr, flaremag, flaremagerr

flarekflux, flarekfluxerr, flarekmag, flarekmagerr = extract_flux('K', 10, kobdata)
flarejflux, flarejfluxerr, flarejmag, flarejmagerr = extract_flux('J', 8, jobdata)
flarehflux, flarehfluxerr, flarehmag, flarehmagerr = extract_flux('H', 6, hobdata)

### Set X data ###
#Jxdata = np.array([18, 28, 29, 30, 33, 41, 42, 49])
#Hxdata = np.array([0,3,16,35,45,61])
#Kxdata = np.array([0, 1, 9, 11, 12, 13, 16, 21, 39, 56])
Jdata = fits.open('jdatetable.fits')[1].data
Hdata = fits.open('hdatetable.fits')[1].data
Kdata = fits.open('kdatetable.fits')[1].data

Jxdata = Jdata['median']
Jxfinish = Jdata['finish']
Jxstart = Jdata['start']
Jxerr = np.array([Jxdata-Jxstart, Jxfinish - Jxdata])

Hxdata = Hdata['median']
Hxfinish = Hdata['finish']
Hxstart = Hdata['start']
Hxerr = np.array([Hxdata-Hxstart, Hxfinish - Hxdata])

Kxdata = Kdata['median']
Kxfinish = Kdata['finish']
Kxstart = Kdata['start']
Kxerr = np.array([Kxdata-Kxstart, Kxfinish - Kxdata])
#plt.figure(figsize=[7,8])
#plt.errorbar(Kxdata, flarekflux, yerr=flarekfluxerr, fmt='ro')
#plt.errorbar(Jxdata, flarejflux, yerr=flarejfluxerr, fmt='bo')
#plt.gca().invert_yaxis()
#plt.tight_layout()

plt.figure(figsize=[7,8])
plt.errorbar(Kxdata, flarekmag, yerr=flarekmagerr, xerr=Kxerr, fmt='ro', label='K')
plt.errorbar(Jxdata, flarejmag, yerr=flarejmagerr, xerr=Jxerr, fmt='bo', label='J')
plt.errorbar(Hxdata, flarehmag, yerr=flarehmagerr, xerr=Hxerr, fmt='ko', label='H')
plt.gca().invert_yaxis()
plt.ylabel('AB Magnitude')
plt.xlabel('Days from first observation')
plt.legend()
plt.tight_layout()

#%% Get SN2010gz lightcurve ###
SNdata = fits.open('SN2010gz_photometry.fits')[1].data

def sort_mag_col(SN_col):
    ### set up arrays ###
    SN_mag = np.zeros(len(SN_col))
    SN_mag_err = np.zeros(len(SN_col))
    
    ### get values from strings ###
    for n, vals in enumerate(SN_col):
        if len(vals) == 0:
            SN_mag[n] = np.nan
            SN_mag_err[n] = np.nan
            continue
        elif vals[0] =='$':
            print('Upper limit detected')
            SN_mag[n] = np.nan
            SN_mag_err[n] = np.nan
            continue
        else:
            SN_mag[n] = vals[0:5]
            SN_mag_err[n] = vals[8:12]
    
    return SN_mag, SN_mag_err

SN_col = SNdata['g']
SN_g_mag, SN_g_mag_err = sort_mag_col(SN_col)
SN_col = SNdata['r']
SN_r_mag, SN_r_mag_err = sort_mag_col(SN_col)
SN_col = SNdata['z']
SN_z_mag, SN_z_mag_err = sort_mag_col(SN_col)

### Sort date stamps ###
SN_dates = SNdata['JD']
start_date = SN_dates[0]
SN_xdata = SN_dates - start_date
#SN_xdata = SN_xdata + 30

SN_z = 0.17

#%% Calculate apparent magnitudes ###
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u

### Define cosmology ###
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

z = 1.51 #approximately
DL = cosmo.luminosity_distance(z)
DL = DL.to(u.pc)
SN_DL = cosmo.luminosity_distance(SN_z)
SN_DL = SN_DL.to(u.pc)

M_J = flarejmag - 5*(np.log10(DL.value)-1)
M_H = flarehmag - 5*(np.log10(DL.value)-1)
M_K = flarekmag - 5*(np.log10(DL.value)-1)

SN_M_g = SN_g_mag - 5*(np.log10(SN_DL.value)-1)
#SN_M_r = SN_r_mag - 5*(np.log10(SN_DL.value)-1)
SN_M_z = SN_z_mag - 5*(np.log10(SN_DL.value)-1)

abinitdates = np.arange(67)/(1+z)
abinitdates = np.round(abinitdates,decimals=2)
abinitdates = abinitdates.astype(str)
abdates = abinitdates
SN_abdates = SN_xdata/(1+z)

### convert xdata to restframe days ###
restJxdata = Jxdata/(1+z)
restJxerr = Jxerr/(1+z)
restHxdata = Hxdata/(1+z)
restHxerr = Hxerr/(1+z)
restKxdata = Kxdata/(1+z)
restKxerr = Kxerr/(1+z)
rest_SN_xdata = SN_xdata/(1+SN_z)
#
#plt.plot(restJxdata, M_J, 'bo', label='J flare mag')
##plt.plot(restHxdata, M_H, 'ko', label='H flare mag')
#plt.plot(restKxdata, M_K, 'ro', label='K flare mag')

#%% calculate rest frame mags

M_g = flarejmag - 5*(np.log10(DL.value)-1) + 2.5*np.log10(1+z)
M_r = flarehmag - 5*(np.log10(DL.value)-1) + 2.5*np.log10(1+z)
M_z = flarekmag - 5*(np.log10(DL.value)-1) + 2.5*np.log10(1+z)

M_g_err = M_g * (flarejmagerr/flarejmag)
M_r_err = M_r * (flarehmagerr/flarehmag)
M_z_err = M_z * (flarekmagerr/flarekmag)

SN_M_g_rest = SN_g_mag - 5*(np.log10(SN_DL.value)-1) + 2.5*np.log10(1+SN_z)
SN_M_r_rest = SN_r_mag - 5*(np.log10(SN_DL.value)-1) + 2.5*np.log10(1+SN_z)
SN_M_z_rest = SN_z_mag - 5*(np.log10(SN_DL.value)-1) + 2.5*np.log10(1+SN_z)

SN_M_g_rest_err = SN_M_g_rest * (SN_g_mag_err/SN_g_mag)
SN_M_r_rest_err = SN_M_r_rest * (SN_r_mag_err/SN_r_mag)
SN_M_z_rest_err = SN_M_z_rest * (SN_z_mag_err/SN_z_mag)

plt.figure(figsize=[7,8])
plt.errorbar(restJxdata, M_g, yerr=M_g_err, xerr=restJxerr, fmt='bo', label='M_g flare mag')
plt.errorbar(restHxdata, M_r, yerr=M_r_err, xerr=restHxerr, fmt='ko', label='M_r flare mag')
plt.errorbar(restKxdata, M_z, yerr=M_z_err, xerr=restKxerr, fmt='ro', label='M_z flare mag')

plt.errorbar(rest_SN_xdata, SN_M_g_rest, yerr=SN_M_g_rest_err, fmt='b*', label='2010gz g mag')
plt.errorbar(rest_SN_xdata, SN_M_r_rest, yerr=SN_M_r_rest_err, fmt='k*', label='2010gz r mag')
plt.errorbar(rest_SN_xdata, SN_M_z_rest, yerr=SN_M_z_rest_err, fmt='r*', label='2010gz z mag')
#
#plt.errorbar(restJxdata, M_g, yerr=M_g_err, xerr=restJxerr, fmt='bo', alpha=0.5)#label='M_g flare mag')
#plt.errorbar(restHxdata, M_r, yerr=M_r_err, xerr=restHxerr, fmt='ko', alpha=0.5)#label='M_r flare mag')
#plt.errorbar(restKxdata, M_z, yerr=M_z_err, xerr=restKxerr, fmt='ro', alpha=0.5)#label='M_z flare mag')
#
#plt.plot(restJxdata, M_g,'bo', label='M_g flare mag')
#plt.plot(restHxdata, M_r, 'ko', label='M_r flare mag')
#plt.plot(restKxdata, M_z, 'ro', label='M_z flare mag')
#
#plt.errorbar(rest_SN_xdata, SN_M_g_rest, yerr=SN_M_g_rest_err, fmt='b*', alpha=0.5)#label='2010gz g mag')
#plt.errorbar(rest_SN_xdata, SN_M_r_rest, yerr=SN_M_r_rest_err, fmt='k*', alpha=0.5)#label='2010gz r mag')
#plt.errorbar(rest_SN_xdata, SN_M_z_rest, yerr=SN_M_z_rest_err, fmt='r*', alpha=0.5)#label='2010gz z mag')
#
#plt.plot(rest_SN_xdata, SN_M_g_rest, 'b*', label='2010gz g mag')
#plt.plot(rest_SN_xdata, SN_M_r_rest, 'k*', label='2010gz r mag')
#plt.plot(rest_SN_xdata, SN_M_z_rest, 'r*', label='2010gz z mag')

plt.legend()
plt.ylabel('Restframe magnitude')
plt.xlabel('Restframe days from first observation')
plt.gca().invert_yaxis()
plt.tight_layout()

