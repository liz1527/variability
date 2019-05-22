#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  8 11:51:25 2018

Code to investigate individual flux curves of stars to see why median flux 
curve looks so wrong

@author: ppxee
"""

### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
#from astropy.table import Table #for handling tables
import numpy as np #for handling arrays
from astropy.stats import median_absolute_deviation
import vari_funcs #my module to help run code neatly
plt.close('all') #close any open plots
#from numpy.lib.recfunctions import append_fields
from photutils import CircularAperture, aperture_photometry

hdr08B = fits.getheader('Images/UDS_08B_K.fits') # random year (same in all)
const = -hdr08B['CD1_1'] # constant that defines unit conversion for FWHM
semesters = ['05B','06B', '07B', '08B', '09B', '10B', '11B', '12B']

### get data and create flux and errors array ###
sdata = fits.open('mag_flux_tables/stars_mag_flux_table.fits')[1].data
sflux = vari_funcs.flux5_stacks(sdata)
smag = vari_funcs.mag5_stacks(sdata)
sflux, sdata = vari_funcs.noneg(sflux, sdata)
smag = vari_funcs.mag5_stacks(sdata)

### set up normalised flux ###
sfluxn = np.copy(sflux)
sfluxn = vari_funcs.normalise_flux(sfluxn)

### Limit to same magnitude range as PSFS ###
avgmag = np.nanmean(smag, axis=1)
mask1 = avgmag > 15
mask2 = avgmag <= 19

### Remove any variable stars from sample ###
mad = median_absolute_deviation(sfluxn, axis=1)
mask3 = mad < 0.01

### Apply mask ###
mask = mask1*mask2*mask3
mask = mask.astype('bool')
sflux = sflux[mask]
sfluxn = sfluxn[mask]
sdata = sdata[mask]
smag = smag[mask]
sfluxerr = vari_funcs.fluxerr5_stacks(sdata)


### Generate a random interger to look at flux and FWHM curves of ###
obnum = np.random.randint(0,1599)
obID = sdata['NUMBER_05B'][obnum]
obflux = sflux[obnum,:]
obfwhm = np.array([sdata['FWHM_WORLD_05B'][obnum], sdata['FWHM_WORLD_06B'][obnum],
          sdata['FWHM_WORLD_07B'][obnum], sdata['FWHM_WORLD_08B'][obnum],
          sdata['FWHM_WORLD_09B'][obnum], sdata['FWHM_WORLD_10B'][obnum],
          sdata['FWHM_WORLD_11B'][obnum], sdata['FWHM_WORLD_12B'][obnum]])*3600
obfr = np.array([sdata['FLUX_RADIUS_05B'][obnum][1], sdata['FLUX_RADIUS_06B'][obnum][1],
          sdata['FLUX_RADIUS_07B'][obnum][1], sdata['FLUX_RADIUS_08B'][obnum][1],
          sdata['FLUX_RADIUS_09B'][obnum][1], sdata['FLUX_RADIUS_10B'][obnum][1],
          sdata['FLUX_RADIUS_11B'][obnum][1], sdata['FLUX_RADIUS_12B'][obnum][1]])*0.1342
obmag = smag[obnum,:]
obfluxerr = sfluxerr[obnum,:]
obfluxn = sfluxn[obnum,:]

### Generate plots ###
t = np.linspace(1,8,num=8)
plt.figure(1, figsize=[9,9])

aperflux = np.empty(8)
plt.figure(2, figsize=[4,8])
for rad in [0.25]:#,0.5,0.75,1,1.25,1.5]:
    ### Set up aperture ### 
    pixelr = (rad/3600) / const
    centre = [29,29]
    aperture = CircularAperture(centre, pixelr)
    
    for n, sem in enumerate(semesters):
        vigtb = fits.open('star_stamps_tables/small_'+sem+'_star_stamps_table.fits')[1].data
        obvig = vigtb['VIGNET'][vigtb['NUMBER']==obID]
        obvig = obvig.reshape([59,59])
        plt.figure(2)
        plt.subplot(4,2,n+1)
        plt.imshow(np.log(obvig))
        vari_funcs.no_ticks()
        plt.title(sem)
        
        ### Determine flux within 3 arcsec apertures ###
        phot = aperture_photometry(obvig, aperture)
        aperflux[n] = phot['aperture_sum'][0]
    
    plt.figure(1)
    plt.subplot(326)
    plt.scatter(t, aperflux)
    plt.xticks(t, semesters)
    plt.xlabel('Semester')
    plt.ylabel('Aperture flux on vignet')
    plt.tight_layout()

plt.subplot(321)
plt.errorbar(t, obflux, obfluxerr,fmt='o')
plt.xticks(t, semesters)
plt.xlabel('Semester')
plt.ylabel('Flux')

#plt.figure()
plt.subplot(322)
plt.scatter(t, obfwhm)
plt.xticks(t, semesters)
plt.xlabel('Semester')
plt.ylabel('FWHM')

#plt.figure()
plt.subplot(323)
plt.scatter(t, obmag)
plt.xticks(t, semesters)
plt.xlabel('Semester')
plt.ylabel('K Band Magnitude')

#plt.figure()
plt.subplot(324)
plt.scatter(t, obfr)
plt.xticks(t, semesters)
plt.xlabel('Semester')
plt.ylabel('0.5 Flux Radius')
plt.tight_layout()

#plt.figure()
plt.subplot(325)
plt.scatter(t, obfluxn)
plt.xticks(t, semesters)
plt.xlabel('Semester')
plt.ylabel('Normalised Flux')
plt.tight_layout()

##plt.figure()
#plt.subplot(326)
#plt.scatter(t, aperflux)
#plt.xticks(t, semesters)
#plt.xlabel('Semester')
#plt.ylabel('3" aperture flux on vignet')
#plt.tight_layout()