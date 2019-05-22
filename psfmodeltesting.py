#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 19 15:16:54 2018

@author: ppxee
"""
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy.stats import norm
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve
from photutils import CircularAperture
from photutils import aperture_photometry
plt.close('all')

def radial_profile(data, center):
    y, x = np.indices((data.shape)) #create coordinate grid
    r = np.sqrt((x - center[0])**2 + (y - center[1])**2) #get radius values for grid
    r = r.astype(np.int)

    tbin = np.bincount(r.ravel(), data.ravel()) # counts number of times value
                                                # of radius occurs in the psf
                                                # weighted by the data
    nr = np.bincount(r.ravel()) # counts number of radii values in psf
    radialprofile = tbin / nr # as weighted is r*data then get profile by 
                              # dividing by unweighted counts of r values.
    return radialprofile 

def FWHM2sigma(FWHM, const):
    ''' Function to convert the FWHM of a distribution into a sigma for that
    distribution. It assumes the distribution is gaussian.
    Input:
        FWHM = Full width half maximum of a distriubtution (in my case usually
                of an object from SExtractor)
    Output:
        sigma = standard deviation value of a guassian distribution with the 
                given FWHM. This roughly equates to the psf of the object. '''
    FWHM /= const
    return FWHM/np.sqrt(8*np.log(2))

def fluxrad2sigma(fluxrad):
    return fluxrad/np.sqrt(8*np.log(2))

def convolve_psf(sem, psf, newpsf, sigmakernel):
    ## Open image ###
    im05B = psf[sem]
    
    ## Convolve Image ###
    print('Convolving ', sem)
    kernel = Gaussian2DKernel(sigmakernel)
    newpsf[sem] = convolve(im05B, kernel, normalize_kernel=True) 

def convolve_one_psf(psf, sigmakernel):
    
    ## Convolve Image ###
    print('Convolving ', sem)
    kernel = Gaussian2DKernel(sigmakernel)
    newpsf = convolve(psf, kernel, normalize_kernel=True) 
    return newpsf
    
sdata = fits.open('mag_flux_tables/stars_mag_flux_table_extra.fits')[1].data
oldsdata = fits.open('mag_flux_tables/stars_mag_flux_table.fits')[1].data
hdr08B = fits.getheader('Images/UDS_08B_K.fits') # random year (same in all)
const = -hdr08B['CD1_1'] # constant that defines unit conversion for FWHM

colname = 'FWHM_WORLD_'
#data = sem05B[colname][:,1]
semesters = ['05B','10B']#['05B', '06B', '07B', '08B', '09B', '10B', '11B', '12B']
centre = [63,63]

avgFWHM = np.zeros(8)
oldavgFWHM = np.zeros(8)
psf = {}
oldpsf = {}
for n, sem in enumerate(semesters):
    # for new
    colnames = colname+sem
    mag = sdata['MAG_APER_'+sem][:,4]
    mask1 = mag > 15 #removes saturated
    mask2 = mag < 20 #removes very faint stars
    mask = mask1 * mask2
    tempsdata = sdata[mask]
    avgFWHM[n] = np.median(tempsdata[colnames]) #* 3600
    if sem == '10B':
        psf[sem] = fits.open('PSFs/limited_'+sem+'_K_PSF.fits')[0].data
    else:
        psf[sem] = fits.open('PSFs/extra_'+sem+'_K_PSF.fits')[0].data
    # for old
    oldmag = sdata['MAG_APER_'+sem][:,4]
    mask1 = mag > 15 #removes saturated
    mask2 = mag < 20 #removes very faint stars
    oldmask = mask1 * mask2
    tempsdata = oldsdata[oldmask]
    oldavgFWHM[n] = np.median(tempsdata[colnames]) #* 3600
    oldpsf[sem] = fits.open('PSFs/limited_'+sem+'_K_PSF.fits')[0].data
    
   
### Find maximum FWHM as this is what all the others willl become ###
aimind = np.argmax(oldavgFWHM)
aimsem = semesters[aimind]
aimpsf = oldpsf[aimsem]

### Convert FWHM into a sigma ###
sigmaold = np.array([FWHM2sigma(fwhm, const) for fwhm in oldavgFWHM])
sigmabroad = sigmaold[aimind]

phot = {}
flux = np.zeros(8)
oldphot = {}
oldflux = np.zeros(8)
aperture = CircularAperture(centre, 1.5)

for n, sem in enumerate(semesters):
    ### Determine flux within 3 arcsec apertures ###
    phot[sem] = aperture_photometry(psf[sem], aperture)
    flux[n] = phot[sem]['aperture_sum'][0]
    oldphot[sem] = aperture_photometry(oldpsf[sem], aperture)
    oldflux[n] = oldphot[sem]['aperture_sum'][0]

plt.figure()
plt.plot(oldflux, 'bs')
plt.plot(flux,'ro')
plt.ylabel('Flux within 3 arcsec aperture')
plt.figure()
ax = plt.subplot(111)
plt.plot(oldavgFWHM, 'bs')
plt.plot(avgFWHM, 'ro')
plt.ylabel('FWHM')
ax.invert_yaxis()
### testing the extra factor method ###
tests = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4]
r = np.arange(0,90,1) * const * 3600 # define radius values

flux10B = oldflux[aimind]

newpsf = {}
newphot = {}
newflux = np.zeros(8)

plt.figure()
for n, sem in enumerate(semesters):
    if sem == semesters[aimind]:
        continue
    
    sigma = sigmaold[n]
    testphot = np.zeros(len(tests))
    
    for m, extra in enumerate(tests):
        
        sigmakernel = np.sqrt(sigmabroad**2 - sigma**2) - extra

        new = convolve_one_psf(psf[sem], sigmakernel)    
        
        ### Determine flux within 3 arcsec apertures ###
        singlephot = aperture_photometry(newpsf[sem], aperture)
        singleflux = singlephot['aperture_sum'][0]
        
#    # find radial profiles
#        radialprofile = radial_profile(newpsf[sem], centre)
#    #    radialprofile = normalise(radialprofile)
#        sqrtrp = np.sqrt(radialprofile)
        
    #    plt.figure()
        plt.plot(oldflux, 'bs')
        plt.plot(newflux,'o', label=extra)
        plt.ylabel('Flux')
        plt.legend()