#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  9 10:06:10 2018

@author: ppxee
"""

### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
import numpy as np #for handling arrays
#from astropy.convolution import Gaussian2DKernel
#from astropy.convolution import convolve
import vari_funcs
#from scipy import ndimage
#import math
#from astropy.stats import median_absolute_deviation

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
    
stars = fits.open('stars_mag_flux_table.fits')
sdata = stars[1].data
#hdr08B = fits.getheader('UDS_08B_K.fits') # random year (same in all)
#const = -hdr08B['CD1_1'] # constant that defines unit conversion for FWHM

colname = 'FLUX_RADIUS_'#'FWHM_WORLD_'
#data = sem05B[colname][:,1]
semesters = ['05B', '06B', '07B', '08B', '09B', '10B', '11B', '12B']
#Extract the flux radii and remove negative values
avgFWHM = np.zeros(8)
for n, sem in enumerate(semesters):
    colnames = colname+sem
    avgFWHM[n] = 2*np.median(sdata[colnames][:,1])

vari_funcs.avg_lightcurve(avgFWHM)
#### Find maximum FWHM as this is what all the others willl become ###
#aimind = np.argmax(avgFWHM)
##
#### Convert FWHM into a sigma ###
#sigmaold = np.array([fluxrad2sigma(fwhm) for fwhm in avgFWHM])
#sigmabroad = sigmaold[aimind]
#
#### Find required sigma ###
## sigker^2 = sigbroad^2 - signar^2
#sigmakernel = np.array([np.sqrt(sigmabroad**2 - sigma**2) for sigma in sigmaold])
