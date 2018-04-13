#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 14:00:59 2017

Code for degrading images - DO NOT RUN ON PPPXEE

@author: ppxee
"""
### Import required libraries ###
#import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
import numpy as np #for handling arrays
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve
import vari_funcs_no06 #my module to help run code neatly
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

sem05B = fits.open('~/Documents/Images/SE_outputs_yearstacks/05B_output.fits')[1].data
sem07B = fits.open('~/Documents/Images/SE_outputs_yearstacks/07B_output.fits')[1].data
sem08B = fits.open('~/Documents/Images/SE_outputs_yearstacks/08B_output.fits')[1].data
sem09B = fits.open('~/Documents/Images/SE_outputs_yearstacks/09B_output.fits')[1].data
sem10B = fits.open('~/Documents/Images/SE_outputs_yearstacks/10B_output.fits')[1].data
sem11B = fits.open('~/Documents/Images/SE_outputs_yearstacks/11B_output.fits')[1].data
sem12B = fits.open('~/Documents/Images/SE_outputs_yearstacks/12B_output.fits')[1].data
hdr08B = fits.getheader('UDS_08B_K.fits') # random year (same in all)
const = -hdr08B['CD1_1'] # constant that defines unit conversion for FWHM

colname = 'FWHM_WORLD'

#Put data in array
avgFWHM = np.array([np.mean(sem05B[colname]), np.mean(sem07B[colname]), 
                 np.mean(sem08B[colname]), np.mean(sem09B[colname]), 
                 np.mean(sem10B[colname]), np.mean(sem11B[colname]),
                 np.mean(sem12B[colname])])
vari_funcs_no06.avg_lightcurve(avgFWHM)

#avgFWHM /= const    
### Find maximum FWHM as this is what all the others willl become ###
aimind = np.argmax(avgFWHM)
#
### Convert FWHM into a sigma ###
sigmaold = np.array([FWHM2sigma(fwhm, const) for fwhm in avgFWHM])
sigmabroad = sigmaold[aimind]

### Find required sigma ###
# sigker^2 = sigbroad^2 - signar^2
sigmakernel = np.array([np.sqrt(sigmabroad**2 - sigma**2) for sigma in sigmaold])
#
#### Open images ###
#im05Bfull = fits.open('UDS_05B_K_bin2x2.fits', memmap=True)
#im05B = im05Bfull[0].data
#hdr5 = im05Bfull[0].header
#im08Bfull = fits.open('UDS_08B_K_bin2x2.fits', memmap=True)
#im08B = im08Bfull[0].data
#hdr8 = im08Bfull[0].header
#im09Bfull = fits.open('UDS_09B_K_bin2x2.fits', memmap=True)
#im09B = im09Bfull[0].data
#hdr9 = im09Bfull[0].header
#im10Bfull = fits.open('UDS_10B_K_bin2x2.fits', memmap=True)
#im10B = im10Bfull[0].data
#hdr10 = im10Bfull[0].header
#im11Bfull = fits.open('UDS_11B_K_bin2x2.fits', memmap=True)
#im11B = im11Bfull[0].data
#hdr11 = im11Bfull[0].header
#im12Bfull = fits.open('UDS_12B_K_bin2x2.fits', memmap=True)
#im12B = im12Bfull[0].data
#hdr12 = im12Bfull[0].header
#
#### Convolve Images - don't need to do 05B or 07B ###
#def convolveimage(im, sigma):
#    kernel = Gaussian2DKernel(sigma)
#    return convolve(im, kernel, normalize_kernel=True)
#print('Convolving 05B')
#newim05B = convolveimage(im05B, sigmakernel[0])
#print('Convolving 08B')
#newim08B = convolveimage(im08B, sigmakernel[2])
#print('Convolving 09B')
#newim09B = convolveimage(im09B, sigmakernel[3])
#print('Convolving 10B')
#newim10B = convolveimage(im10B, sigmakernel[4])
#print('Convolving 11B')
#newim11B = convolveimage(im11B, sigmakernel[5])
#print('Convolving 12B')
#newim12B = convolveimage(im12B, sigmakernel[6])
#
#### Save the file ###
#hdu = fits.PrimaryHDU(newim05B, header=hdr5)
#hdulist = fits.HDUList([hdu])
#hdulist.writeto('new_UDS_05B_K_bin2x2.fits', overwrite=True)
#hdu = fits.PrimaryHDU(newim08B, header=hdr8)
#hdulist = fits.HDUList([hdu])
#hdulist.writeto('new_UDS_08B_K_bin2x2.fits', overwrite=True)
#hdu = fits.PrimaryHDU(newim09B, header=hdr9)
#hdulist = fits.HDUList([hdu])
#hdulist.writeto('new_UDS_09B_K_bin2x2.fits', overwrite=True)
#hdu = fits.PrimaryHDU(newim10B, header=hdr10)
#hdulist = fits.HDUList([hdu])
#hdulist.writeto('new_UDS_10B_K_bin2x2.fits', overwrite=True)
#hdu = fits.PrimaryHDU(newim11B, header=hdr11)
#hdulist = fits.HDUList([hdu])
#hdulist.writeto('new_UDS_11B_K.fits_bin2x2', overwrite=True)
#hdu = fits.PrimaryHDU(newim12B, header=hdr12)
#hdulist = fits.HDUList([hdu])
#hdulist.writeto('new_UDS_12B_K.fits_bin2x2', overwrite=True)
#
#### CLose the file ###
#im05Bfull.close()
#del im05Bfull[0].data
#im08Bfull.close()
#del im08Bfull[0].data
#im09Bfull.close()
#del im09Bfull[0].data
#im10Bfull.close()
#del im10Bfull[0].data
#im11Bfull.close()
#del im11Bfull[0].data
#im12Bfull.close()
#del im12Bfull[0].data
















