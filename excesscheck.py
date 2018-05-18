#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 15 14:34:11 2018

testing effect of chi square correction on the value of excess variance

@author: ppxee
"""

### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
#from astropy.table import Table #for handling tables
import numpy as np #for handling arrays
import vari_funcs #my module to help run code neatly
plt.close('all') #close any open plots
#from numpy.lib.recfunctions import append_fields

newmag = np.load('magorderedarray.npy')
newmagerr = np.load('newmagerrarray.npy')
newchanmag = np.load('chanmagorderedarray.npy')
newchanmagerr = np.load('channewmagerrarray.npy')
newsmag = np.load('smagorderedarray.npy')
newsmagerr = np.load('snewmagerrarray.npy')

### Open the fits files and get data ###
tbdata = fits.open('mag_flux_tables/mag_flux_table_best.fits')[1].data

### extract magnitude array ###
allmag = vari_funcs.mag5_stacks(tbdata)

### remove 99s ###
allmag, tbdata = vari_funcs.no99(allmag, tbdata)
#allmag = vari_funcs.normalise_mag(allmag)
### get error array ###
magerr = vari_funcs.magerr5_stacks(tbdata)

### get correction info ###
corrections = np.load('chisquarecorrections.npy')

obmag = allmag[5200,:]
oberr = magerr[5200,:]
meanmag = np.nanmean(obmag)
newoberr =  np.sqrt(np.square(oberr) + np.square(corrections[1,39]))

def normsigmasq(flux, baseerr):
    ''' Function that calculates the excess varience value for each row in an 
    array 
    Inputs:
        flux = array of fluxes from objects in a number of epochs 
        baseerr = array of errors that the mean error should be calculated from
    Output:
        sig = array of excess variance values for every object '''
    avgflux = np.nanmean(flux)
    meanerr = np.nanmean(baseerr)
    meanerrsq = np.square(meanerr)
    N = np.size(flux)
#    numobs = np.size(flux, axis=0)
    sig = ((flux- avgflux)**2 - meanerrsq)# 
    sigsum = np.nansum(sig)
    normsig = sigsum/N
    return normsig

def min_chi_sq(mag, magerr):
    avgmag = np.nanmedian(mag, axis=1) #use median mags as a start for the expected model
    # create array that checks model around the median value
    testexpect = np.tile(np.array([avgmag]).transpose(), [1,50])
    testchanges = np.linspace(-0.1, 0.1) 
    testexpect += testchanges
    # Find chi-squared values with all possible models
    chisq = np.array([(np.square(mag-testexpect[:,n,None]))/np.square(magerr) for n in range(50)])
    chisq = np.nansum(chisq, axis=2) #sum to complete chi squared calculation
    chisqmin = np.nanmin(chisq, axis=0) #find the minimum chi-squared value
    return chisqmin
obsigma1= normsigmasq(obmag, oberr)
obsigma2 = normsigmasq(obmag, newoberr)
#
#sigma1= min_chi_sq(allmag, magerr)
#sigma2 = min_chi_sq(newmag, newmagerr)
#chansigma2 = min_chi_sq(newchanmag, newchanmagerr)
#ssigma2 = min_chi_sq(newsmag, newsmagerr)
#sigma3= vari_funcs.normsigmasq(allmag, magerr)
#sigma4 = vari_funcs.normsigmasq(newmag, newmagerr)

#plt.hist(sigma1[~np.isnan(sigma1)],bins=1000)
#plt.hist(sigma2[~np.isnan(sigma2)],bins=1000)
#plt.figure()
#plt.hist(sigma3[~np.isnan(sigma3)],bins=100)
#plt.hist(sigma4[~np.isnan(sigma4)],bins=100)
#vari_funcs.avg_lightcurve(obmag, oberr)
#vari_funcs.avg_lightcurve(obmag, newoberr)
#meanmag= np.nanmean(allmag, axis=1)
#meanmagnew = np.nanmean(newmag, axis=1)
#meanchanmagnew = np.nanmean(newchanmag, axis=1)
#meansmagnew = np.nanmean(newsmag, axis=1)
#plt.scatter(meanmag, sigma1)
#plt.yscale('log')
#plt.xlabel('Mean Magnitude')
#plt.ylabel('Chi squared value')
#plt.figure()
#plt.plot(meansmagnew, ssigma2, 'm*', markersize=10, mfc='None')
#plt.plot(meanmagnew, sigma2, 'b+')
#plt.plot(meanchanmagnew, chansigma2, 'ro', markersize=10, mfc='None')
#plt.yscale('log')
#plt.xlabel('Mean Magnitude')
#plt.ylabel('Chi squared value')
#
#
##plt.figure()
##plt.scatter(meanmag[~np.isnan(sigma3)], sigma3[~np.isnan(sigma3)])
##plt.yscale('log')
##plt.ylim(ymin=np.nanmin(sigma3))
##plt.xlabel('Mean Magnitude')
##plt.ylabel('Excess Variance')
##plt.figure()
##plt.scatter(meanmagnew[~np.isnan(sigma4)], sigma4[~np.isnan(sigma4)])
##plt.yscale('log')
##plt.ylim(ymin=np.nanmin(sigma3))
##plt.xlabel('Mean Magnitude')
##plt.ylabel('Excess Variance')
#fig = vari_funcs.flux_variability_plot(newmag, newchanmag, 'excess',
#                                       fluxerr = newmagerr, 
#                                       starfluxerr = newsmagerr,
#                                       starflux=newsmag, stars=True, 
#                                       chanerr = newchanmagerr,
#                                       normalised=True)