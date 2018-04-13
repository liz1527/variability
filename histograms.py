#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 14:32:45 2017

Code that create histrograms to see if certain catagories of sources vary more 
than others

@author: ppxee
"""


### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
#from astropy.table import Table #for handling tables
import numpy as np #for handling arrays
#import math
from scipy import stats
from astropy.stats import median_absolute_deviation
import vari_funcs_no06
plt.close('all') #close any open plots

### Open the fits files and get data ###
combined = fits.open('non_xray_mag_flux_table_best.fits')
tbdata = combined[1].data
chandra = fits.open('xray_mag_flux_table_best.fits')
chandata = chandra[1].data
#chanposstars = fits.open('chandra_possible_stars.fits')
#pstardata = chanposstars[1].data

### Restrict objects to those in the Chandra field ###
tbdata = vari_funcs_no06.chandra_only(tbdata)

### Create arrays of flux values ###

flux = vari_funcs_no06.mag5_stacks(tbdata) # for UDS objects
fluxchan = vari_funcs_no06.mag5_stacks(chandata) # for chandra non-stellar objects
#fluxpstar = vari_funcs_no06.mag5_stacks(pstardata) # for chandra objects identified as stars

### remove values that are +/-99 ###
flux, tbdata = vari_funcs_no06.no99(flux, tbdata)
fluxchan, chandata = vari_funcs_no06.no99(fluxchan, chandata)
#fluxpstar = vari_funcs_no06.no99(fluxpstar)

### Multiply all flux values in a yearstack by the correct constant ###
flux = vari_funcs_no06.psf_correct(flux, flux, 'median') 
fluxchan = vari_funcs_no06.psf_correct(flux, fluxchan, 'median') 
#fluxpstar = vari_funcs_no06.psf_correct(flux, fluxpstar, 'mean') 

#calculate average flux for each object
#avgflux = np.mean(flux, axis=1) #for UDS
#avgfluxchan = np.mean(fluxchan, axis=1) #for non-stellar chandra
#avgfluxpstar = np.mean(fluxpstar, axis=1) #for stellar chandra

### Define bins for histograms and average light curves ###
fluxbin1, tbdata1 = vari_funcs_no06.fluxbin(15.0, 17.5, flux, tbdata)
fluxbin2, tbdata2 = vari_funcs_no06.fluxbin(17.5, 20.0, flux, tbdata)
fluxbin3, tbdata3 = vari_funcs_no06.fluxbin(20.0, 22.5, flux, tbdata)
fluxbin4, tbdata4 = vari_funcs_no06.fluxbin(22.5, 25.0, flux, tbdata)
#fluxbin5 = vari_funcs_no06.fluxbin(1e4, 1e5, flux)
#fluxbin6 = vari_funcs_no06.fluxbin(1e6, 1e7, flux)

### Calculate mad of each bin ###
mad = median_absolute_deviation(flux, axis=1) #for UDS
mad1 = median_absolute_deviation(fluxbin1, axis=1) #for UDS
mad2 = median_absolute_deviation(fluxbin2, axis=1) #for UDS
mad3 = median_absolute_deviation(fluxbin3, axis=1) #for UDS
mad4 = median_absolute_deviation(fluxbin4, axis=1) #for UDS
#mad5 = median_absolute_deviation(fluxbin5, axis=1) #for UDS
#mad6 = median_absolute_deviation(fluxbin6, axis=1) #for UDS

### Define bins for histograms and average light curves chandra ###
fluxbinchan1, chandata1 = vari_funcs_no06.fluxbin(15.0, 17.5, fluxchan, chandata)
fluxbinchan2, chandata2 = vari_funcs_no06.fluxbin(17.5, 20.0, fluxchan, chandata)
fluxbinchan3, chandata3 = vari_funcs_no06.fluxbin(20.0, 22.5, fluxchan, chandata)
fluxbinchan4, chandata4 = vari_funcs_no06.fluxbin(22.5, 25.0, fluxchan, chandata)
#fluxbinchan5 = vari_funcs_no06.fluxbin(1e4, 1e5, fluxchan)
#fluxbinchan6 = vari_funcs_no06.fluxbin(1e6, 1e7, fluxchan)

### Calculate mad of each bin chandra ###
madchan = median_absolute_deviation(fluxchan, axis=1)
madchan1 = median_absolute_deviation(fluxbinchan1, axis=1)
madchan2 = median_absolute_deviation(fluxbinchan2, axis=1)
madchan3 = median_absolute_deviation(fluxbinchan3, axis=1)
madchan4 = median_absolute_deviation(fluxbinchan4, axis=1)
#madchan5 = median_absolute_deviation(fluxbinchan5, axis=1)
#madchan6 = median_absolute_deviation(fluxbinchan6, axis=1)

### Find outliers ###
#
#outliers = vari_funcs_no06.find_outliers(mad2, threshold=3.5)
#varys = fluxbin2[outliers]
#chanoutliers = vari_funcs_no06.find_outliers(madchan2, threshold=3.5)
#varymad = np.append(mad2[outliers], madchan2[chanoutliers])

#### Plot histograms of MAD; one for chandra, one for UDS ###
colors = ['red','blue','green']
bins = np.linspace(0,0.08)

_ = plt.hist(madchan1, bins, normed = True, 
                            color=colors[0], alpha=0.5, histtype = 'stepfilled', 
                            label = 'Chandra')
_ = plt.hist(mad1, bins, normed = True, 
                            color=colors[1], alpha=0.5, histtype = 'step', 
                            label = 'UDS')
plt.xlabel('MAD magnitude')
plt.ylabel('Number')
plt.ylim(0,130)
plt.title('Variability distribution of sources with magnitudes between 15.0 and 17.5')
plt.legend()
plt.figure()

bins = np.linspace(0, 0.13)

_ = plt.hist(madchan2, bins, normed = True, 
                            color=colors[0], alpha=0.5, histtype = 'stepfilled', 
                            label = 'Chandra')
_ = plt.hist(mad2, bins, normed = True, 
                            color=colors[1], alpha=0.5, histtype = 'step', 
                            label = 'UDS')
#_ = plt.hist(varymad, bins, normed = True, 
#                            color=colors[2], alpha=0.5, histtype = 'step', 
#                            label = 'UDS')
plt.xlabel('MAD magnitude')
plt.ylabel('Number')
plt.ylim(0,60)
plt.title('Variability distribution of sources with magnitudes between 17.5 and 20.0')
plt.legend()
plt.figure()

bins = np.linspace(0, 0.40)
_ = plt.hist(madchan3, bins, normed = True, 
                            color=colors[0], alpha=0.5, histtype = 'stepfilled', 
                            label = 'Chandra')
_ = plt.hist(mad3, bins, normed = True, 
                            color=colors[1], alpha=0.5, histtype = 'step', 
                            label = 'UDS')
plt.xlabel('MAD magnitude')
plt.ylabel('Number')
plt.ylim(0,16)
plt.title('Variability distribution of sources with fluxes between 20.0 and 22.5')
plt.legend()
#plt.figure()
#
#_ = plt.hist(madchan4, normed = True, 
#                            color=colors[0], alpha=0.5, histtype = 'stepfilled', 
#                            label = 'Chandra')
#_ = plt.hist(mad4, normed = True, 
#                            color=colors[1], alpha=0.5, histtype = 'step', 
#                            label = 'UDS')
#plt.xlabel('MAD Flux')
#plt.ylabel('Number')
#plt.title('Variability distribution of sources with fluxes between 22.5 and 25.0')
#plt.legend()


### Plot mean variability v flux bin ###
#meanmad = [np.mean(mad1), np.mean(mad2), np.mean(mad3), np.mean(mad4), 
#           np.mean(mad5), np.mean(mad6)]
#meanmad = np.array(meanmad)
#meanflux = [np.mean(fluxbin1), np.mean(fluxbin2), 
#           np.mean(fluxbin3), np.mean(fluxbin4), 
#           np.mean(fluxbin5), np.mean(fluxbin6)]
#meanflux = np.array(meanflux)
#
#meanmadchan = [np.mean(madchan1), np.mean(madchan2), np.mean(madchan3), 
#               np.mean(madchan4), np.mean(madchan5), np.mean(madchan6)]
#meanmadchan = np.array(meanmadchan)
#meanfluxchan = [np.mean(fluxbinchan1), np.mean(fluxbinchan2), 
#           np.mean(fluxbinchan3), np.mean(fluxbinchan4), 
#           np.mean(fluxbinchan5), np.mean(fluxbinchan6)]
#meanfluxchan = np.array(meanfluxchan)
#
#plt.plot(meanflux, meanmad, 'b+', label = 'UDS mean')
#plt.plot(meanfluxchan, meanmadchan, 'ro', label = 'Chandra mean')
#plt.xscale('log')
#plt.yscale('log')
#plt.xlabel('Mean Flux')
#plt.ylabel('MAD Flux')
#plt.legend()

### Do KS Test to check is distributions are different ###
[D2, p2] = stats.ks_2samp(mad2, madchan2)
[D3, p3] = stats.ks_2samp(mad3, madchan3)
#[D5, p5] = stats.ks_2samp(mad5, madchan5)









combined.close()
chandra.close()
#chanposstars.close()







