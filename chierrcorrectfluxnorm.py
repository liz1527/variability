#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 18 09:46:30 2018

OLD CODE - NOT UPDATED WITH RESTRUCTURING 

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

def add_chi_err(mag, magerr):
    testerrchange = np.logspace(1,4,1000) #create array of values to change err by
    allmed = np.array([]) #empty array to contain median chi-square values
    n=0 #counter
    for errchange in testerrchange:
        newmagerr = np.sqrt(np.square(magerr) + np.square(errchange)) #change error
        chisq = min_chi_sq(mag, newmagerr) #calulate chisquare values for all object in bin
        med = np.nanmedian(chisq) # find the median of these values
        # need to find the err change that gives a median value closest to known chi-square value
        if med < 6.346: # 6.346 is 0.5 chi square value for 7 dof and med 
                        # values decrease so stop when its less than this
            # need to check whether the value just above or just below is closer
            possmed = np.array([med,allmed[-1]]) 
            diff = np.abs(possmed - 6.35)
            if diff[0]<diff[1]:
                sigma = errchange
            else:
                sigma = testerrchange[n-1]
            break #quit loop as no point continuing after value is found
        allmed = np.append(allmed, med) #keep value on array
        n=n+1 #increase counter
    # calculate new error array with correct correction
    newmagerr = np.sqrt(np.square(magerr) + np.square(sigma))
    return sigma, newmagerr #output the required err change and the new array

def times_chi_err(mag, magerr):
    #same as above but multiply change rather than add in quadriture
    testerrchange = np.linspace(1,50,5000)
    allmed2 = np.array([])
    n=0
    for errchange in testerrchange:
        newmagerr2 = errchange*magerr
        chisq = min_chi_sq(mag, newmagerr2)
        med = np.nanmedian(chisq)
        if med < 6.35:
            possmed = np.array([med,allmed2[-1]])
            diff = np.abs(possmed - 6.35)
            if diff[0]<diff[1]:
                sigma2 = errchange
            else:
                sigma2 = testerrchange[n-1]
            break
        allmed2 = np.append(allmed2, med)
        n=n+1
    newmagerr2 = sigma2*magerr
    return sigma2, newmagerr2

def add_chi_err_norm(mag, magerr): 
    testerrchange = np.linspace(0,1,1000) #create array of values to change err by
    allmed = np.array([]) #empty array to contain median chi-square values
    n=0 #counter
    for errchange in testerrchange:
        newmagerr = np.sqrt(np.square(magerr) + np.square(errchange)) #change error
        chisq = min_chi_sq(mag, newmagerr) #calulate chisquare values for all object in bin
        med = np.nanmedian(chisq) # find the median of these values
        # need to find the err change that gives a median value closest to known chi-square value
        if med < 6.346: # 6.346 is 0.5 chi square value for 7 dof and med 
                        # values decrease so stop when its less than this
            # need to check whether the value just above or just below is closer
            possmed = np.array([med,allmed[-1]]) 
            diff = np.abs(possmed - 6.35)
            if diff[0]<diff[1]:
                sigma = errchange
            else:
                sigma = testerrchange[n-1]
            break #quit loop as no point continuing after value is found
        allmed = np.append(allmed, med) #keep value on array
        n=n+1 #increase counter
    # calculate new error array with correct correction
    newmagerr = np.sqrt(np.square(magerr) + np.square(sigma))
    return sigma, newmagerr #output the required err change and the new array

### Open the fits files and get data ###
tbdata = fits.open('mag_flux_tables/mag_flux_table_best.fits')[1].data
chandata = fits.open('mag_flux_tables/xray_mag_flux_table_best.fits')[1].data
sdata = fits.open('mag_flux_tables/stars_mag_flux_table.fits')[1].data

### extract magnitude arrays ###
allmag = vari_funcs.flux5_stacks(tbdata)
allchanmag = vari_funcs.flux5_stacks(chandata) # for chandra non-stellar objects
allsmag = vari_funcs.flux5_stacks(sdata)

### remove 99s ###
allmag, tbdata = vari_funcs.noneg(allmag, tbdata)
allchanmag, chandata = vari_funcs.noneg(allchanmag, chandata)
allsmag, sdata = vari_funcs.noneg(allsmag, sdata)

# need to run on all fluxbins, like in find_outliers
#bins = np.array([13, 15])
#bins = np.append(bins, np.arange(16,24,0.2))
#bins = np.append(bins, [24, 25, 26])
bins = np.array([13, 15])
bins = np.append(bins, np.arange(16,24,0.2))
bins = np.append(bins, [24])

bins = 10**((30-bins)/2.5)
bins = np.flip(bins, axis=0)

#mag, bindata = vari_funcs.fluxbin(bins[25], bins[26], allmag, tbdata)
#magerr = vari_funcs.magerr5_stacks(bindata)
#errchange, newmagerr = add_chi_err(mag, magerr)

### Bin data ###
allerrchange = np.array([])
newallmag = np.array([])
newallmagerr = np.array([])
newallchanmag = np.array([])
newallchanmagerr = np.array([])
newallsmag = np.array([])
newallsmagerr = np.array([])
for n, binedge in enumerate(bins):
    print(binedge)
    if n==np.size(bins)-1:
        break
    mag, bindata = vari_funcs.fluxbin(binedge, bins[n+1], allmag, tbdata) #bindata
    magerr = vari_funcs.fluxerr5_stacks(bindata) #make error array
    
    meanflux = np.mean(mag, axis=1)
    mag, magerr = vari_funcs.normalise_flux_and_errors(mag, magerr)
    
    errchange, newmagerr = times_chi_err(mag, magerr) #find correction
    allerrchange = np.append(allerrchange, errchange) #create array of corrections
    
    ### Apply correction tp stars and X-ray ###
    chanmag, chanbin = vari_funcs.fluxbin(binedge, bins[n+1], allchanmag, chandata)
    chanmagerr = vari_funcs.fluxerr5_stacks(chanbin)
    chanmeanmag = np.mean(chanmag, axis=1)
    chanmag, chanmagerr = vari_funcs.normalise_flux_and_errors(chanmag, chanmagerr)
    newchanmagerr = chanmagerr*errchange #np.sqrt(np.square(chanmagerr) + np.square(errchange)) #change error
    smag, sbin = vari_funcs.fluxbin(binedge, bins[n+1], allsmag, sdata)
    smagerr = vari_funcs.fluxerr5_stacks(sbin)
    smeanmag = np.mean(smag, axis=1)
    smag, smagerr = vari_funcs.normalise_flux_and_errors(smag, smagerr)
    newsmagerr = smagerr*errchange #np.sqrt(np.square(smagerr) + np.square(errchange)) #change error
    
    #create new arrays
    if n == 0:
        newallmag = mag
        meanmagnew = meanflux
        oldallmagerr = magerr
        newallmagerr = newmagerr
        newallchanmag = chanmag
        newallchanmagerr = newchanmagerr
        meanchanmagnew = chanmeanmag
        newallsmag = smag
        newallsmagerr = newsmagerr
        meansmagnew = smeanmag
    else:
        newallmag = np.vstack([newallmag, mag])
        meanmagnew = np.append(meanmagnew, meanflux)
        oldallmagerr = np.vstack([oldallmagerr, magerr])
        newallmagerr = np.vstack([newallmagerr, newmagerr])
        newallchanmag = np.vstack([newallchanmag, chanmag])
        newallchanmagerr = np.vstack([newallchanmagerr, newchanmagerr])
        meanchanmagnew = np.append(meanchanmagnew, chanmeanmag)
        newallsmag = np.vstack([newallsmag, smag])
        newallsmagerr = np.vstack([newallsmagerr, newsmagerr])
        meansmagnew = np.append(meansmagnew, smeanmag)

#%% plot excess variance with new arrays

excess = vari_funcs.sigmasq(newallmag, newallmagerr)
chanexcess = vari_funcs.sigmasq(newallchanmag, newallchanmagerr)
sexcess = vari_funcs.sigmasq(newallsmag, newallsmagerr)

plt.figure()
plt.plot(meansmagnew, sexcess, 'm*', markersize=10, mfc='None')
plt.plot(meanmagnew, excess, 'b+')
plt.plot(meanchanmagnew, chanexcess, 'ro', markersize=10, mfc='None')
plt.yscale('log')
plt.xscale('log')
plt.xlabel('Mean Flux')
plt.ylabel('Excess Variance')
plt.title('With Normalised Multiplicative Correction')

#%% plot errchange with magnitude
plt.figure()
plt.scatter(bins[0:42], allerrchange)
plt.xscale('log')
plt.xlabel('Flux')
plt.ylabel('Multiplicative Correction')

chisq = min_chi_sq(newallmag, newallmagerr)
chanchisq = min_chi_sq(newallchanmag, newallchanmagerr)
schisq = min_chi_sq(newallsmag, newallsmagerr)

plt.figure()
plt.plot(meansmagnew, schisq, 'm*', markersize=10, mfc='None')
plt.plot(meanmagnew, chisq, 'b+')
plt.plot(meanchanmagnew, chanchisq, 'ro', markersize=10, mfc='None')
plt.yscale('log')
plt.xscale('log')
plt.xlabel('Mean Flux')
plt.ylabel('Chi squared value')
plt.title('With Normalised Multiplicative Correction')

##%% Do with times instead of add
#### Bin data ###
#allerrchange = np.array([])
#newallmag = np.array([])
#newallmagerr = np.array([])
#newallchanmag = np.array([])
#newallchanmagerr = np.array([])
#newallsmag = np.array([])
#newallsmagerr = np.array([])
#for n, binedge in enumerate(bins):
#    print(binedge)
#    if n==np.size(bins)-1:
#        break
#    mag, bindata = vari_funcs.fluxbin(binedge, bins[n+1], allmag, tbdata) #bindata
#    magerr = vari_funcs.fluxerr5_stacks(bindata) #make error array
#    errchange, newmagerr = times_chi_err(mag, magerr) #find correction
#    allerrchange = np.append(allerrchange, errchange) #create array of corrections
#    
#    ### Apply correction tp stars and X-ray ###
#    chanmag, chanbin = vari_funcs.fluxbin(binedge, bins[n+1], allchanmag, chandata)
#    chanmagerr = vari_funcs.fluxerr5_stacks(chanbin)
#    newchanmagerr = np.sqrt(np.square(chanmagerr) + np.square(errchange)) #change error
#    smag, sbin = vari_funcs.fluxbin(binedge, bins[n+1], allsmag, sdata)
#    smagerr = vari_funcs.fluxerr5_stacks(sbin)
#    newsmagerr = np.sqrt(np.square(smagerr) + np.square(errchange)) #change error
#    
#    #create new arrays
#    if n == 0:
#        newallmag = mag
#        oldallmagerr = magerr
#        newallmagerr = newmagerr
#        newallchanmag = chanmag
#        newallchanmagerr = newchanmagerr
#        newallsmag = smag
#        newallsmagerr = newsmagerr
#    else:
#        newallmag = np.vstack([newallmag, mag])
#        oldallmagerr = np.vstack([oldallmagerr, magerr])
#        newallmagerr = np.vstack([newallmagerr, newmagerr])
#        newallchanmag = np.vstack([newallchanmag, chanmag])
#        newallchanmagerr = np.vstack([newallchanmagerr, newchanmagerr])
#        newallsmag = np.vstack([newallsmag, smag])
#        newallsmagerr = np.vstack([newallsmagerr, newsmagerr])
#
##%% plot excess variance with new arrays
##_, excess = vari_funcs.flux_variability_plot(newallmag, newallchanmag, 'excess',
##                                       fluxerr = newallmagerr, 
##                                       starfluxerr = newallsmagerr,
##                                       starflux=newallsmag, stars=True, 
##                                       chanerr = newallchanmagerr,
##                                       normalised=True)
#_, excess2 = vari_funcs.flux_variability_plot(newallmag, newallchanmag, 'excess',
#                                       fluxerr = newallmagerr, 
#                                       starfluxerr = newallsmagerr,
#                                       starflux=newallsmag, stars=True, 
#                                       chanerr = newallchanmagerr)
##%% plot errchange with magnitude
#plt.figure()
#plt.scatter(bins[0:42], allerrchange)
#plt.xscale('log')
#plt.xlabel('Flux')
#plt.ylabel('Multiplicative Error Change')
#
##%% Plot new chisquared graph
#meanmagnew = np.mean(newallmag, axis=1)
#meanchanmagnew = np.mean(newallchanmag, axis=1)
#meansmagnew = np.mean(newallsmag, axis=1)
#chisq = min_chi_sq(newallmag, newallmagerr)
#chanchisq = min_chi_sq(newallchanmag, newallchanmagerr)
#schisq = min_chi_sq(newallsmag, newallsmagerr)
#
#plt.figure()
#plt.plot(meansmagnew, schisq, 'm*', markersize=10, mfc='None')
#plt.plot(meanmagnew, chisq, 'b+')
#plt.plot(meanchanmagnew, chanchisq, 'ro', markersize=10, mfc='None')
#plt.yscale('log')
#plt.xscale('log')
#plt.xlabel('Mean Flux')
#plt.ylabel('Chi squared value')
#plt.title('With Multiplicative Correction')