#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 15:25:00 2019

Module containing functions for the variability statisitics used

@author: ppxee
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.stats import median_absolute_deviation

import flux_funcs #for binning/editting flux arrays

def sigmasq(flux, baseerr):
    ''' Function that calculates the excess varience value for each row in an 
    array 
    Inputs:
        flux = array of fluxes from objects in a number of epochs 
        baseerr = array of errors that the mean error should be calculated from
    Output:
        normsig = array of excess variance values for every object '''
    avgflux = np.nanmean(flux, axis=1)
#    meanerr = np.nanmean(baseerr, axis=1)
#    meanerrsq = np.square(meanerr)
    N = np.size(flux, axis=1)
    numobs = np.size(flux, axis=0)
    sig = [((flux[n, :]- avgflux[n])**2 - (baseerr[n,:])**2) for n in range(numobs)]# 
    sigsum = np.nansum(sig, axis=1)
    normsig = sigsum/(N)
    return normsig

def normsigmasq(flux, baseerr):
    ''' Function that calculates the excess varience value for each row in an 
    array 
    Inputs:
        flux = array of fluxes from objects in a number of epochs 
        baseerr = array of errors that the mean error should be calculated from
    Output:
        normsig = array of excess variance values for every object '''
    if np.shape(flux) == (0,8):
        return np.array([])
    avgflux = np.nanmean(flux, axis=1)
    N = np.size(flux, axis=1)
    numobs = np.size(flux, axis=0)
    sig = [((flux[n, :]- avgflux[n])**2 - (baseerr[n,:])**2) for n in range(numobs)]# 
    sigsum = np.nansum(sig, axis=1)
    normsig = sigsum/(N)
    return normsig

def maximum_likelihood_fig(testmag, testmagerr, meanmag, posvar):
    ''' Function that calculates the maximum likelihood variance for a single 
    source and plots the corresponding likelihood curce with error bars marked.
    
    Inputs:
        testmag = array of magnitudes/fluxes (i.e. the light curve) to test
        testmagerr = array of corresponding errors
        meanmag = mean of the lightcurve/test value for light curve being flat
        posvar = array of sigmas to test for maximum likelihood
        
    Outputs:
        sig = maximum likelihood sigma for the inputted light curve
        err = the error on sig according to the likelihood curve
    '''
    # Calculate likelihood curve
    L = np.array([np.nanprod((np.exp((-0.5*((testmag - meanmag)**2))/(
            testmagerr**2 + testsig**2)))/(((2*np.pi)**0.5)*
            (testmagerr**2 + testsig**2)**0.5)) for testsig in posvar])
    sig = float(posvar[L==np.nanmax(L)][0]) #sigma value at max L
    err = np.sqrt(np.average((posvar-np.average(posvar, weights=L))**2, weights=L))
    plt.figure()
    plt.plot(posvar, L)
    plt.vlines(sig, np.min(L), np.max(L))
    plt.vlines(sig+err, np.min(L), np.max(L))
    return sig, err

def maximum_likelihood(testmag, testmagerr, meanmag, posvar, n=None, printn=10):
    ''' Function that calculates the maximum likelihood variance for a single 
    source. Has the capability to print a counter as it progresses as is slow 
    to run over a full catalogue.
    
    Inputs:
        testmag = array of magnitudes/fluxes (i.e. the light curve) to test
        testmagerr = array of corresponding errors
        meanmag = mean of the lightcurve/test value for light curve being flat
        posvar = array of sigmas to test for maximum likelihood
        n = submitted counter for how many times the function has run. Default
            is None as do not usually want the counter
        printn = at what multiple of n to print a counter. Default is 10
    Outputs:
        sig = maximum likelihood sigma for the inputted light curve
        err = the error on sig according to the likelihood curve
    '''
    if n != None and n%printn == 0:
        print(n)
    # Calculate likelihood curve
    L = np.array([np.nanprod((np.exp((-0.5*((testmag - meanmag)**2))/(
            testmagerr**2 + testsig**2)))/(((2*np.pi)**0.5)*
            (testmagerr**2 + testsig**2)**0.5)) for testsig in posvar])
    sig = float(posvar[L==np.nanmax(L)][0]) #sigma value at max L
    if np.sum(L) == 0:
        return sig, np.nan
    else:
        err = np.sqrt(np.average((posvar-np.average(posvar, weights=L))**2, weights=L))
        return sig, err

def my_chisquare_err(flux, fluxerr):
    ''' Function that finds the chi square value of each light curve when passed
    an array of flux values and corresponding array of flux errors. It assumes
    each row is a single lightcurves and tests against the null hypothesis that 
    the light curve is flat with a value equal to the mean value of the actual 
    lightcurve.
    Inputs:
        flux = a 2D array of flux values where each row is a lightcurve for a 
               single object
        fluxerr = a 2D array of flux errors that correspond to the flux array
    Outputs:
        chi = a 1D array of chi sq values for each row in the flux arrays
    '''
#    flux, fluxerr = normalise_flux_and_errors(flux, fluxerr)
    meanflux = np.nanmean(flux, axis=1)
    top = np.square(flux-meanflux[:,None])
    bot = np.square(fluxerr)
    chi = np.nansum(top/bot, axis=1)
    return chi

def mod_z_score(arr):
    '''Function to find the modified z score of a given array, used to find
    variables in my first pass at this project
    Inputs:
        arr = array to find mod-z of
    Outputs:
        zvalues = array of z-values for that array
    '''
    medx = np.median(arr)
    mad = median_absolute_deviation(arr)
    zvalues = np.array([(0.6745*(x-medx))/mad for x in arr])
    return zvalues

def find_outliers(flux, tbdata, bins, threshold=6):
    '''Function used to find outliers when using MAD to select variables. It 
    splits the data into flux bins, calulates all the MAD values in that bin
    and then uses the modified z score to find which objects had disproportionately
    high MAD values for that bin. Anything above a given threshold was said to
    be an outlier.
    Inputs:
        flux = a 2D array of flux values where each row is a lightcurve for a 
               single object
        tbdata = original data table for those flux values
        bins = array of bin edges
        threshold = threshold at which everything with a mod-z score above that 
                    value is said to be an outlier. Default is 6
    Outputs:
        outliers = bool array defining which objects are outliers
        tbnew = table of data in the same order as outliers and allmodz so easy
                to compare/apply the boolean
        allmodz = array of z values for all the objects
    '''
    ### Bin data ###
    allmodz = []
    tbnew = np.recarray([0], dtype=tbdata.dtype, formats=tbdata.formats)
    for n, binedge in enumerate(bins):
        if n==np.size(bins)-1:
            break
        fluxbin1, tbbin1 = flux_funcs.fluxbin(binedge, bins[n+1], flux, tbdata)
        #calulate mad values in the bins
        mad1 = median_absolute_deviation(fluxbin1, axis=1) #for UDS
        modz = mod_z_score(mad1)
        tbnew = np.hstack([tbnew, tbbin1])
        allmodz = np.append(allmodz, modz)
    outliers = allmodz>=threshold
    return outliers, tbnew, allmodz
