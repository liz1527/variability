#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 15:25:00 2019

Module containing functions for the variability statisitics used

@author: ppxee
"""

import numpy as np
import matplotlib.pyplot as plt


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
