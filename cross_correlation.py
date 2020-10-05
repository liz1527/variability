#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  3 12:15:56 2020

Module to contain the functions to run a cross correlation analysis on the 
J and K band lightcurves.

@author: ppxee
"""

### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
from astropy.table import Table #for handling tables
import numpy as np #for handling arrays
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
import random

def mean_subtract_normalise(flux):
    ''' Function to create mean subtracted and normalised lightcurves for use 
    in cross correlation 
    Inputs:
        flux = array of flux values where rows are light curves and columns are
               values for each epoch
    Outputs:
        newflux = array of mean subtracted, normalised light curves the same 
                  shape as flux.
    '''
    meanflux = np.nanmean(flux, axis=1)
    newflux = (flux - meanflux[:,None])/meanflux[:,None]
    return newflux
def weighted_mean_subtract_normalise(flux, fluxerr):
    ''' Function to create mean subtracted and normalised lightcurves for use 
    in cross correlation but including errors
    Inputs:
        flux = array of flux values where rows are light curves and columns are
               values for each epoch
        fluxerr = array of flux errors the same size/shape as flux
    Outputs:
        newflux = array of mean subtracted, normalised light curves the same 
                  shape as flux, where the errors were used to weight the 
                  average.
    '''
    meanflux = np.average(flux, axis=1, weights=1/(fluxerr**2))
    newflux = (flux - meanflux[:,None])/meanflux[:,None]
    return newflux

def weighted_mean_subtract_normalise_errs(flux, fluxerr):
    ''' Function to create mean subtracted and normalised lightcurves for use 
    in cross correlation but including errors
    Inputs:
        flux = array of flux values where rows are light curves and columns are
               values for each epoch
        fluxerr = array of flux errors the same size/shape as flux
    Outputs:
        newflux = array of mean subtracted, normalised light curves the same 
                  shape as flux, where the errors were used to weight the 
                  average.
        newfluxerr = array of corresponding flux errors the same size/shape as 
        flux
    '''
    meanflux = np.average(flux, axis=1, weights=1/(fluxerr**2))
    meanfluxerr = np.sqrt(1/np.sum(1/(fluxerr**2), axis=1))
    newflux1 = (flux - meanflux[:,None])
    newflux1err = fluxerr+meanfluxerr[:,None]
    newflux = newflux1/meanflux[:,None]
    newfluxerr = newflux * np.sqrt((newflux1err/newflux1)**2 + (meanfluxerr[:,None]/meanflux[:,None])**2)
    return newflux, newfluxerr

def weighted_mean_subtract_normalise_single(flux, fluxerr):
    ''' Function to create mean subtracted and normalised lightcurve for use 
    in cross correlation but including errors
    Inputs:
        flux = 1D array of flux values where columns are
               values for each epoch
        fluxerr = array of flux errors the same size/shape as flux
    Outputs:
        newflux = array of mean subtracted, normalised light curves the same 
                  shape as flux, where the errors were used to weight the 
                  average.
    '''
    meanflux = np.average(flux, weights=1/(fluxerr**2))
    newflux = (flux - meanflux)/meanflux
    return newflux

def make_corr_arrays(flux, mask, full_months):
    ''' Function to create arrays for the cross correlation, where each column
    is a month and the months with no flux values are set to np.nan. This means
    the epochs are spaced correctly in time.
    Inputs:
        flux = array of mean subtracted, normalised flux values where rows are 
               light curves and columns are values for each epoch
        mask = array that indicates which months there are data for, should be 
               the same length as full_months
        full_months = array containing all the month names in the range where 
                      the data was taken (e.g. sep05 -> jan13)
    Outputs:
        newflux = array of light curves with the correct time spacing. Should 
                  have the same number of rows as flux and columns as the length
                  of full_months
    '''
    newflux = np.empty([len(flux), len(full_months)]) # create empty array the right shape
    newflux[:,mask] = flux
    newflux[:,~mask] = np.nan
    return newflux

def make_corr_arrays_single(flux, mask, full_months):
    ''' Function to create arrays for the cross correlation, where each column
    is a month and the months with no flux values are set to np.nan. This means
    the epochs are spaced correctly in time.
    Inputs:
        flux = array of mean subtracted, normalised flux values where columns 
        are values for each epoch
        mask = array that indicates which months there are data for, should be 
               the same length as full_months
        full_months = array containing all the month names in the range where 
                      the data was taken (e.g. sep05 -> jan13)
    Outputs:
        newflux = array of light curves with the correct time spacing. Should 
                  have the same number of columns as the length
                  of full_months
    '''
    newflux = np.empty([1,len(full_months)]) # create empty array the right shape
    newflux[:,mask] = flux
    newflux[:,~mask] = np.nan
    return newflux

def cross_correlate(kflux, jflux, tau, type='ccf'):
    ''' Function to run the stacked cross-correlation analysis at a specfic lag 
    tau between the J and K band lightcurves given.
    Inputs:
        kflux = array of k-band light curves with the correct time spacing. 
        jflux = array of j-band light curves with the correct time spacing. 
        tau = the lag to caluclate the ccf value at
        type = 'ccf' or 'dcf' sepcifies what type of cross-correlation to run.
               default is 'ccf'
    Outputs:
        ccf = the cross-correlation value for this tau 
        err = the bootstrapped error on this ccf value
        len(t_k) = the total number of pairs used to calculate this value.
    '''
    ### set up time arrays from J and K ###
    _, monthnum = np.shape(kflux)
    t_k = np.arange(monthnum)
    t_j = t_k - tau #get array of inds that are tau from k inds
    
    ### limit to values with valid t_j values ###
    t_k = t_k[t_j<monthnum] # t_j should not be larger than total num of months
    t_j = t_j[t_j<monthnum]
    t_k = t_k[t_j>=0] # t_j should not be negative
    t_j = t_j[t_j>=0]
    
    ### limit to values where there are not nans in J or K at those indicies ###
    for k, j in zip(t_k, t_j):
        if np.isnan(kflux[0,k]) or np.isnan(jflux[0,j]): # can test on one row as all the same
            ### set index of both to 999 if there is no flux in one position ###
            t_k[t_k == k] = 999
            t_j[t_j == j] = 999
    
    t_k = t_k[t_k!=999]
    t_j = t_j[t_j!=999]
    
#    ### check lags are correct ###
#    checktau = t_k - t_j
#    print('Tau = '+str(tau))
#    print(checktau)
#    
    ### limit flux arrays to just these months ###
    pair_k_flux = kflux[:,t_k]
    pair_j_flux = jflux[:,t_j]
    
    if type == 'ccf':
        ccf = calculate_ccf(pair_k_flux, pair_j_flux)
    elif type == 'dcf':
        ccf = calculate_dcf(pair_k_flux, pair_j_flux)
        
    err = bootstrapping(pair_k_flux, pair_j_flux, type=type)
    
    return ccf, err, len(t_k)

def cross_correlate_de_z(kflux, jflux, tau, z, type='dcf'):
    ''' Function to run the stacked cross-correlation analysis at a specfic lag 
    tau between the J and K band lightcurves given when the light curves have 
    been de-redshifted so you cannot assume all the rows have the same light 
    curve spacings. In this case the ccf is done in bins which run from tau-0.5
    to tau+0.5 so that tau is in the centre of the bin
    Inputs:
        kflux = array of k-band light curves with the correct time spacing. 
        jflux = array of j-band light curves with the correct time spacing. 
        tau = the lag to caluclate the ccf value at
        z = redshifts of the objects
        type = 'ccf' or 'dcf' sepcifies what type of cross-correlation to run.
               default is 'dcf'
    Outputs:
        ccf = the cross-correlation value for this tau 
        err = the bootstrapped error on this ccf value
        len(t_k) = the total number of pairs used to calculate this value.
    '''
    ### limit to values where there are not nans in J or K at indicies in month bins ###
    pair_k_flux = np.array([])
    pair_j_flux = np.array([])
    for n in range(len(kflux)):
        ### set up time arrays from J and K ###
        _, monthnum = np.shape(kflux)
        full_t = np.arange(monthnum)
        t_k = full_t/(1+z[n]) # month times of j and k de redshifted
        t_arr = full_t/(1+z[n]) # need array that doesn't change for indexing
        
        for k in t_k:
            ### Set t_j ###
            t_j = full_t/(1+z[n])
            
            ### find range of t_j values that are valid for current tau bin ###
            mask1 = t_j > k - (tau+0.5) # t_j is within tau + 0.5 months of k
            mask2 = t_j <= k - (tau-0.5) # t_j is not more than tau-0.5 months from k
            mask = mask1*mask2.astype(bool)
            t_j_match = t_j[mask] # get values that match for loop
            t_j[~mask] = 999 # set non matches to null
            
            ### check if months that match have values ###
            for j in t_j_match:
                if np.isnan(kflux[n,t_arr==k]) or np.isnan(jflux[n,t_arr==j]):
                    ### set j index to null if no value ###
                    t_j[t_j==j] = 999 
            
            ### restrict to values with non-null indices ###
            t_j = t_j[t_j!=999]
            
            ### if there are no matches, set k to null and move on ###
            if len(t_j) == 0:
                t_k[t_k==k] = 999 #this probably isn't actually necessary in this version
                continue
        
            ### if there are matches add t_k to the array with each paired t_j ###
            for j in t_j:
                pair_k_flux = np.append(pair_k_flux, kflux[n,t_k==k])
                pair_j_flux = np.append(pair_j_flux, jflux[n,t_arr==j])

    if type == 'ccf':
        ccf = calculate_ccf(pair_k_flux, pair_j_flux)
    elif type == 'dcf':
        ccf = calculate_dcf(pair_k_flux, pair_j_flux)
        
    err = bootstrapping(pair_k_flux, pair_j_flux, repeats=1000)
#    err=np.nan
    return ccf, err, len(pair_k_flux)

def cross_correlate_shifted(kflux, jflux, tau, type='ccf'):
    ''' Function to run the stacked cross-correlation analysis at a specfic lag 
    tau between the J and K band lightcurves given when some of the light curves
    have been shifted so you cannot assume all the rows have the same light 
    curve spacings.
    Inputs:
        kflux = array of k-band light curves with the correct time spacing. 
        jflux = array of j-band light curves with the correct time spacing. 
        tau = the lag to caluclate the ccf value at
        type = 'ccf' or 'dcf' sepcifies what type of cross-correlation to run.
               default is 'ccf'
    Outputs:
        ccf = the cross-correlation value for this tau 
        err = the bootstrapped error on this ccf value
        len(t_k) = the total number of pairs used to calculate this value.
    '''
    ### set up time arrays from J and K ###
    _, monthnum = np.shape(kflux)
    t_k = np.arange(monthnum)
    t_j = t_k - tau #get array of inds that are tau from k inds
    
    ### limit to values with valid t_j values ###
    t_k = t_k[t_j<monthnum] # t_j should not be larger than total num of months
    t_j = t_j[t_j<monthnum]
    t_k = t_k[t_j>=0] # t_j should not be negative
    t_j = t_j[t_j>=0]
    
    ### limit to values where there are not nans in J or K at those indicies ###
    pair_k_flux = np.array([])
    pair_j_flux = np.array([])
    for n in range(len(kflux)):
        for k, j in zip(t_k, t_j):
            if np.isnan(kflux[n,k]) or np.isnan(jflux[n,j]): 
                ### set index of both to 999 if there is no flux in one position ###
                t_k[t_k == k] = 999
                t_j[t_j == j] = 999
        
        t_k = t_k[t_k!=999]
        t_j = t_j[t_j!=999]
    
        ### limit flux arrays to just these months ###
        pair_k_flux = np.append(pair_k_flux, kflux[n,t_k])
        pair_j_flux = np.append(pair_j_flux, jflux[n,t_j]) 
        

    if type == 'ccf':
        ccf = calculate_ccf(pair_k_flux, pair_j_flux)
    elif type == 'dcf':
        ccf = calculate_dcf(pair_k_flux, pair_j_flux)#    err = bootstrapping(pair_k_flux, pair_j_flux, repeats=1000)
    err=np.nan
    return ccf, err, len(pair_k_flux)

def calculate_ccf(pair_k_flux, pair_j_flux):
    ''' Function to calculate the ccf value for the paired flux arrays given. 
    Needs to be a separate function so that I can call it from the bootstrapping
    function as well as the cross-correlation function.
    Inputs:
        pair_k_flux = 1D array of k-band flux values 
        pair_j_flux = paired 1D array of j-band flux values 
    Outputs:
        ccf = the output of a stacked cross-correlation of these pairs of flux
              values
    '''
    ### multiply arrays ###
    k_flux_test_j_flux = pair_k_flux * pair_j_flux
    
    ### find the sum of this multiplied array and divide by tot number of pairs ###
    multi_sum = np.sum(k_flux_test_j_flux)
    tot_num_pairs = np.size(k_flux_test_j_flux)
    ccf = multi_sum/tot_num_pairs
    
    return ccf

def calculate_dcf(pair_k_flux, pair_j_flux):
    ''' Function to calculate the dcf value for the paired flux arrays given. 
    The dcf is the version used in papers, as oppose to the ccf that Mike sent.
    Needs to be a separate function so that I can call it from the bootstrapping
    function as well as the cross-correlation function.
    Inputs:
        pair_k_flux = 1D array of k-band flux values 
        pair_j_flux = paired 1D array of j-band flux values 
    Outputs:
        ccf = the output of a stacked cross-correlation of these pairs of flux
              values
    '''
    ### multiply arrays ###
    udcf = ((pair_k_flux - np.mean(pair_k_flux)) * (pair_j_flux - np.mean(pair_j_flux)))/(np.std(pair_k_flux)*np.std(pair_j_flux))
    
    ### find the sum of this multiplied array and divide by tot number of pairs ###
    multi_sum = np.sum(udcf)
    tot_num_pairs = np.size(udcf)
    ccf = multi_sum/tot_num_pairs
    
    return ccf    
    
def bootstrapping(pair_k_flux, pair_j_flux, repeats=5000, type='ccf'):
    ''' Function to calculate the bootstrapped error on the ccf value given by
    the paired flux arrays
    Inputs:
        pair_k_flux = 1D array of k-band flux values 
        pair_j_flux = paired 1D array of j-band flux values 
        repeats = number of resampling repeats to do during bootstrapping, 
                  default is 5000
        type = 'ccf' or 'dcf' sepcifies what type of cross-correlation to run.
               default is 'ccf'
    Outputs:
        sig = the standard deviation of the ccf values found through bootstapping,
              this is the standard error on the ccf value.
    '''
    ccf = np.zeros(repeats) #create array to store ccf values in
    n=1 #start counter
    while n<repeats: # repeats until specified number reached
        ### resample data ###
        inds = random.choices(range(len(pair_k_flux)), k=len(pair_k_flux)) # get random inds with same length
        
        ### run ccf ###
        if type == 'ccf':
            ccf[n] = calculate_ccf(pair_k_flux[inds], pair_j_flux[inds])
        elif type == 'dcf':
            ccf[n] = calculate_dcf(pair_k_flux[inds], pair_j_flux[inds])
        
        ### increase counter ##
        n+=1
    
    ### calculate sigma on ccf values ###
    sig = np.std(ccf)
    return sig