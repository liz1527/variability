#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 16:31:17 2019

Module for functions that act on the flux or magnitude light curve arrays

@author: ppxee
"""
import numpy as np

def normalise_flux(flux):
    ''' Normalise each objects flux so its average flux is 1
    Input:
        flux = array of object flux values 
    Output:
        array of object flux values normalised to an average flux of 1
    '''
    avgflux = np.mean(flux, axis=1)
    return flux / avgflux[:,None]

def normalise_flux_and_errors(flux, fluxerr):
    ''' Normalise each objects flux so its average flux is 1
    Input:
        flux = array of object flux values 
        fluxerr = array of flux errors
    Output:
        flux = array of object flux values normalised to an average flux of 1
        fluxerr = array of flux errors scales to match
    '''
    avgflux = np.mean(flux, axis=1)
    flux = flux/avgflux[:,None]
    fluxerr = fluxerr/avgflux[:,None]
    return flux, fluxerr

def normalise_mag(mag):
    ''' Normalise each objects mag so its average magnitude is 1
    Input:
        mag = array of object magnitude values 
    Output:
        array of object magnitude values normalised to an average mag of 1 
    
    *** NOT SURE THIS FUNCTION IS IN ANY WAY CORRECT ***    
    '''
    avgflux = np.mean(mag, axis=1)
    diff = avgflux - 1
    return mag - diff[:,None]

def no99(fluxn, tbdata):
    ''' Function to remove any objects that have a magnitude of 99 in any epoch
    of their light curve from the analysis
    
    Inputs:
        fluxn = array of magnitude values
        tbdata = original table of data (should be the same length as fluxn)
    
    Outputs:
        fluxn = new array of magnitude values where lines that contained 99s have 
                been removed
        tbdata = new table of data where any objects masked out in the magnitude 
                 array have also been removed here so that the arrays still 
                 correspond to eachother
    '''
    fluxn[fluxn == 99] = np.nan
    mask = ~np.isnan(fluxn).any(axis=1)
    fluxn = fluxn[mask]
    tbdata = tbdata[mask]
    return fluxn, tbdata

def noneg(fluxn, tbdata):
    ''' Function to remove any objects that have a negative flux in any epoch
    of their light curve from the analysis
    
    Inputs:
        fluxn = array of flux values
        tbdata = original table of data (should be the same length as fluxn)
    
    Outputs:
        fluxn = new array of flux values where lines that contained negatives  
                have been removed
        tbdata = new table of data where any objects masked out in the flux 
                 array have also been removed here so that they still 
                 correspond to eachother
    '''
    fluxn[fluxn <= 0] = np.nan
    mask = ~np.isnan(fluxn).any(axis=1)
    fluxn = fluxn[mask]
    tbdata = tbdata[mask]
    return fluxn, tbdata

def fluxlim(fluxn, tbdata, lim=3527):
    ''' Function that imposes a flux limit on the catalogue as defined by lim.
    any object with an average flux less than this value is removed from the 
    catalogue.
    
    Inputs:
        fluxn = array of flux values
        tbdata = original table of data (should be the same length as fluxn)
        lim = defines the flux limit to be imposed
    
    Outputs:
        fluxn = new array of flux values where objects below the flux limit  
                have been removed
        tbdata = new table of data where any objects masked out in the flux 
                 array have also been removed here so that they still 
                 correspond to eachother
    '''
    fluxn, tbdata = noneg(fluxn, tbdata) #remove negative values first
    avgflux = np.mean(fluxn, axis=1) #find average so only remove those whose average is less than limit
    fluxn[avgflux <= lim] = np.nan
    mask = ~np.isnan(fluxn).any(axis=1)
    fluxn = fluxn[mask]
    tbdata = tbdata[mask]
    return fluxn, tbdata

def semfluxlim(fluxn, tbdata):
    ''' Function that imposes a semester flux limit on the catalogue.
    Any object with an flux less than the detection depth in any semester
    is removed from the catalogue.
    
    Inputs:
        fluxn = array of flux values
        tbdata = original table of data (should be the same length as fluxn)
    
    Outputs:
        fluxn = new array of flux values where objects below the flux limit  
                have been removed
        tbdata = new table of data where any objects masked out in the flux 
                 array have also been removed here so that they still 
                 correspond to eachother
    '''
#    fluxn, tbdata = noneg(fluxn, tbdata) #remove negative values first
    limits = [2920.74, 6917.62, 4023.65, 5613.71, 
              1985.52, 2725.65, 2111.91, 1915.96]
    for n, lim in enumerate(limits):
        fluxn[:,n][fluxn[:,n] < lim] = np.nan        
    mask = ~np.isnan(fluxn).any(axis=1)
    fluxn = fluxn[mask]
    tbdata = tbdata[mask]
    return fluxn, tbdata

def fluxbin(min, max, flux, tbdata):
    ''' Separate a flux array into a bin depending on the average flux of the
    object 
    Inputs:
        min = the minimum flux value for the bin
        max = the maximum flux value for the bin
        flux = the array for fluxes for the object
        tbdata = original table of data
    Output:
        fluxbin = the array of fluxes for objects whose average flux is within 
                    the limits of the bin
        tbdata = table of data for that bin                
    '''
    avgflux = np.mean(flux, axis=1)
    bina = avgflux >= min
    binb = avgflux < max
    bin = bina*binb
    fluxbin = flux[bin]
    tbdata = tbdata[bin]
    return fluxbin, tbdata

def fluxbinerr(min, max, flux, fluxerr):
    ''' Separate a flux array into a bin depending on the average flux of the
    object with errors returned as well
    Inputs:
        min = the minimum flux value for the bin
        max = the maximum flux value for the bin
        flux = the array for fluxes for the object
        fluxerr = the array for flux errors for the object
    Output:
        fluxbin = the array of fluxes for objects whose average flux is within 
                    the limits of the bin 
        fluxerrbin = the array of flux errors for objects whose average flux is 
                    within the limits of the bin '''
    avgflux = np.mean(flux, axis=1)
    bina = avgflux >= min
    binb = avgflux < max
    bin = bina*binb
    fluxbin = flux[bin]
    fluxerrbin = fluxerr[bin]
    return fluxbin, fluxerrbin