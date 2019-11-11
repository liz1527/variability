#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 14:27:16 2019

Module containing functions used to correct for seeing changes and test whether
these change have worked

@author: ppxee
"""
import numpy as np

def psf_correct_flux(baseflux, initflux, avgtype):
    ''' Function that applies a fractional PSF correction to initflux, based on
    the constant required to correct an epochs average flux to the overall 
    average flux in the baseflux array.
    i.e. if doing a correction based on all objects for objects in the chandra
    data then baseflux = flux and initflux = fluxchan.
    If doing a correction based on the stellar fluxes for all objects in the 
    field then baseflux = sflux and initflux = flux. 
    Basetype must be defined as either 'all' or 'star' so that the mean value
    can be used for all objects and the median value can be used for stellar 
    objects 
    
    *** OUTDATED METHOD - FUNCTION KEPT FOR REFERENCE ***
    
    Inputs:
        baseflux = flux array to base the corrections on
        initflux = flux array to be corrected
        basetype = either 'mean' or 'median' which dictate if the mean or 
                    median value is used for the correction 
    Output:
        Flux array with values crudely corrected for differences in seeing
        (average flux should now be the same for each epoch). '''
    if avgtype == 'mean':
        avgfluxperepoch = np.mean(baseflux, axis=0)#for UDS
        avgflux = np.mean(baseflux)
        const = avgflux/avgfluxperepoch
    elif avgtype == 'median':
        avgfluxperepoch = np.median(baseflux, axis=0)#for UDS
        avgflux = np.median(baseflux)
        const = avgflux/avgfluxperepoch
    else:
        print('Invalid basetype')
        return
    return initflux * const[None,:]

def psf_correct_mag(basemag, initmag, avgtype):
    ''' Function that applies a fractional PSF correction to initflux, based on
    the constant required to correct an epochs average flux to the overall 
    average flux in the baseflux array.
    i.e. if doing a correction based on all objects for objects in the chandra
    data then baseflux = flux and initflux = fluxchan.
    If doing a correction based on the stellar fluxes for all objects in the 
    field then baseflux = sflux and initflux = flux. 
    Basetype must be defined as either 'all' or 'star' so that the mean value
    can be used for all objects and the median value can be used for stellar 
    objects 
    
    *** OUTDATED METHOD - FUNCTION KEPT FOR REFERENCE ***
    
    Inputs:
        baseflux = flux array to base the corrections on
        initflux = flux array to be corrected
        basetype = either 'mean' or 'median' which dictate if the mean or 
                    median value is used for the correction 
    Output:
        Flux array with values crudely corrected for differences in seeing
        (average flux should now be the same for each epoch). '''
    if avgtype == 'mean':
        avgfluxperepoch = np.mean(basemag, axis=0)#for UDS
        avgflux = np.mean(basemag)
        const = avgflux-avgfluxperepoch
    elif avgtype == 'median':
        avgfluxperepoch = np.median(basemag, axis=0)#for UDS
        avgflux = np.median(basemag)
        const = avgflux-avgfluxperepoch
    else:
        print('Invalid basetype')
        return
    return initmag + const[None,:]

def err_correct(flux, fluxerr, fluxnew):
    ''' Function that applies a correction to the array of error values that 
    matches the correction applied to the corresponding array of fluxes.
    
    *** OUTDATED METHOD - FUNCTION KEPT FOR REFERENCE ***
    
    Inputs:
        flux = initial flux array before any corrections were applied
        fluxerr = initial flux err array
        fluxcorr = array of fluxes that have been corrected
    Output:
        Flux error array with values crudely corrected '''
        
    return fluxnew * (fluxerr/flux)

def err_correct_flux(oldflux, fluxerr):
    ''' Function that applies a correction to the array of error values that 
    matches the correction applied to the corresponding array of fluxes.
    
    *** OUTDATED METHOD - FUNCTION KEPT FOR REFERENCE ***
    
    Inputs:
        flux = initial flux array before any corrections were applied
        fluxerr = initial flux err array
        fluxcorr = array of fluxes that have been corrected
    Output:
        Flux error array with values crudely corrected '''
    avgflux = np.mean(oldflux, axis=1)
    return fluxerr/avgflux[:,None]