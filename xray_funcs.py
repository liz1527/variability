#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 15:41:52 2019

Module containing functions that define X-ray properties of variables.

@author: ppxee
"""
import numpy as np
import vari_funcs
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
from astropy import constants
# Set up cosmology
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

def get_xray_L(tbdata, Xray=True, band='Full'):
    ''' Function to find the X-ray luminosity in a variety of flux bands. 
    Inputs:
        tbdata = data table for all the objects
        Xray = whether the object has an X-ray counterpart. If it does not 
                then the survey limit is used as an upperlimit. Default is True
        band = what band of flux should be used:
                - 'Full' 0.5-10keV (default)
                - 'Hard' 2-10keV 
                - 'Soft' 0.2-5keV
                - 'Uhrd' 5-10keV (not recommended)
    Outputs:
        flux = broad band flux in units of ergs/cm^2/s
        L = broad band luminosity in units of ergs/s
        z = redshift
    '''
    z = vari_funcs.get_z(tbdata)
#    z[z<0.5] = np.nan
#    z[z>4] = np.nan
    
    ### get luminosity distance ###
    DL = cosmo.luminosity_distance(z) # need to sort out units
    DL = DL.to(u.cm)
    
    if band=='Hard': 
        if Xray == True: # if it is an X-ray source, get flux from catalogue
            flux = tbdata['Hard_flux'] 
            flux# Units of erg cm**-2 s**-1
        else: # if it is non X-ray - use the upper limit
            flux = np.zeros(len(tbdata))
            flux += 6.5e-16 # Units of erg cm**-2 s**-1
    elif band=='Full': 
        if Xray == True:
            flux = tbdata['Full_flux'] # Units of erg cm**-2 s**-1
        else:
            flux = 4.4e-16 # Units of erg cm**-2 s**-1
    elif band=='Soft': 
        if Xray == True:
            flux = tbdata['Soft_flux'] # Units of erg cm**-2 s**-1
        else:
            flux = 1.4e-16 # Units of erg cm**-2 s**-1
    elif band=='Uhrd': 
        if Xray == True:
            flux = tbdata['Uhrd_flux'] # Units of erg cm**-2 s**-1
        else:
            flux = 9.2e-15 # Units of erg cm**-2 s**-1
            
    ### Add units ###
    flux = flux* (u.erg) * (u.cm)**-2 * (u.s)**-1 
    
    ### get luminosity ###
    L = flux*4*np.pi*(DL**2)
    
    return flux, L, z

def get_xray_L_mixedtable(tbdata, band='Full'):
    ''' Function to find the X-ray luminosity in a variety of flux bands. 
    Inputs:
        tbdata = data table for all the objects, must have x-ray boolean incl
        band = what band of flux should be used:
                - 'Full' 0.5-10keV (default)
                - 'Hard' 2-10keV 
                - 'Soft' 0.2-5keV
                - 'Uhrd' 5-10keV (not recommended)
    Outputs:
        flux = broad band flux in units of ergs/cm^2/s
        L = broad band luminosity in units of ergs/s
        z = redshift
    '''
    z = vari_funcs.get_z(tbdata)
#    z[z<0.5] = np.nan
#    z[z>4] = np.nan
    
    ### get luminosity distance ###
    DL = cosmo.luminosity_distance(z) # need to sort out units
    DL = DL.to(u.cm)
    
    flux = np.zeros(len(tbdata))
    for n in range(len(tbdata)):
        if band=='Hard': 
            if tbdata['X-ray'][n] == True: # if it is an X-ray source, get flux from catalogue
                flux[n] = tbdata['Hard_flux'] 
#                flux# Units of erg cm**-2 s**-1
            else: # if it is non X-ray - use the upper limit
                flux[n] = 6.5e-16 # Units of erg cm**-2 s**-1
        elif band=='Full': 
            if tbdata['X-ray'][n] == True:
                flux[n] = tbdata['Full_flux'] # Units of erg cm**-2 s**-1
            else:
                flux[n] = 4.4e-16 # Units of erg cm**-2 s**-1
        elif band=='Soft': 
            if tbdata['X-ray'][n] == True:
                flux[n] = tbdata['Soft_flux'] # Units of erg cm**-2 s**-1
            else:
                flux[n] = 1.4e-16 # Units of erg cm**-2 s**-1
        elif band=='Uhrd': 
            if tbdata['X-ray'][n] == True:
                flux = tbdata['Uhrd_flux'] # Units of erg cm**-2 s**-1
            else:
                flux[n] = 9.2e-15 # Units of erg cm**-2 s**-1
            
    ### Add units ###
    flux = flux* (u.erg) * (u.cm)**-2 * (u.s)**-1 
    
    ### get luminosity ###
    L = flux*4*np.pi*(DL**2)
    
    return flux, L, z

def get_xray_L_2(tbdata, Xray=True, band='Hard'):
    ''' Function to find the monochromatic X-ray luminosity at 2keV in a 
    variety of flux bands. This assumes a power law for flux density and then
    finds the constant in this equation. The flux density is then used to find
    the monochromatic flux density which is in turn used to calulate the 
    monochromatic luminosity.
    Inputs:
        tbdata = data table for all the objects
        X-ray = whether the object has an X-ray counterpart. If it does not 
                then the survey limit is used as an upperlimit. Default is True
        band = what band of flux should be used:
                - 'Hard' 2-10keV (default)
                - 'Full' 0.5-10keV
                - 'Soft' 0.2-5keV
                - 'Uhrd' 5-10keV (not recommended)
    Outputs:
        L_2 = monochromatic luminosity at 2 keV in units of ergs/s/Hz
        F_2 = monochromatic flux at 2 keV in units of ergs/keV/s/cm^2
        flux = broad band flux in units of ergs/cm^2/s
        z = redshift
    '''
    z = vari_funcs.get_z(tbdata)
    z[z<0.5] = np.nan
    z[z>4] = np.nan
    
    ### get luminosity distance ###
    DL = cosmo.luminosity_distance(z) # need to sort out units
    DL = DL.to(u.cm)
    
    if band=='Hard': 
        upplim = 10 ## set band limits in keV
        lowlim = 2
        if Xray == True: # if it is an X-ray source, get flux from catalogue
            flux = tbdata['Hard_flux'] 
            flux# Units of erg cm**-2 s**-1
        else: # if it is non X-ray - use the upper limit
            flux = np.zeros(len(tbdata))
            flux += 6.5e-16 # Units of erg cm**-2 s**-1
    elif band=='Full': 
        upplim = 10
        lowlim = 0.5
        if Xray == True:
            flux = tbdata['Full_flux'] # Units of erg cm**-2 s**-1
        else:
            flux = 4.4e-16 # Units of erg cm**-2 s**-1
    elif band=='Soft': 
        upplim = 2
        lowlim = 0.5
        if Xray == True:
            flux = tbdata['Soft_flux'] # Units of erg cm**-2 s**-1
        else:
            flux = 1.4e-16 # Units of erg cm**-2 s**-1
    elif band=='Uhrd': 
        upplim = 10
        lowlim = 5
        if Xray == True:
            flux = tbdata['Uhrd_flux'] # Units of erg cm**-2 s**-1
        else:
            flux = 9.2e-15 # Units of erg cm**-2 s**-1
            
    ### Add units ###
    flux = flux* (u.erg) * (u.cm)**-2 * (u.s)**-1 
    upplim = upplim * u.keV
    lowlim = lowlim * u.keV
    
    ### redshift limits ###
#    upplim = upplim/(1+z)
#    lowlim = lowlim/(1+z)
    
    ### get integrated flux density ###
    denom = ((upplim**(0.1))/(0.1)) - ((lowlim**(0.1))/(0.1))
    print(denom)
    
    ### use this and flux value to get the power law constant ###
    const = flux / denom
    
    ### calculate flux density ###
    nu = 2 * u.keV # 2kev is value to evaluate at
    F_2 = const * (nu**(-0.9))
    
    
    ### calculate luminosity density ###
    L_2 = 4 * np.pi * (DL**2) * F_2#const * (nu**-0.9)
    
    L_2 = L_2.to((u.erg) * (u.s)**-1 * (u.Hz)**-1, equivalencies=u.spectral())
        
    L_2[L_2==0] = np.nan
    
    return L_2, F_2, flux, z #L_2_w_Hz


def get_L_O(tbdata):
    ''' Function to get the monochromatic luminosity at 2500 A
    Inputs:
        tbdata = table of data for the objects
    Outputs:
        L_O = monochromatic luminosity at 2500A in units of ergs/s/Hz
    '''
    z = vari_funcs.get_z(tbdata)
    z[z<0.5] = np.nan
    z[z>4] = np.nan
    
    ### get luminosity distance ###
    DL = cosmo.luminosity_distance(z) # need to sort out units
    DL = DL.to(u.cm)
    
    ### get U band magnitude ###
    umag = tbdata['M_U_z_p']
    
    L_U = 10**((34.1-umag)/2.5)
    L_U = L_U * u.W * (u.Hz)**-1 # Units of W/Hz
    
    L_U = L_U.to((u.erg) * (u.s)**-1 * (u.Hz)**-1, equivalencies=u.spectral())
    
    ### Use frequency ratios to find Omag ###
    Ofreq = constants.c/(2500 * u.AA)
    Ufreq = constants.c/(3743 * u.AA)
    
    L_O = L_U * (Ufreq/Ofreq)     # assuming flux propotional to freq**-1
#    F_O = F_O.to((u.erg) * (u.s)**-1 * (u.cm)**2 * (u.Hz)**-1, equivalencies=u.spectral()) # units of ergs/s/cm**"/Hz
    
    
#    L_O = L_O.to((u.erg) * (u.s)**-1 * (u.Hz)**-1, equivalencies=u.spectral())
    
    return L_O


def calc_alpha_Ox(L_O, L_2, xband=2*u.keV, optband=2500*u.AA):
    ''' Function to calculate the alpha_ox value of any given objects
    Inputs:
        L_O = monochromatic luminosity in optical band defined in optband
        L_2 = monochromatic luminosity in xray band defined in xbad
        xband = value at which L_2 was evaluated. default is 2 keV
        optband = value at which L_O was evaluated. default is 2500 A
        
    Outputs:
        alpha = alpha_ox value
    '''
    
    ### convert units ###
    xband = xband.to(u.um, equivalencies=u.spectral())
    
    optband = optband.to(u.um, equivalencies=u.spectral())
    
#    numer = np.log(optband.value * L_O.value) - np.log(xband.value * L_2.value)
#    denom = np.log(optband.value) - np.log(xband.value)
#    
#    alpha = -(numer/denom) + 1
    
    alpha = -0.3838 * (np.log10(L_2.value/L_O.value))
    
    mask = np.isinf(alpha)
    alpha[mask] = np.nan
    return alpha

def calc_alpha_kx(L_k, L_2, xband=2*u.keV, optband=1.6*u.um):
    ''' Function to calculate the alpha_kx value of any given objects
    Inputs:
        L_k = monochromatic luminosity in k band as defined in optband
        L_2 = monochromatic luminosity in xray band defined in xbad
        xband = value at which L_2 was evaluated. default is 2 keV
        optband = value at which L_O was evaluated. default is 1.6 um
        
    Outputs:
        alpha = alpha_kx value
    '''
    
    ### convert units ###
    xband = xband.to(u.um, equivalencies=u.spectral())
    
    optband = optband.to(u.um, equivalencies=u.spectral())
    
    numer = np.log(optband.value * L_k.value) - np.log(xband.value * L_2.value)
    denom = np.log(optband.value) - np.log(xband.value)
    
    alpha = -(numer/denom) + 1
    
#    alpha2 = -0.3838 * (np.log(L_2.value)/np.log(L_k.value))
    
    mask = np.isinf(alpha)
    alpha[mask] = np.nan
    return alpha

def get_L_k(tbdata):
    '''Function to get the monochromatic luminosity at 1.6um
    Inputs:
        tbdata = table of data for the objects
    Outputs:
        L_O = monochromatic luminosity at 1.6um in units of ergs/s/Hz
    '''
    ### get magnitude ###
    kmag = tbdata['M_K_z_p']
    
    ### convert to luminosity ###
    L_k = 10**((34.1-kmag)/2.5)
    L_k = L_k * u.W * (u.Hz)**-1 # Units of W/Hz
    
    L_k = L_k.to((u.erg) * (u.s)**-1 * (u.Hz)**-1, equivalencies=u.spectral())
    
    return L_k