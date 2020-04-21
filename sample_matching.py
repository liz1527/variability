#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 12:18:07 2020

Code to create matched samples of X-ray and non-X-ray variables

@author: ppxee
"""


import time
start = time.time()
#print(start)

### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
from astropy.table import Table #for handling tables
import numpy as np #for handling arrays
#import math
#from astropy.stats import median_absolute_deviation
import vari_funcs #my module to help run code neatly
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
plt.close('all') #close any open plots
#from numpy.lib.recfunctions import append_fields
def get_differences(ob, noxsize, noxlum, noxz, abso=True):
    ### Extract arrays ###
    obsize = ob['FWHM_WORLD']*3600 # FWHM good estimate of size (converted to arcsecs)
    oblum = ob['M_K_z_p'] # Absolute mag = luminosity
    obz = ob['z_p'] # best photometric z so consistant with lum calcs
    
#    print('FWHM = '+str(obsize))
#    print('abso mag = '+str(oblum))
#    print('z = '+str(obz))

    ### Find differences ###
    sizediff = noxsize - obsize
    lumdiff = noxlum - oblum
    zdiff = noxz - obz
    
    if abso == True:
        ### Find absolute differences as we don't mind which is bigger ###
        abssizediff = np.abs(sizediff)
        abslumdiff = np.abs(lumdiff)
        abszdiff = np.abs(zdiff)
        return abssizediff, abslumdiff, abszdiff
    else:
        return sizediff, lumdiff, zdiff
    
def find_matches(ob, lumrange, zrange, noxsize, noxlum, noxz):
    print(ob['ID'])
#    ### Extract arrays ###
#    obsize = ob['FWHM_WORLD']*3600 # FWHM good estimate of size (converted to arcsecs)
#    oblum = ob['M_K_z_p'] # Absolute mag = luminosity
#    obz = ob['z_p'] # best photometric z so consistant with lum calcs
#    
##    print('FWHM = '+str(obsize))
##    print('abso mag = '+str(oblum))
##    print('z = '+str(obz))
#
#    ### Find differences ###
#    sizediff = noxsize - obsize
#    lumdiff = noxlum - oblum
#    zdiff = noxz - obz
#    
#    ### Find absolute differences as we don't mind which is bigger ###
#    abssizediff = np.abs(sizediff)
#    abslumdiff = np.abs(lumdiff)
#    abszdiff = np.abs(zdiff)
    ### Find absolute differences as we don't mind which is bigger ###
    abssizediff, abslumdiff, abszdiff = get_differences(ob, noxsize, noxlum, noxz)
    
    ### Sort by size match ###
    ''' argsort gives us a way of working down the list without changing order
    so we can work from closest size outwards and still identify the correct
    matched source in noxvarydata '''
    sortedinds = np.argsort(abssizediff)
    
    ### Iterate over sortedinds list to find objects that fit the criteria ###
    goodinds = np.array([])
    for ind in sortedinds:
        ### check luminosity is in range ###
        if abslumdiff[ind] < lumrange:
            lummatch = 1 
#            print('luminosity match')
        else:
            lummatch = 0
            continue # quit as 1 already out of range
        
        ### check z is in range ###
        if abszdiff[ind] < zrange:
            zmatch = 1
#            print('z match')
        else:
            zmatch = 0
        
        ### find if both in range ###
        match = lummatch * zmatch
        if match == 1:
#            print('Match found!')
            goodinds = np.append(goodinds, ind)
    goodinds = goodinds.astype(int)
    ### append total number of matches ###
    num = len(goodinds)
    if num == 0:
        print('No matches found')
        ### save inds as nans to indicate no matches ###
        closeind = np.nan
        bestsizeind = np.nan
    else:
        print(str(num)+' matches found')
    
        ### find closest overall match out of the good inds ###
        totdiff = abssizediff[goodinds] + abslumdiff[goodinds] + abszdiff[goodinds]
        closeind = goodinds[np.argmin(totdiff)]
        
        ### find closest in size ###
        bestsizeind = goodinds[np.argmin(abssizediff[goodinds])]
    return goodinds, closeind,  bestsizeind


#%% Open the fits files and split by X-ray detection ###
varydata = Table.read('variable_tables/J_and_K_variables_varystats_DR11data.fits')
xvarydata = varydata[varydata['X-ray']==True]
noxvarydata = varydata[varydata['X-ray']==False]

#%% get data arrays so don't extract on every loop ###
noxsize = noxvarydata['FWHM_WORLD']*3600 # FWHM good estimate of size (converted to arcsecs)
noxlum = noxvarydata['M_K_z_p'] # Absolute mag = luminosity
noxz = noxvarydata['z_p'] # best photometric z so consistant with lum calcs

#%% Define acceptable ranges ###
''' these are the max differences there can be from the match to the object '''
lumrange = 0.5
zrange = 0.2

#%% Iterate over objects ###
allcloseinds = np.array([])
allbestsizeinds = np.array([])
matches = np.array([])
for ob in xvarydata:
    goodinds, closeind, bestsizeind = find_matches(ob, lumrange, zrange, 
                                                   noxsize, noxlum, noxz)
    allcloseinds = np.append(allcloseinds, closeind)
    allbestsizeinds = np.append(allbestsizeinds, bestsizeind)
    
#    break # quit after one loop for when checking
print(str(len(allbestsizeinds[np.isnan(allbestsizeinds)])) + ' with no matches')

#%% Check each index to see if it is repeated ###
allbestinds2 = np.zeros(len(xvarydata)) + 9999 # so 9999 indicates those that no longer have matches
matchedinds = np.array([])
for n, test in enumerate(allbestsizeinds):
    if np.isnan(test) == True:
        allbestinds2[n] = test # so stays as no matches
        continue
    elif np.isin(test, matchedinds) == True: #i.e. the no x object already has a match
        continue # skip loop as already found best match for these objects
        
    test = test.astype(int) # as int so can use as index

    ### get values for the nox matched object ###
    testnoxsize = noxsize[test]
    testnoxlum = noxlum[test]
    testnoxz = noxz[test]
    
    ### get the matched x ray ob data ###
    mask = np.isin(allbestsizeinds, test)
    testobs = xvarydata[mask]
    
    ### create index array for ease of locating ###
    indarr = np.arange(len(xvarydata))
    
    ### Test which object is closest ###
    if len(testobs) == 1:
        allbestinds2[n] = test
        continue # as this means there is only one object matched to it
    else:
        testsizediff = np.array([])
        for ob in testobs:
            abssizediff, abslumdiff, abszdiff = get_differences(ob, testnoxsize, 
                                                                testnoxlum, testnoxz)
            testsizediff = np.append(testsizediff, abssizediff)
            
        ### get minimum size difference and isolate its data ###
        smalldiff = np.argmin(testsizediff)
        bestob = testobs[smalldiff]
        
        ### get original index in xvarydata ###
        mask = np.isin(xvarydata['ID'], bestob['ID'])
        
        ### set value in allbestinds2 to test at that index ###
        if allbestinds2[mask] == 9999:
            allbestinds2[mask] = test
        else:
            print('Error: object already matched')
            break # quit so I can build in capability for this
        
    ### save test in an array that indicates this has already been matched ###
    matchedinds = np.append(matchedinds, test)
        
#%% Rerun matching to get matches for those that were previously repeats ###
    
allbestinds3 = np.array([])
for n, ob in enumerate(xvarydata):
    if allbestinds2[n] != 9999: #ie it has a match or has no matches
        ### Set index and continue ###
        allbestinds3 = np.append(allbestinds3, allbestinds2[n])
        continue
    goodinds, closeind, bestsizeind = find_matches(ob, lumrange, zrange, 
                                                   noxsize, noxlum, noxz)
    while np.isin(bestsizeind, matchedinds) == True:
        ### This means the best match is taken so need to find next best ###
        goodinds = goodinds[goodinds!=bestsizeind] # delete best option
        testobs = noxvarydata[goodinds] # get nox data for other options
        if len(testobs) == 0:
            print('No further matches available')
            allbestinds3 = np.append(allbestinds3, np.nan) # set to nan as no match
            break # quit loop
        else:
            testnoxsize = noxsize[goodinds]
            testnoxlum = noxlum[goodinds]
            testnoxz = noxz[goodinds]
            goodinds2, closeind2, bestsizeind2 = find_matches(ob, lumrange, zrange, 
                                                           testnoxsize, 
                                                           testnoxlum, 
                                                           testnoxz)
            ### apply mask to goodinds to get original inds not inds of minitable ###
            goodinds = goodinds[goodinds2] 
            closeind = goodinds[closeind2]
            bestsizeind = goodinds[bestsizeind2]
    
    if len(allbestinds3) != n+1: #to test if nan has been returned
        ### set new ind value and append to matched inds ###
        allbestinds3 = np.append(allbestinds3, bestsizeind)
        matchedinds = np.append(matchedinds, bestsizeind)
    
print(str(len(allbestinds3[np.isnan(allbestinds3)])) + ' with no matches')
end = time.time()
print('time (s) = '+str(end-start))
    
    
    
    