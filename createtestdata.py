#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 17 11:49:11 2018

OLD CODE - NOT UPDATED WITH RESTRUCTURING 

code to create a mock catalogue that can be used to test my variability code.

@author: ppxee
"""

### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
from astropy.table import Table
#from astropy.table import Table #for handling tables
import numpy as np #for handling arrays
#import math
from astropy.stats import median_absolute_deviation
import vari_funcs #my module to help run code neatly
plt.close('all') #close any open plots

def add_variables(mag, num, sigma, indicies):
#        print(str(num)+' '+str(pervary)+'% variable sources added to data')
    for n in range(num):
        vals = sigma * np.random.randn(8) + 1 #gives normally distributed noise with sigma and mean=1
        idx = np.random.randint(0, size)
        while np.isin(idx, indicies):
            idx = np.random.randint(0, size)
    #        print('index already variable')
        mag[idx,:] = vals
        indicies = np.append(indicies, idx)
    return mag, indicies

noise_data = Table.read('median_noise_info.fits')
noise_sigmas = noise_data['median_noise_sigma']
noise_bins = noise_data['Low_mag_lim']

#go = True
count = 1
failcount = 0
passcount = 0
total =1# 10000

while count <= total: #go == True:
#    print(count)
    ### Create uniform flux data set ###
    size = 42000 # random number chosen 
    mag = np.ones([size,8]) 
    meanmag = np.empty(size)
    
    ### Check that this gives all 0 values of MAD ##
    mad = median_absolute_deviation(mag, axis=1)
    if np.sum(mad) != 0:
        print('FAIL: non zero MAD found in uniform data set')
        print(count)
        count += 1
        failcount += 1
        continue
    
    ### Add noise to the data ###        
    noisemag = np.empty([size,8])
    medsigma = np.empty(len(noise_sigmas))
    for l, signoise in enumerate(noise_sigmas):
        noisesize = 1000
        startind = l*noisesize
        for n in range(noisesize):
            noise = signoise * np.random.randn(8) #gives noise with specific sigma
            noisemag[startind+n,:] = mag[startind+n,:] + noise
            if l == 41:
                meanmag[startind+n] = np.random.randint(noise_bins[l],3e6)
            else:
                meanmag[startind+n] = np.random.randint(noise_bins[l],noise_bins[l+1])
        #get median std within bin to plot later
        sigmabin = np.std(noisemag[startind:startind+n,:], axis=1)
        medsigma[l] = np.nanmedian(sigmabin)
    indicies = np.array([]) # empty array of index values

    ### Add in variety variable sources ###
#    
#    pervary1 = 20 #percentage variability (50% = 0.5 peak to trough difference)
#    num1 = 100 # number of variable sources to be added
#    mag, indicies = add_variables(mag, num1, pervary1, indicies)
    
#    pervary2 = 40
#    num2 = 50
#    mag, indicies = add_variables(mag, num2, pervary2, indicies)
    
#    num = num1 #+ num2
            
    ### Find MAD and std values with noise ###
    noiseMAD = median_absolute_deviation(noisemag, axis=1)
    sigma = np.std(noisemag, axis=1)
    plt.figure(figsize=[9,6])
    plt.scatter(meanmag,sigma,marker='+')
    plt.xscale('log')
    plt.ylim(ymax=1.75)
    plt.xlim(xmin=1e2,)
    plt.ylabel('sigma')
    plt.xlabel('Mean Flux')
    plt.title('Test data with noise, no variable sources added')
    plt.plot(noise_bins, medsigma, 'k--')
    
    ### Use modified z-score ###
    modz = vari_funcs.mod_z_score(noiseMAD)
    thresh = 5 #set mod z threshold
    vary = modz[modz >= thresh]
    numnoise = np.size(vary)
#    print('Found '+str(numnoise)+' sources with mod-z >= '+str(thresh))
#    if numnoise > num:
#        print('FAIL: found too many variable source')
#        print(count)
#        count += 1
#        failcount += 1
#        continue
#    elif numnoise < num:
#        print('FAIL: found too few variables sources')
#        print(count)
#        count += 1
#        failcount += 1
#        continue
        
    count += 1
    passcount += 1
#    go = False #breaks the loop at end of tests
    
    
### Output results if run multiple times ###
successrate = (passcount/total)*100
print('Pass count = '+str(passcount))
print('Fail count = '+str(failcount))
print('% sucessful = ' + str(successrate))















