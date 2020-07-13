#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 13:58:17 2020

Code to run a stacked cross correlation on J and K month lightcurves

@author: ppxee
"""

import time
start = time.time()
print(start)

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
import random
plt.close('all') #close any open plots
#from numpy.lib.recfunctions import append_fields

def getdata(data):
    j_flux = data['Month_Flux_J']
    j_fluxerr = data['Month_Fluxerr_J']
    k_flux = data['Month_Flux_K']
    k_fluxerr = data['Month_Fluxerr_K']
    
    return j_flux, j_fluxerr, k_flux, k_fluxerr

#%% Import data for variables selected in J and K + define constants ###
xvarydata = Table.read('variable_tables/J_and_K_variables_month_varystats_DR11data_chan_xray.fits')
noxvarydata = Table.read('variable_tables/J_and_K_variables_month_varystats_DR11data_chan_notxray.fits')
#varydata = varydata[varydata['X-ray']==True]

### Set constants for sample and bins ###
min_chi = 100 # chi cut required to create sample
key = '$L_{X}$' # column to split on
unit = '$erg s^{-1} cm^{-2}$' #for label
num_bins = 2 # how many bins to create
tau_arr = np.arange(-50,50) # tau values to evaluate ccf at 
#tau_arr = np.arange(-90,90)
#log = True
log = False

### Get monthly_numbers as this contains all possible months ###
month_info = fits.open('Images/Convolving_Images/monthly_numbers.fits')[1].data #get month count data
full_months = month_info['Month'] #extract month nanes

### Create J and K month arrays so can match where the data should go ###
kmonths = ['sep05','oct05','nov05','dec05', 'jan06', #'dec06', 
          'jan07', 'aug07', 'sep07', 'oct07', 'sep08', 'oct08', 'nov08', 
          'jul09', 'aug09', 'sep09', 'oct09', 'nov09', 'dec09', 'jan10', 
          'feb10', 'aug10', 'sep10', 'oct10', 'nov10', 'dec10', 'jan11', #'feb11', 
          'aug11', 'sep11', 'oct11', 'nov11', 'dec11', 'jan12', 'feb12', 
          'jul12', 'aug12', 'sep12', 'oct12', 'nov12']
kmask = np.isin(full_months, kmonths)

jmonths = ['sep05', 'oct05', 'nov05', 'dec05', 'jan06', 'oct06', 'nov06',
          'dec06', 'aug07', 'sep07', 'oct07', 'oct08', 'nov08', 'aug09',
          'sep09', 'oct09', 'nov09', 'dec09', 'aug10', 'sep10', 'oct10',
          'nov10', 'dec10', 'jan11', 'aug11', 'sep11', 'oct11', 'nov11',
          'dec11', 'jan12', 'jul12', 'aug12', 'sep12', 'oct12', 'nov12']
jmask = np.isin(full_months, jmonths)

#%% Create samples###
xmask = xvarydata['Chi_K'] > min_chi
xvarydata = xvarydata[xmask]
noxmask = noxvarydata['Chi_K'] > min_chi
noxvarydata = noxvarydata[noxmask]

#%% Split xray data data into bins by number in bin ###
### sort into increasing according to xray lum ###
_, x_xlum, _ = vari_funcs.xray_funcs.get_xray_L(xvarydata)
xvarydata = xvarydata[x_xlum!=0] # remove any with 0 as this is likely to be invalid
x_xlum = x_xlum[x_xlum!=0] # remove any with 0 as this is likely to be invalid
inds = np.argsort(x_xlum)
x_xlum = x_xlum[inds].value
xvarydata = xvarydata[inds] 

num = int(round(len(xvarydata)/num_bins,0)) # calculate number in each bin
if log == True:
    bins = np.logspace(np.log10(x_xlum[0]),
                       np.log10(x_xlum[-1]),30)
else:    
    bins = np.linspace(x_xlum[0], x_xlum[-1],30)

for n in range(num_bins+1):
    if n == 2:
        bindata = noxvarydata # make nox its own bin
        label = 'Not X-ray Detected'
    else:
        bindata = xvarydata[n*num:(n+1)*num] # define bin according to number in bin
        ### get bin edges ###
        bin_min = round(np.min(x_xlum[n*num:(n+1)*num]),2)
        bin_max = round(np.max(x_xlum[n*num:(n+1)*num]),2)
        ### set label ###
        label = '{:.2e}'.format(bin_min)+unit+r' $\leq$ '+key+r' $<$ ' + '{:.2e}'.format(bin_max)+unit
    print(len(bindata))
    
    plt.figure(1)
    plt.hist(x_xlum[n*num:(n+1)*num], bins, histtype='step')
    if log == True:
        plt.xscale('log')   
    plt.xlabel(key)
    plt.tight_layout()
    
    #%% Get lightcurves with errors ###
    j_flux, j_fluxerr, k_flux, k_fluxerr = getdata(bindata)

    #%% Mean subract and normalise ###
    
    test_j_flux = vari_funcs.cross_correlation.weighted_mean_subtract_normalise(
            j_flux, j_fluxerr)
    test_k_flux = vari_funcs.cross_correlation.weighted_mean_subtract_normalise(
            k_flux, k_fluxerr)
    
    #%% Create correlation arrays ###
    ''' Need arrays that have a space for every possible month so that the values 
    can be separated by the correct time periods'''

    ### Assign new arrays ###
    corr_test_j_flux = vari_funcs.cross_correlation.make_corr_arrays(test_j_flux, 
                                                                     jmask, 
                                                                     full_months)
    corr_test_k_flux = vari_funcs.cross_correlation.make_corr_arrays(test_k_flux, 
                                                                     kmask,
                                                                     full_months)

    #%% Calculate the CCF at various tau values ###
        
    out = np.array([vari_funcs.cross_correlation.cross_correlate(
            corr_test_k_flux, corr_test_j_flux, tau) for tau in tau_arr])

    ### Unpack values ###
    ccf = out[:,0]
    ccf_err = out[:,1]


    #%% Make plots ###
    plt.figure(2,figsize=[10,10])
    #plt.subplot(211)
    #plt.plot(tau_arr, ccf,'o')
    plt.errorbar(tau_arr, ccf, yerr=ccf_err, fmt='o',
                 label=label)
    
    plt.xlabel('Lag (months)')
    plt.ylabel('Cross-Correlation Function')
    plt.ylim(-0.015,0.025)
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    
    
    plt.figure(n+3,figsize=[10,10])
    #plt.subplot(211)
    #plt.plot(tau_arr, ccf,'o')
    plt.errorbar(tau_arr, ccf, yerr=ccf_err, fmt='o', color='C'+str(n),
                 label=label)
    plt.xlabel('Lag (months)')
    plt.ylabel('Cross-Correlation Function')
    plt.ylim(-0.015,0.025)
    plt.grid(True)
    plt.legend()
    plt.tight_layout()


end = time.time()
print(end-start)
    
    
    
    
    



