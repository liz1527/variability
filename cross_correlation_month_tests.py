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
plt.close('all') #close any open plots
#from numpy.lib.recfunctions import append_fields

#%% Open the fits files and get data ###
### Import data for variables selected in K ###
varydata = Table.read('variable_tables/J_and_K_variables_month_varystats_DR11data.fits')
#mask = np.nanmean(varydata['Month_Flux_K'], axis=1) > 1e4
mask = varydata['Chi_K'] > 100
varydata = varydata[mask]
xvarydata = varydata[varydata['X-ray']==True]
noxvarydata = varydata[varydata['X-ray']==False]


#%% Get lightcurves with errors ###
def getdata(data):
    j_flux = data['Month_Flux_J']
    j_fluxerr = data['Month_Fluxerr_J']
    k_flux = data['Month_Flux_K']
    k_fluxerr = data['Month_Fluxerr_K']
    
    return j_flux, j_fluxerr, k_flux, k_fluxerr

j_flux, j_fluxerr, k_flux, k_fluxerr = getdata(varydata)
x_j_flux, x_j_fluxerr, x_k_flux, x_k_fluxerr = getdata(xvarydata)
nox_j_flux, nox_j_fluxerr, nox_k_flux, nox_k_fluxerr = getdata(noxvarydata)

mask = np.nanmean(k_flux, axis=1) > 1e4
#mask = varydata['z'] > 1.5
high_j_flux = j_flux[mask]
high_k_flux = k_flux[mask]
low_j_flux = j_flux[~mask]
low_k_flux = k_flux[~mask]
high_j_fluxerr = j_fluxerr[mask]
high_k_fluxerr = k_fluxerr[mask]
low_j_fluxerr = j_fluxerr[~mask]
low_k_fluxerr = k_fluxerr[~mask]

#%% Mean subract and normalise ###
#
#test_j_flux = vari_funcs.cross_correlation.mean_subtract_normalise(x_j_flux)
#test_k_flux = vari_funcs.cross_correlation.mean_subtract_normalise(x_k_flux)

test_j_flux = vari_funcs.cross_correlation.weighted_mean_subtract_normalise(x_j_flux, 
                                                                            x_j_fluxerr)
test_k_flux = vari_funcs.cross_correlation.weighted_mean_subtract_normalise(x_k_flux, 
                                                                            x_k_fluxerr)
#%% Create correlation arrays ###
''' Need arrays that have a space for every possible month so that the values 
can be separated by the correct time periods'''

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

### Assign new arrays ###
corr_test_j_flux = vari_funcs.cross_correlation.make_corr_arrays(test_j_flux, 
                                                                 jmask, 
                                                                 full_months)
corr_test_k_flux = vari_funcs.cross_correlation.make_corr_arrays(test_k_flux, 
                                                                 kmask, 
                                                                 full_months)

    
tau_arr = np.arange(-50,50)
#tau_arr = np.arange(-90,90)

out = np.array([vari_funcs.cross_correlation.cross_correlate(corr_test_k_flux, 
                                                             corr_test_j_flux, 
                                                             tau) for tau in tau_arr])

ccf = out[:,0]
ccf_err = out[:,1]
num_pairs = out[:,2]

plt.figure()
#plt.subplot(211)
#plt.plot(tau_arr, ccf,'o')
plt.errorbar(tau_arr, ccf, yerr=ccf_err, fmt='o')
plt.xlabel('Lag (months)')
plt.ylabel('Cross-Correlation Function')
plt.grid(True)
plt.tight_layout()

    
#plt.subplot(212)
plt.figure()
plt.plot(tau_arr, num_pairs,'o')
plt.xlabel('Lag (months)')
plt.ylabel('Pairs per Light Curve')
plt.tight_layout()

#%% Test how average lag changes with CCF Value ###
#round_ccf = np.round(ccf, decimals=4)
#test_values = np.arange(-0.008, 0.009,0.0001)
#test_values = np.round(test_values, decimals=4)
round_ccf = np.round(ccf, decimals=3)
test_values = np.arange(-0.008, 0.009,0.001)
test_values = np.round(test_values, decimals=3)
mean_tau = np.array([])
for val in test_values:
    ### Find places where ccf rounds to the test value ###
    mask = np.isin(round_ccf, val)
    
    ### Find Tau at these and average to find centre of the peak at that ccf value ###
    taus = tau_arr[mask]
    print(len(taus))
    if len(taus) <= 1:
        mean_tau = np.append(mean_tau, np.nan)
    else:
        mean_tau = np.append(mean_tau, np.nanmean(taus))

plt.figure()
plt.plot(tau_arr, round_ccf,'o')
plt.xlabel('Mean lag (months)')
plt.ylabel('Rounded CCF Value')
plt.grid(True)

#plt.figure()
plt.plot(mean_tau, test_values,'o')
#plt.xlabel('Mean lag (months)')
#plt.ylabel('CCF Value')
##plt.grid(True)
    
#%% Test shifting some of the lightcurves artifically ###
#mask = [10,20,30,40,50,60,70]
import random
def make_shifted_data(corr_test_k_flux, lag, percentage):
    ### Use percentage to calculate the number of lightcurves to shift ###
    length = len(corr_test_k_flux)
    no_steps = int(length*percentage/100) # number of lcs to shift based on %
    mask = random.sample(range(0, length), no_steps)
    shift = np.arange(0,lag,1, dtype=int)
    ### Make test array for J out of K ###
    nanarr = np.empty([len(mask),lag])#, dtype=int)
    nanarr[:,:] = np.nan
    corr_test_j_flux_new = np.copy(corr_test_k_flux)
    corr_test_j_flux_moved = np.delete(corr_test_j_flux_new[mask,:], shift, axis=1)
    corr_test_j_flux_moved = np.append(corr_test_j_flux_moved, nanarr, axis=1)
    
    corr_test_j_flux_new[mask,:] = corr_test_j_flux_moved
    return corr_test_j_flux_new

plt.figure(figsize=[10,10])
#plt.plot(full_months, corr_test_k_flux[0,:],'o')
#plt.plot(full_months, corr_test_j_flux[0,:],'o')
#for lag in [10]:#[0,5,10,15,20,25,30,35,40,45]:
#    corr_test_j_flux_new = make_shifted_data(corr_test_k_flux, lag=lag, percentage=100)
#    
#    test_out = np.array([
#            vari_funcs.cross_correlation.cross_correlate_shifted(corr_test_k_flux,
#                                                                 corr_test_j_flux_new, 
#                                                                 tau) for tau in tau_arr])
#    
#    test_ccf = test_out[:,0]
#    test_ccf_err = test_out[:,1]
#    test_num_pairs = test_out[:,2]
#    
#    plt.plot(tau_arr, test_ccf,'-', label='20% with lag = '+str(lag)+' months')
#
##    plt.errorbar(tau_arr, test_ccf,yerr=ccf_err, fmt='-', label='20% with lag = '+str(lag)+' months')
#
##plt.figure()
##plt.subplot(211)
#plt.plot(tau_arr, ccf,'o', label='Real data')
##plt.figure()
##plt.subplot(211)
##plt.plot(tau_arr, test_ccf,'o')
#plt.xlabel('Lag (months)')
#plt.ylabel('Cross-Correlation Function')
#plt.grid(True)
#plt.legend()
#plt.tight_layout()

plt.figure(figsize=[10,10])
plt.errorbar(tau_arr, ccf, yerr=ccf_err, fmt='o', label='Real data')
for percent in [5]:#[5,10,15,20,25,30,35,40,45]:
    corr_test_j_flux_new = make_shifted_data(corr_test_k_flux, lag=10, percentage=percent)

    test_out = np.array([
            vari_funcs.cross_correlation.cross_correlate_shifted(corr_test_k_flux,
                                                                 corr_test_j_flux_new, 
                                                                 tau) for tau in tau_arr])    
    test_ccf = test_out[:,0]
    test_num_pairs = test_out[:,1]
    
    plt.plot(tau_arr, test_ccf, 'o', label=str(percent)+'% with lag = 10 months')
#    plt.errorbar(tau_arr, test_ccf,yerr=ccf_err, fmt='-', 
#                 label=str(percent)+'% with lag = 10 months')

#plt.figure()
#plt.subplot(211)
#plt.figure()
#plt.subplot(211)
#plt.plot(tau_arr, test_ccf,'o')
plt.xlabel('Lag (months)')
plt.ylabel('Cross-Correlation Function')
plt.grid(True)
plt.legend()
plt.tight_layout()
##plt.subplot(212)
#plt.figure()
#plt.plot(tau_arr, test_num_pairs,'o')
#plt.xlabel('Lag (months)')
#plt.ylabel('Pairs per Light Curve')
#plt.tight_layout()














