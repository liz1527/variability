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
import weightedstats as ws
import random
from scipy.stats import skew
import scipy.optimize
plt.close('all') #close any open plots
#from numpy.lib.recfunctions import append_fields

def getdata(data):
    j_flux = data['Month_Flux_J']
    j_fluxerr = data['Month_Fluxerr_J']
    k_flux = data['Month_Flux_K']
    k_fluxerr = data['Month_Fluxerr_K']
    
    return j_flux, j_fluxerr, k_flux, k_fluxerr

def parabola(x, a, b, c): #for curve fitting
    return a*x**2 + b*x + c

def weighted_mean_and_err(tau, ccf, cut_off=0.5):
    
    ### Limit to range of ccf function ##
    max_ccf = np.max(ccf)
    lim = cut_off * max_ccf
    tau = tau[ccf>lim]
    ccf = ccf[ccf>lim]
    
    ### find weighted mean and standard error on the mean ###
    tau_mean = ws.numpy_weighted_mean(tau, weights=ccf)
    tau_var = ws.numpy_weighted_mean((tau-tau_mean)**2, weights=ccf)
    tau_std = np.sqrt(tau_var)
    tau_SE = tau_std/len(tau)
    
    return tau_mean, tau_SE

#%% Import data for variables selected in J and K + define constants ###
varydata = Table.read('variable_tables/J_and_K_variables_month_varystats_DR11data.fits')
#varydata = varydata[varydata['X-ray']==True]
#varydata = Table.read('UDS_catalogues/chandra_month_varystats_noneg_DR11data.fits')

### Set constants for sample and bins ###
min_chi = 100 # chi cut required to create sample
dez = input('Should the light curves be de-redshifted? (Y or N) ')

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

#%% Create sample and x/nox subsamples ###
mask = varydata['Chi_K'] > min_chi
varydata = varydata[mask]
#varydata = varydata[varydata['X-ray']==True]
#varydata.sort(keys=key) #sort into increasing according to key
#varydata = varydata[varydata[key]!=0] # remove any with 0 as this is likely to be invalid
#varydata = varydata[varydata['z_use']<4] # remove any z>4 as this is likely to be invalid
#varydata = varydata[varydata[key]<4] # remove any >4 as this is likely to be invalid

#xvarydata = varydata[varydata['X-ray']==True]
#noxvarydata = varydata[varydata['X-ray']==False]


tau_arr = np.arange(-24,25) # tau values to evaluate ccf at 

#%% Get lightcurves with errors ###
j_flux, j_fluxerr, k_flux, k_fluxerr = getdata(varydata)

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

if dez == 'Y' or dez == 'y':
    ### do ccf with deredshifted light curves ####
    z = varydata['z_use']    
    out = np.array([vari_funcs.cross_correlation.cross_correlate_de_z(
            corr_test_k_flux, corr_test_j_flux, tau, z, type='dcf') for tau in tau_arr])
else:
    out = np.array([vari_funcs.cross_correlation.cross_correlate(
            corr_test_k_flux, corr_test_j_flux, tau, type='dcf') for tau in tau_arr])

### Unpack values ###
ccf = out[:,0]
ccf_err = out[:,1]

#%% Find weighted mean and skew of ccf ###
mean_lag, mean_lag_err = weighted_mean_and_err(tau_arr, ccf)
median_lag = ws.numpy_weighted_median(tau_arr, weights=ccf)
ccf_skew = skew(ccf)
max_lag = tau_arr[np.argmax(ccf)]

#%% Fit a parabola for those points around the centre of the ccf function ###
sub_tau = np.arange(-5,6)
test_ccf = ccf[np.isin(tau_arr, sub_tau)]
fit_params, pcov = scipy.optimize.curve_fit(parabola, sub_tau, test_ccf)
plot_tau = np.linspace(-5,6, 30)
ccf_fit = parabola(plot_tau, *fit_params)
max_lag_fit = plot_tau[np.argmax(ccf_fit)]
    
#%% Make plots ###
plt.figure(2,figsize=[10,10])
#plt.subplot(211)
#plt.plot(tau_arr, ccf,'o')
plt.errorbar(tau_arr, ccf, yerr=ccf_err, fmt='o')
#plt.vlines(mean_lag, -0.015,0.02, color='k', linestyle='dashed',
#           label='mean = '+str(mean_lag))
#plt.vlines(max_lag, -0.015,0.02, color='k', linestyle='dotted',
#           label='max = '+str(max_lag))
#plt.vlines(max_lag_fit, -0.015,0.02, color='k', linestyle='dotted',
#           label='fitted max = '+str(max_lag_fit))
#plt.plot(plot_tau, ccf_fit)
plt.xlabel('Lag (months)')
plt.ylabel('Cross-Correlation Function')
#    plt.ylim(-0.015,0.02)
plt.grid(True)
#plt.legend(loc='lower center')
plt.tight_layout()
    
end = time.time()
print(end-start)
    
    

    
    



