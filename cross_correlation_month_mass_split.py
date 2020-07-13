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

#%% Open the fits files and get data ###
### Import data for variables selected in K ###
varydata = Table.read('variable_tables/J_and_K_variables_month_varystats_DR11data.fits')
#mask = np.nanmean(varydata['Month_Flux_K'], axis=1) > 1e4
mask = varydata['Chi_K'] > 100
varydata = varydata[mask]
xvarydata = varydata[varydata['X-ray']==True]
noxvarydata = varydata[varydata['X-ray']==False]

##%% Split data into 3 mass bins by value ###
#def make_mass_bin(varydata, low, upp):
#    mask1 = varydata['Mstar_z_p'] >= low
#    mask2 = varydata['Mstar_z_p'] < upp
#    mask = mask1*mask2.astype(bool)
#    bindata = varydata[mask]
#    return bindata
#lowbindata = make_mass_bin(varydata, low=5e8, upp=5e9)
#midbindata = make_mass_bin(varydata, low=5e9, upp=5e10)
#highbindata = make_mass_bin(varydata, low=5e10, upp=5e11)

#%% Split data into 3 mass bins by number in bin ###
varydata.sort(keys='Mstar_z_p') #sort into increasing mass
varydata = varydata[varydata['Mstar_z_p']!=0]
num = int(round(len(varydata)/3,0))
lowbindata = varydata[0:num]
midbindata = varydata[num:2*num]
highbindata = varydata[2*num:]

plt.figure()
bins = np.logspace(np.log10(varydata['Mstar_z_p'][0]),np.log10(varydata['Mstar_z_p'][-1]),25)
plt.hist([lowbindata['Mstar_z_p'],midbindata['Mstar_z_p'],highbindata['Mstar_z_p']], 
         bins, histtype='step')
plt.xscale('log')
#%% Get lightcurves with errors ###
def getdata(data):
    j_flux = data['Month_Flux_J']
    j_fluxerr = data['Month_Fluxerr_J']
    k_flux = data['Month_Flux_K']
    k_fluxerr = data['Month_Fluxerr_K']
    
    return j_flux, j_fluxerr, k_flux, k_fluxerr

low_j_flux, low_j_fluxerr, low_k_flux, low_k_fluxerr = getdata(lowbindata)
mid_j_flux, mid_j_fluxerr, mid_k_flux, mid_k_fluxerr = getdata(midbindata)
high_j_flux, high_j_fluxerr, high_k_flux, high_k_fluxerr = getdata(highbindata)

#%% Mean subract and normalise ###

test_low_j_flux = vari_funcs.cross_correlation.weighted_mean_subtract_normalise(
        low_j_flux, low_j_fluxerr)
test_low_k_flux = vari_funcs.cross_correlation.weighted_mean_subtract_normalise(
        low_k_flux, low_k_fluxerr)

test_mid_j_flux = vari_funcs.cross_correlation.weighted_mean_subtract_normalise(
        mid_j_flux, mid_j_fluxerr)
test_mid_k_flux = vari_funcs.cross_correlation.weighted_mean_subtract_normalise(
        mid_k_flux, mid_k_fluxerr)

test_high_j_flux = vari_funcs.cross_correlation.weighted_mean_subtract_normalise(
        high_j_flux, high_j_fluxerr)
test_high_k_flux = vari_funcs.cross_correlation.weighted_mean_subtract_normalise(
        high_k_flux, high_k_fluxerr)
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
corr_test_low_j_flux = vari_funcs.cross_correlation.make_corr_arrays(test_low_j_flux, 
                                                                 jmask, 
                                                                 full_months)
corr_test_low_k_flux = vari_funcs.cross_correlation.make_corr_arrays(test_low_k_flux, 
                                                                 kmask,
                                                                 full_months)

corr_test_mid_j_flux = vari_funcs.cross_correlation.make_corr_arrays(test_mid_j_flux, 
                                                                 jmask, 
                                                                 full_months)
corr_test_mid_k_flux = vari_funcs.cross_correlation.make_corr_arrays(test_mid_k_flux, 
                                                                 kmask,
                                                                 full_months)

corr_test_high_j_flux = vari_funcs.cross_correlation.make_corr_arrays(test_high_j_flux, 
                                                                 jmask, 
                                                                 full_months)
corr_test_high_k_flux = vari_funcs.cross_correlation.make_corr_arrays(test_high_k_flux, 
                                                                 kmask,
                                                                 full_months)
#### Make test array for J out of K ###
#nanarr = np.empty([len(corr_test_k_flux),4])#, dtype=int)
#nanarr[:,:] = np.nan
#corr_test_j_flux = np.copy(corr_test_k_flux)
#corr_test_j_flux = np.delete(corr_test_j_flux, [0,1,2,3], axis=1)
#corr_test_j_flux = np.append(corr_test_j_flux, nanarr, axis=1)
#
#plt.figure()
#plt.plot(full_months, corr_test_k_flux[0,:],'o')
#plt.plot(full_months, corr_test_j_flux[0,:],'o')


#%% Calculate the CCF at various tau values ###
    
tau_arr = np.arange(-50,50)
#tau_arr = np.arange(-90,90)

low_out = np.array([vari_funcs.cross_correlation.cross_correlate(
        corr_test_low_k_flux, corr_test_low_j_flux, tau) for tau in tau_arr])

mid_out = np.array([vari_funcs.cross_correlation.cross_correlate(
        corr_test_mid_k_flux, corr_test_mid_j_flux, tau) for tau in tau_arr])

high_out = np.array([vari_funcs.cross_correlation.cross_correlate(
        corr_test_high_k_flux, corr_test_high_j_flux, tau) for tau in tau_arr])

### Unpack values ###
low_ccf = low_out[:,0]
low_ccf_err = low_out[:,1]

mid_ccf = mid_out[:,0]
mid_ccf_err = mid_out[:,1]

high_ccf = high_out[:,0]
high_ccf_err = high_out[:,1]

### get bin edges ###
low_min = round(np.min(lowbindata['Mstar_z_p']),0)
low_max = round(np.max(lowbindata['Mstar_z_p']),0)
mid_min = round(np.min(midbindata['Mstar_z_p']),0)
mid_max = round(np.max(midbindata['Mstar_z_p']),0)
high_min = round(np.min(highbindata['Mstar_z_p']),0)
high_max = round(np.max(highbindata['Mstar_z_p']),0)

#%% Make plots ###
plt.figure(figsize=[10,10])
#plt.subplot(211)
#plt.plot(tau_arr, ccf,'o')
plt.errorbar(tau_arr, low_ccf, yerr=low_ccf_err, fmt='o',
             label='{:.2e}'.format(low_min)+r' $ M\odot \leq M_{star} < $' +
             '{:.2e}'.format(low_max)+r'$ M\odot$')
plt.errorbar(tau_arr, mid_ccf, yerr=mid_ccf_err, fmt='o',
             label='{:.2e}'.format(mid_min)+r' $ M\odot \leq M_{star} < $' +
             '{:.2e}'.format(mid_max)+r'$ M\odot$')
plt.errorbar(tau_arr, high_ccf, yerr=high_ccf_err, fmt='o',
             label='{:.2e}'.format(high_min)+r' $ M\odot \leq M_{star} < $' +
             '{:.2e}'.format(high_max)+r'$ M\odot$')

plt.xlabel('Lag (months)')
plt.ylabel('Cross-Correlation Function')
plt.ylim(-0.015,0.025)
plt.grid(True)
plt.legend()
plt.tight_layout()


plt.figure(figsize=[10,10])
#plt.subplot(211)
#plt.plot(tau_arr, ccf,'o')
plt.errorbar(tau_arr, low_ccf, yerr=low_ccf_err, fmt='o', color='C0',
             label='{:.2e}'.format(low_min)+r' $ M\odot \leq M_{star} < $' +
             '{:.2e}'.format(low_max)+r'$ M\odot$')
plt.xlabel('Lag (months)')
plt.ylabel('Cross-Correlation Function')
plt.ylim(-0.015,0.025)
plt.grid(True)
plt.legend()
plt.tight_layout()


plt.figure(figsize=[10,10])
#plt.subplot(211)
#plt.plot(tau_arr, ccf,'o')
plt.errorbar(tau_arr, mid_ccf, yerr=mid_ccf_err, fmt='o', color='C1',
             label='{:.2e}'.format(mid_min)+r' $ M\odot \leq M_{star} < $' +
             '{:.2e}'.format(mid_max)+r'$ M\odot$')

plt.xlabel('Lag (months)')
plt.ylabel('Cross-Correlation Function')
plt.ylim(-0.015,0.025)
plt.grid(True)
plt.legend()
plt.tight_layout()


plt.figure(figsize=[10,10])
#plt.subplot(211)
#plt.plot(tau_arr, ccf,'o')
plt.errorbar(tau_arr, high_ccf, yerr=high_ccf_err, fmt='o', color='C2',
             label='{:.2e}'.format(high_min)+r' $ M\odot \leq M_{star} < $' +
             '{:.2e}'.format(high_max)+r'$ M\odot$')

plt.xlabel('Lag (months)')
plt.ylabel('Cross-Correlation Function')
plt.ylim(-0.015,0.025)
plt.grid(True)
plt.legend()
plt.tight_layout()
            
end = time.time()
print(end-start)
    
    
    
    
    



