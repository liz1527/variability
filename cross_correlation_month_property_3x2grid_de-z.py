#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 13:58:17 2020

Code to run a stacked cross correlation on J and K month lightcurves in mass 
and redshift bins.

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
from scipy.interpolate import splev, splrep
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

def double_parabola(x, a, b, c, d, e, f): #for curve fitting
    return (a*x**2 + b*x + c) * (d*x**2 + e*x + f)
#%% Import data for variables selected in J and K + define constants ###
varydata = Table.read('variable_tables/J_and_K_variables_month_varystats_DR11data.fits')
#varydata = varydata[varydata['X-ray']==True]

### Set constants for sample and bins ###
min_chi = 100 # chi cut required to create sample
    
#save = input('Save plots? (True or False) ')
#if save == 'True':
#    save = True
#else:
#    save = False
#
#if save == True:
#    subfolder = input('subfolder? (leave blank if none, include / at end) ')
#    filepath = 'plots/new_catalogue/JK_lightcurve_comp/dcf/grid/'+subfolder+''

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

testmonths = ['nov05', 'nov06', 'sep07', 'oct08', 'oct09', 'oct10', 'nov11', 'sep12']
testmask = np.isin(full_months, testmonths)

### set up month tick details ###
month_info = fits.open('Images/Convolving_Images/monthly_numbers.fits')[1].data #get month count data
full_months = month_info['Month'] #extract month nanes
tick_inds = np.load('Images/Convolving_Images/tick_inds_K.npy') #load tick locations
inds = np.arange(len(full_months)) #load tick locations
mask = np.zeros(len(full_months)) #set up mask
mask[tick_inds] = 1
mask = mask.astype(bool)
month_ticks = np.copy(full_months)
#month_ticks = month_ticks[mask]#retrieve tick details
month_ticks[~mask] = ''#set labels to blank

x = np.arange(0, len(month_info['Frames in v11']))
kmask = np.isin(full_months, kmonths)
Kx_months = x[kmask]
jmask = np.isin(full_months, jmonths)
Jx_months = x[jmask]
testmask = np.isin(full_months, testmonths)
testx_months = x[testmask]

#%% Create sample and split into grid bins ###
### get sample that has hi chi values and no anomalous masses or z
varydata = varydata[varydata['Chi_K'] > min_chi]
varydata = varydata[varydata['Mstar_z_p']!=0] # remove any with 0 as this is likely to be invalid
varydata = varydata[varydata['z_p']<5] # remove any z>5 as this is likely to be invalid


### first sort and split on mass ###
key = 'Mstar_z_p'
unit = '$M_{\odot}$' #for label
zkey = 'z_use'
zunit = '' #for label
varydata.sort(keys=key) #sort into increasing according to key#varydata = varydata[varydata[key]<4] # remove any >4 as this is likely to be invalid

num_m_bins = 3
num = int(round(len(varydata)/num_m_bins,0)) # calculate number in each bin

num_bins=9
def make_bins(varydata, mbins, zbins):
    
    ### top left ###
    mask1 = varydata[key] >= mbins[1] # high mass
    mask2 = varydata[zkey] < zbins[0] # low z
    mask = mask1*mask2.astype(bool)
    bindata1 = varydata[mask]
    plt.scatter(bindata1[zkey], bindata1[key], label='top left')
    print(len(bindata1))
    
    ### top middle ###
    mask1 = varydata[key] >= mbins[1] # high mass
    mask2 = varydata[zkey] >= zbins[0] # mid z
    mask3 = varydata[zkey] < zbins[1] # mid z
    mask = mask1*mask2*mask3.astype(bool)
    bindata2 = varydata[mask]
    plt.scatter(bindata2[zkey], bindata2[key], label='top middle')
    print(len(bindata2))
    
    ### middle left ###
    mask1 = varydata[key] < mbins[1] # mid mass
    mask2 = varydata[key] >= mbins[0] # mid mass
    mask3 = varydata[zkey] < zbins[0] # low z
    mask = mask1*mask2*mask3.astype(bool)
    bindata4 = varydata[mask]
    plt.scatter(bindata4[zkey], bindata4[key], label='middle left')
    print(len(bindata4))
    
    ### middle middle ###
    mask1 = varydata[key] < mbins[1] # mid mass
    mask2 = varydata[key] >= mbins[0] # mid mass
    mask3 = varydata[zkey] >= zbins[0] # mid z
    mask4 = varydata[zkey] < zbins[1] # mid z
    mask = mask1*mask2*mask3*mask4.astype(bool)
    bindata5 = varydata[mask]
    plt.scatter(bindata5[zkey], bindata5[key], label='middle')
    print(len(bindata5))

    
    ### bottom left ###
    mask1 = varydata[key] < mbins[0] # low mass
    mask2 = varydata[zkey] < zbins[0] # low z
    mask = mask1*mask2.astype(bool)
    bindata7 = varydata[mask]
    plt.scatter(bindata7[zkey], bindata7[key], label='bottom left')
    print(len(bindata7))
    
    ### bottom middle ###
    mask1 = varydata[key] < mbins[0] # low mass
    mask2 = varydata[zkey] >= zbins[0] # mid z
    mask3 = varydata[zkey] < zbins[1] # mid z
    mask = mask1*mask2*mask3.astype(bool)
    bindata8 = varydata[mask]
    plt.scatter(bindata8[zkey], bindata8[key], label='bottom middle')
    print(len(bindata8))
    
    plt.hlines(mbins, 0, 2.5)
    plt.vlines(zbins, 2e9, 5e11)
#    plt.legend()
    return bindata1, bindata2, bindata4, bindata5, bindata7, bindata8
    
### Set up Mstar z plot so see bin ranges ###
plt.figure(1)
#plt.scatter(varydata[zkey], varydata[key])
plt.yscale('log')
plt.ylabel('Stellar Mass')
plt.xlabel('Redshift')

### Define bin boundaries ###
mbins = np.array([2.3e10, 5.65e10])
zbins = np.array([1.11, 2])

### get bins ### 
bindatas = make_bins(varydata, mbins, zbins)

#%% Run ccf on bins ###

### Set up arrays for mean vs mass plot ###
max_lag = np.zeros(num_bins)
max_lag_fit = np.zeros(num_bins)
spline_max_lag = np.zeros(num_bins)
spline_max_lag_fit = np.zeros(num_bins)
tau_arr = np.arange(-24,25) # tau values to evaluate ccf at 

### set up coord arrays for text (set up this way to make it clear where they are for) ###
xcoords = [0.4, 1.3,
           0.4, 1.3,
           0.4, 1.3]
ycoords = [1.5e11, 1.5e11,
           3.25e10, 3.25e10,
           1e10, 1e10]

### Set up Mstar z plot so see ccf values ###
plt.figure(2)
plt.hlines(mbins, 0, 5)
plt.vlines(zbins, 2e9, 5e11)
plt.yscale('log')
plt.ylabel('Stellar Mass')
plt.xlabel('Redshift')
plt.title('Real data, raw max')

plt.figure(3)
plt.hlines(mbins, 0, 5)
plt.vlines(zbins, 2e9, 5e11)
plt.yscale('log')
plt.ylabel('Stellar Mass')
plt.xlabel('Redshift')
plt.title('Real data, fitted max')

plt.figure(4)
plt.hlines(mbins, 0, 5)
plt.vlines(zbins, 2e9, 5e11)
plt.yscale('log')
plt.ylabel('Stellar Mass')
plt.xlabel('Redshift')
plt.title('Spline data, raw max')

plt.figure(5)
plt.hlines(mbins, 0, 5)
plt.vlines(zbins, 2e9, 5e11)
plt.yscale('log')
plt.ylabel('Stellar Mass')
plt.xlabel('Redshift')
plt.title('Spline data, fitted max')

plt.figure(6, figsize=[10,10])
plt.suptitle('Raw data')

plt.figure(7, figsize=[10,10])
plt.suptitle('Splines')

plt.figure(8, figsize=[10,10])
plt.suptitle('Pair Numbers')

plt.figure(9, figsize=[10,10])
plt.suptitle('Pair Numbers Splines')

labels = np.array(['top left', 'top right',
                   'middle left', 'middle right',
                   'bottom left', 'bottom right'])
for n, bindata in enumerate(bindatas):
    print(n)
    print('running ccf using raw data')
    #%% Get lightcurves with errors ###
    j_flux, j_fluxerr, k_flux, k_fluxerr = getdata(bindata)

    #%% Mean subract and normalise ###
    
    test_j_flux, test_j_fluxerr = vari_funcs.cross_correlation.weighted_mean_subtract_normalise_errs(
            j_flux, j_fluxerr)
    test_k_flux, test_k_fluxerr = vari_funcs.cross_correlation.weighted_mean_subtract_normalise_errs(
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
    #%% De-redshift light curves ###
    ''' Need to shift the arrays into the restframe using individual redshifts'''
    
    #%% Calculate the CCF at various tau values on de-redshifted liightcurves ###
    z = bindata['z_use']    
    out = np.array([vari_funcs.cross_correlation.cross_correlate_de_z(
            corr_test_k_flux, corr_test_j_flux, tau, z, type='dcf') for tau in tau_arr])

    ### Unpack values ###
    ccf = out[:,0]
    ccf_err = out[:,1]
    pair_num = out[:,2]

    #%% Find weighted mean and skew of ccf ###
    max_lag[n] = tau_arr[np.argmax(ccf)]
    
    #%% Fit a parabola for those points around the centre of the ccf function ###
    sub_tau = np.arange(-7,8)
    test_ccf = ccf[np.isin(tau_arr, sub_tau)]
    fit_params, pcov = scipy.optimize.curve_fit(parabola, sub_tau, test_ccf)
    plot_tau = np.linspace(-7,7, 100)
    ccf_fit = parabola(plot_tau, *fit_params)
    max_lag_fit[n] = plot_tau[np.argmax(ccf_fit)]
    
    #%% Make grid plots ###
    plt.figure(2)
    plt.text(xcoords[n], ycoords[n], str(max_lag[n]))
    plt.figure(3)
    plt.text(xcoords[n], ycoords[n], str(round(max_lag_fit[n],2)))
    
    #%% Plot CCF and fit ###
    plt.figure(6)
    plt.subplot(3,2,n+1)
    plt.errorbar(tau_arr, ccf, yerr=ccf_err, fmt='o', zorder=0)
    plt.plot(plot_tau, ccf_fit, zorder=1)
    plt.vlines(max_lag[n], np.min(ccf), np.max(ccf),linestyle='dashed', 
               label='max = '+str(max_lag[n]))
    plt.vlines(max_lag_fit[n], np.min(ccf), np.max(ccf),linestyle='dashdot', 
               label='fitted max = '+str(round(max_lag_fit[n],2)))
    plt.legend()
    plt.xlabel('Lag (months)')
    plt.ylabel('CCF')
    plt.title(labels[n])
    
    #%% Plot Number of Pairs ###
    plt.figure(8)
    plt.subplot(3,2,n+1)
    plt.plot(tau_arr, pair_num, 'o')
    plt.legend()
    plt.xlabel('Lag (months)')
    plt.ylabel('Number of Pairs')
    plt.title(labels[n])
    
    #%% Repeat with splines ###
    print('running ccf using splines')
    j_splines = []
    k_splines = []
    test_x = np.linspace(0,len(full_months)-4,100)
    knots = testx_months[1:-1]
    new_test_j_flux = np.zeros(np.shape(test_j_flux))
    new_test_k_flux = np.zeros(np.shape(test_k_flux))
    for m in range(len(bindata)):
        j_spl = splrep(Jx_months, test_j_flux[m,:], w=1/test_j_fluxerr[m,:], t=knots)
        j_splines.append(j_spl)
        k_spl = splrep(Kx_months, test_k_flux[m,:], w=1/test_k_fluxerr[m,:], t=knots)
        k_splines.append(k_spl)
        
        ###  use splines to make test arrays ###
        new_test_j_flux[m,:] = splev(Jx_months, j_spl)
        new_test_k_flux[m,:] = splev(Kx_months, k_spl)
    
    spline_test_j_flux = new_test_j_flux
    spline_test_k_flux = new_test_k_flux
    
    #%% Create correlation arrays ###
    ''' Need arrays that have a space for every possible month so that the values 
    can be separated by the correct time periods'''

    ### Assign new arrays ###
    spline_corr_test_j_flux = vari_funcs.cross_correlation.make_corr_arrays(spline_test_j_flux, 
                                                                     jmask, 
                                                                     full_months)
    spline_corr_test_k_flux = vari_funcs.cross_correlation.make_corr_arrays(spline_test_k_flux, 
                                                                     kmask,
                                                                     full_months)

    #%% Calculate the CCF at various tau values ###
    spline_out = np.array([vari_funcs.cross_correlation.cross_correlate_de_z(
            spline_corr_test_k_flux, spline_corr_test_j_flux, tau, z, type='dcf') for tau in tau_arr])

    ### Unpack values ###
    spline_ccf = spline_out[:,0]
    spline_ccf_err = spline_out[:,1]
    spline_pair_num = spline_out[:,2]
    
    #%% Find weighted mean and skew of ccf ###
    spline_max_lag[n] = tau_arr[np.argmax(spline_ccf)]
    
    #%% Fit a parabola for those points around the centre of the ccf function ###
    spline_test_ccf = spline_ccf[np.isin(tau_arr, sub_tau)]
    fit_params, pcov = scipy.optimize.curve_fit(parabola, sub_tau, spline_test_ccf)
    spline_ccf_fit = parabola(plot_tau, *fit_params)
    spline_max_lag_fit[n] = plot_tau[np.argmax(spline_ccf_fit)]
    
    #%% Plot Spline CCF and fit ###
    plt.figure(7)
    plt.subplot(3,2,n+1)
    plt.errorbar(tau_arr, spline_ccf, yerr=spline_ccf_err, fmt='o', zorder=0)
    plt.plot(plot_tau, spline_ccf_fit, zorder=1)
    plt.vlines(spline_max_lag[n], np.min(spline_ccf), np.max(spline_ccf),linestyle='dashed', 
               label='max = '+str(spline_max_lag[n]))
    plt.vlines(spline_max_lag_fit[n], np.min(spline_ccf), np.max(spline_ccf),linestyle='dashdot', 
               label='fitted max = '+str(round(spline_max_lag_fit[n],2)))
    plt.legend()
    plt.xlabel('Lag (months)')
    plt.ylabel('CCF')
    plt.title(labels[n])
    
    #%% Plot Number of Pairs ###
    plt.figure(9)
    plt.subplot(3,2,n+1)
    plt.plot(tau_arr, spline_pair_num, 'o')
    plt.legend()
    plt.xlabel('Lag (months)')
    plt.ylabel('Number of Pairs')
    plt.title(labels[n])
    
    #%% Make plots ###
    plt.figure(4)
    plt.text(xcoords[n], ycoords[n], str(spline_max_lag[n]))
    plt.figure(5)
    plt.text(xcoords[n], ycoords[n], str(round(spline_max_lag_fit[n],2)))

plt.figure(6)
plt.tight_layout()
plt.subplots_adjust(top=0.9)

plt.figure(7)
plt.tight_layout()
plt.subplots_adjust(top=0.9)
    
end = time.time()
print(end-start)
    
    

    
    



