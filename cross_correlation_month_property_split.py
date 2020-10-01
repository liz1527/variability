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
n = input('Mass(M) or Redshift (Z) bins? ')
if n == 'M':
    key = 'Mstar_z_p' # column to split on
    unit = '$M_{\odot}$' #for label
elif n=='Z':
    key = 'z_use' # column to split on
    unit = '' #for label
else:
    print('Invalid Key Entered')
    
num_bins = int(input('Number of bins: '))#3 # how many bins to create
tau_arr = np.arange(-24,25) # tau values to evaluate ccf at 
#tau_arr = np.arange(-10,10)
log = input('log scale? (True or False) ')
if log == 'True':
    log = True
else:
    log = False
save = input('Save plots? (True or False) ')
if save == 'True':
    save = True
else:
    save = False

if save == True:
    subfolder = input('subfolder? (leave blank if none) ')
    if key == 'z_p' or key == 'z' or key == 'z_use':
        filepath = 'plots/new_catalogue/JK_lightcurve_comp/dcf/z_split/'+subfolder+'/'
    elif key == 'Mstar_z_p':
        filepath = 'plots/new_catalogue/JK_lightcurve_comp/dcf/mass_split/'+subfolder+'/'

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
varydata.sort(keys=key) #sort into increasing according to key
varydata = varydata[varydata[key]!=0] # remove any with 0 as this is likely to be invalid
#varydata = varydata[varydata['z_use']<4] # remove any z>4 as this is likely to be invalid
#varydata = varydata[varydata[key]<4] # remove any >4 as this is likely to be invalid

#xvarydata = varydata[varydata['X-ray']==True]
#noxvarydata = varydata[varydata['X-ray']==False]

##%% Set up Mstar z plot so see bin ranges ###
#plt.figure(20)
#plt.plot(varydata['z_p'], varydata['Mstar_z_p'], 'o')
#plt.yscale('log')
#plt.ylabel('Stellar Mass')
#plt.xlabel('Redshift')

#%% Limit mass range for z study or vice versa ###
if key == 'z_p' or key == 'z' or key == 'z_use':
    user_max = input('Enter Maximum Stellar Mass: ')
    user_min = input('Enter Minimum Stellar Mass: ')
    varydata = varydata[varydata['Mstar_z_p']<float(user_max)] #1e11
    varydata = varydata[varydata['Mstar_z_p']>float(user_min)] #5e10
elif key == 'Mstar_z_p':
    user_max = input('Enter Maximum redshift: ')
    user_min = input('Enter Minimum redshift: ')
    varydata = varydata[varydata['z_p']<float(user_max)] #3.5
    varydata = varydata[varydata['z_p']>float(user_min)] #2.5
    

#%% Split data into bins by number in bin ###
num = int(round(len(varydata)/num_bins,0)) # calculate number in each bin
if log == True:
    bins = np.logspace(np.log10(varydata[key][0]),
                       np.log10(varydata[key][-1]),30)
else:    
    bins = np.linspace(varydata[key][0], varydata[key][-1],30)

### Set up arrays for mean vs mass plot ###
mean_lag = np.zeros(num_bins)
mean_norm_lag = np.zeros(num_bins)
mean_lag_err = np.zeros(num_bins)
mean_norm_lag_err = np.zeros(num_bins)
median_lag = np.zeros(num_bins)
median_norm_lag = np.zeros(num_bins)
ccf_skew = np.zeros(num_bins)
norm_ccf_skew = np.zeros(num_bins)
all_bin_edges = np.zeros([2, num_bins])
all_bin_mean = np.zeros(num_bins)
max_lag = np.zeros(num_bins)
max_norm_lag = np.zeros(num_bins)
max_lag_fit = np.zeros(num_bins)
max_norm_lag_fit = np.zeros(num_bins)

for n in range(num_bins):
    bindata = varydata[n*num:(n+1)*num] # define bin according to number in bin
    print(len(bindata))
    ### get bin edges ###
    all_bin_edges[0,n] = np.min(bindata[key])
    all_bin_edges[1,n] = np.max(bindata[key])
    all_bin_mean[n] = np.mean(bindata[key])
    bin_min = round(all_bin_edges[0,n],2)
    bin_max = round(all_bin_edges[1,n],2)
    label='{:.2e}'.format(bin_min)+unit+r' $\leq$ '+key+r' $<$ '+'{:.2e}'.format(bin_max)+unit
    plt.figure(20)
    plt.hlines(bin_min, 0, 6, color='C'+str(n))
    plt.hlines(bin_max, 0, 6, color='C'+str(n))
    
    plt.figure(1)
    plt.hist(bindata[key], bins, histtype='step')
    if log == True:
        plt.xscale('log')   
    plt.xlabel(key)
    plt.tight_layout()
    if save==True:
        if key == 'z_p' or key == 'z' or key == 'z_use':
            plt.savefig(filepath+'all_z_hist.png', overwrite='True')
        elif key == 'Mstar_z_p':
            plt.savefig(filepath+'all_mass_hist.png', overwrite='True')

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
        
    if dez == 'Y' or dez == 'y':
        ### do ccf with deredshifted light curves ####
        z = bindata['z_use']    
        out = np.array([vari_funcs.cross_correlation.cross_correlate_de_z(
                corr_test_k_flux, corr_test_j_flux, tau, z, type='dcf') for tau in tau_arr])
    else:
        out = np.array([vari_funcs.cross_correlation.cross_correlate(
                corr_test_k_flux, corr_test_j_flux, tau, type='dcf') for tau in tau_arr])

    ### Unpack values ###
    ccf = out[:,0]
    ccf_err = out[:,1]

    #%% Find weighted mean and skew of ccf ###
    mean_lag[n], mean_lag_err[n] = weighted_mean_and_err(tau_arr, ccf)
    median_lag[n] = ws.numpy_weighted_median(tau_arr, weights=ccf)
    ccf_skew[n] = skew(ccf)
    max_lag[n] = tau_arr[np.argmax(ccf)]
    
    #%% Fit a parabola for those points around the centre of the ccf function ###
    sub_tau = np.arange(-5,6)
    test_ccf = ccf[np.isin(tau_arr, sub_tau)]
    fit_params, pcov = scipy.optimize.curve_fit(parabola, sub_tau, test_ccf)
    plot_tau = np.linspace(-5,6, 30)
    ccf_fit = parabola(plot_tau, *fit_params)
    max_lag_fit[n] = plot_tau[np.argmax(ccf_fit)]
    
    #%% Make plots ###
    plt.figure(2,figsize=[10,10])
    #plt.subplot(211)
    #plt.plot(tau_arr, ccf,'o')
    plt.errorbar(tau_arr, ccf, yerr=ccf_err, fmt='o', color = 'C'+str(n),
                 label=label)
#    plt.vlines(mean_lag[n], -0.015,0.02, color='C'+str(n), linestyle='dashed')
#    plt.vlines(median_lag[n], -0.015,0.02, color='C'+str(n), linestyle='dotted')
#    plt.plot(plot_tau, ccf_fit, 'C'+str(n))
    plt.xlabel('Lag (months)')
    plt.ylabel('Cross-Correlation Function')
    plt.ylim(-0.7,0.9)
    plt.grid(True)
    plt.legend(loc='lower center')
    plt.tight_layout()
    if save==True:
        if key == 'z_p' or key == 'z' or key == 'z_use':
            plt.savefig(filepath+'all_z.png', overwrite='True')
        elif key == 'Mstar_z_p':
            plt.savefig(filepath+'all_mass.png', overwrite='True')
    
    plt.figure(40,figsize=[10,10])
    #plt.subplot(211)
    #plt.plot(tau_arr, ccf,'o')
    plt.vlines(max_lag_fit[n], -0.015,0.02, color='C'+str(n), linestyle='dashed')
#    plt.vlines(median_lag[n], -0.015,0.02, color='C'+str(n), linestyle='dotted')
    plt.plot(plot_tau, ccf_fit, 'C'+str(n),label=label)
    plt.xlabel('Lag (months)')
    plt.ylabel('Cross-Correlation Function')
    plt.xlim(-25,25)
    plt.ylim(-0.7,0.9)
    plt.grid(True)
    plt.legend(loc='lower center')
    plt.tight_layout()
    if save==True:
        if key == 'z_p' or key == 'z':
            plt.savefig(filepath+'fit_z.png', overwrite='True')
        elif key == 'Mstar_z_p':
            plt.savefig(filepath+'fit_mass.png', overwrite='True')
    
    plt.figure(n+4,figsize=[10,10])
    #plt.subplot(211)
    #plt.plot(tau_arr, ccf,'o')
    plt.errorbar(tau_arr, ccf, yerr=ccf_err, fmt='o', color='C'+str(n),
                 label=label)
    plt.vlines(mean_lag[n], -0.015,0.02, color='C'+str(n), linestyle='dashed',
               label='Mean = '+str(round(mean_lag[n],2)))
    plt.vlines(median_lag[n], -0.015,0.02, color='C'+str(n), linestyle='dotted',
               label='Median = '+str(round(median_lag[n],2)))
    plt.plot(plot_tau, ccf_fit, 'C'+str(n))
    plt.xlabel('Lag (months)')
    plt.ylabel('Cross-Correlation Function')
    plt.ylim(-0.7,0.9)
    plt.grid(True)
    plt.legend(loc='lower center')
    plt.tight_layout()
    if save==True:
        if key == 'z_p' or key == 'z' or key == 'z_use':
            plt.savefig(filepath+str(n)+'_z.png', overwrite='True')
        elif key == 'Mstar_z_p':
            plt.savefig(filepath+str(n)+'_mass.png', overwrite='True')
    
    #%% Normalise so same area under each ccf ###
    area = np.trapz(x=tau_arr, y=ccf)
    norm_ccf = ccf/area
    norm_ccf_err = ccf_err/area

    #%% Find weighted mean and skew of normalised ccf ###
    mean_norm_lag[n], mean_norm_lag_err[n] = weighted_mean_and_err(tau_arr, ccf)
    median_norm_lag[n] = ws.numpy_weighted_median(tau_arr, weights=ccf)
    norm_ccf_skew[n] = skew(norm_ccf)
    max_norm_lag[n] = tau_arr[np.argmax(norm_ccf)]
    
    #%% Fit a parabola for those points around the centre of the ccf function ###
    sub_tau = np.arange(-5,6)
    test_ccf = norm_ccf[np.isin(tau_arr, sub_tau)]
    fit_params, pcov = scipy.optimize.curve_fit(parabola, sub_tau, test_ccf)
    plot_tau = np.linspace(-5,6,100)
    ccf_fit = parabola(plot_tau, *fit_params)
    max_norm_lag_fit[n] = plot_tau[np.argmax(ccf_fit)]
    
    #%% Create plots ###
    plt.figure(3,figsize=[10,10])
    #plt.subplot(211)
    #plt.plot(tau_arr, ccf,'o')
    plt.errorbar(tau_arr, norm_ccf, yerr=norm_ccf_err, fmt='o',
                 label=label)
#    plt.vlines(mean_norm_lag[n], -0.05,0.125, color='C'+str(n), linestyle='dashed')
#    plt.vlines(median_lag[n], -0.05,0.125, color='C'+str(n), linestyle='dotted')
    plt.plot(plot_tau, ccf_fit, 'C'+str(n))
    plt.xlabel('Lag (months)')
    plt.ylabel('Normalised Cross-Correlation Function')
    plt.ylim(-0.05,0.125)
    plt.grid(True)
    plt.legend(loc='lower center')
    plt.tight_layout()
    if save==True:
        if key == 'z_p' or key == 'z' or key == 'z_use':
            plt.savefig(filepath+'normalised/all_z.png', overwrite='True')
        elif key == 'Mstar_z_p':
            plt.savefig(filepath+'normalised/all_mass.png', overwrite='True')
#    
    plt.figure(41,figsize=[10,10])
    #plt.subplot(211)
    #plt.plot(tau_arr, ccf,'o')
    plt.vlines(max_norm_lag_fit[n], -0.05,0.125, color='C'+str(n), linestyle='dashed')
#    plt.vlines(median_lag[n], -0.015,0.02, color='C'+str(n), linestyle='dotted')
    plt.plot(plot_tau, ccf_fit, 'C'+str(n),label=label)
    plt.xlabel('Lag (months)')
    plt.ylabel('Normalised Cross-Correlation Function')
    plt.xlim(-25,25)
    plt.ylim(-0.05,0.125)
    plt.grid(True)
    plt.legend(loc='lower center')
    plt.tight_layout()
    if save==True:
        if key == 'z_p' or key == 'z' or key == 'z_use':
            plt.savefig(filepath+'normalised/fit_z.png', overwrite='True')
        elif key == 'Mstar_z_p':
            plt.savefig(filepath+'normalised/fit_mass.png', overwrite='True')
            
    plt.figure(n+4+num_bins,figsize=[10,10])
    #plt.subplot(211)
    #plt.plot(tau_arr, ccf,'o')
    plt.errorbar(tau_arr, norm_ccf, yerr=norm_ccf_err, fmt='o', color='C'+str(n),
                 label=label)
    plt.vlines(mean_norm_lag[n], -0.05,0.125, color='C'+str(n), linestyle='dashed',
               label='Mean = '+str(round(mean_norm_lag[n],2)))
    plt.vlines(median_lag[n], -0.05,0.125, color='C'+str(n), linestyle='dotted',
               label='Median = '+str(round(median_norm_lag[n],2)))
    plt.plot(plot_tau, ccf_fit, 'C'+str(n))
    plt.xlabel('Lag (months)')
    plt.ylabel('Cross-Correlation Function')
    plt.ylim(-0.05,0.125)
    plt.grid(True)
    plt.legend(loc='lower center')
    plt.tight_layout()
    if save==True:
        if key == 'z_p' or key == 'z' or key == 'z_use':
            plt.savefig(filepath+'/normalised/'+str(n)+'_z.png', overwrite='True')
        elif key == 'Mstar_z_p':
            plt.savefig(filepath+'/normalised/'+str(n)+'_mass.png', overwrite='True')
    
#%% Plot mean lag vs mass ###
all_bin_errors = np.abs(all_bin_edges - all_bin_mean[None,:])
plt.figure()
#plt.errorbar(all_bin_mean, mean_norm_lag, xerr=all_bin_errors, fmt='o')
plt.errorbar(all_bin_mean, mean_lag, xerr=all_bin_errors, yerr=mean_lag_err, fmt='o')
plt.xlabel(key)
plt.ylabel('Weighted Mean Lag (months)')
if log==True:
    plt.xscale('log')
plt.tight_layout()
if save == True:
    if key == 'z_p' or key == 'z' or key == 'z_use':
        plt.savefig(filepath+'/mean_lag_vs_z.png')
    elif key == 'Mstar_z_p':
        plt.savefig(filepath+'/mean_lag_vs_mass.png')
#
#plt.figure()
#plt.errorbar(mean_norm_lag, all_bin_mean, yerr=all_bin_errors, fmt='o')
#plt.ylabel(key)
#plt.xlabel('Weighted Mean Lag (months)')
#if log==True:
#    plt.yscale('log')
#plt.tight_layout()

#%% Plot median lag vs mass ###
plt.figure()
#plt.errorbar(all_bin_mean, median_norm_lag, xerr=all_bin_errors, fmt='o')
plt.errorbar(all_bin_mean, median_lag, xerr=all_bin_errors, fmt='o')
plt.xlabel(key)
plt.ylabel('Weighted Median Lag (months)')
if log==True:
    plt.xscale('log')
plt.tight_layout()
if save == True:
    if key == 'z_p' or key == 'z' or key == 'z_use':
        plt.savefig(filepath+'/median_lag_vs_z.png')
    elif key == 'Mstar_z_p':
        plt.savefig(filepath+'/median_lag_vs_mass.png')

#plt.figure()
#plt.errorbar(median_norm_lag, all_bin_mean, yerr=all_bin_errors, fmt='o')
#plt.ylabel(key)
#plt.xlabel('Weighted Median Lag (months)')
#if log==True:
#    plt.yscale('log')
#plt.tight_layout()


#%% Plot max lag vs mass ###
plt.figure()
#plt.errorbar(all_bin_mean, max_norm_lag, xerr=all_bin_errors, fmt='o')
plt.errorbar(all_bin_mean, max_lag, xerr=all_bin_errors, fmt='o')
plt.xlabel(key)
plt.ylabel('Max Lag (months)')
if log==True:
    plt.xscale('log')
plt.tight_layout()
if save == True:
    if key == 'z_p' or key == 'z' or key == 'z_use':
        plt.savefig(filepath+'/max_lag_vs_z.png')
    elif key == 'Mstar_z_p':
        plt.savefig(filepath+'/max_lag_vs_mass.png')

#%% Plot max lag vs mass ###
plt.figure()
#plt.errorbar(all_bin_mean, max_norm_lag_fit, xerr=all_bin_errors, fmt='o')
plt.errorbar(all_bin_mean, max_lag_fit, xerr=all_bin_errors, fmt='o')
plt.xlabel(key)
plt.ylabel('Max Fitted Lag (months)')
if log==True:
    plt.xscale('log')
plt.tight_layout()
if save == True:
    if key == 'z_p' or key == 'z' or key == 'z_use':
        plt.savefig(filepath+'/max_lag_fit_vs_z.png')
    elif key == 'Mstar_z_p':
        plt.savefig(filepath+'/max_lag_fit_vs_mass.png')
        
end = time.time()
print(end-start)
    
    

    
    



