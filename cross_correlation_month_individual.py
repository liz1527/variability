#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 13:58:17 2020

Code to run a individual cross correlation on J and K month lightcurves

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
import matplotlib.colors as colors
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

def prewhiten(flux, bindata, months):
    ### Get semester data ###
    model_j_flux = bindata['Flux_J']
    model_j_fluxerr = bindata['Fluxerr_J']
    
    ### Normalise ###
    model_j_flux = vari_funcs.cross_correlation.weighted_mean_subtract_normalise(
            model_j_flux, model_j_fluxerr)
    
    ### Subtract to get residuals ###
    newflux = np.zeros(np.shape(flux))
    for n, mon in enumerate(months):
        if mon in ('jul05', 'aug05', 'sep05','oct05','nov05','dec05', 'jan06', 'feb06'):
            ### This means the month is in semester 05B ##
            newflux[:,n] = flux[:,n] - model_j_flux[:,0]
        elif mon in ('jul06', 'aug06', 'sep06','oct06','nov06','dec06', 'jan07', 'feb07'):
            ### This means the month is in semester 06B ##
            newflux[:,n] = flux[:,n] - model_j_flux[:,1]
        elif mon in ('jul07', 'aug07', 'sep07','oct07','nov07','dec07', 'jan08', 'feb08'):
            ### This means the month is in semester 07B ##
            newflux[:,n] = flux[:,n] - model_j_flux[:,2]
        elif mon in ('jul08', 'aug08', 'sep08','oct08','nov08','dec08', 'jan09', 'feb09'):
            ### This means the month is in semester 08B ##
            newflux[:,n] = flux[:,n] - model_j_flux[:,3]
        elif mon in ('jul09', 'aug09', 'sep09','oct09','nov09','dec09', 'jan10', 'feb10'):
            ### This means the month is in semester 09B ##
            newflux[:,n] = flux[:,n] - model_j_flux[:,4]
        elif mon in ('jul10', 'aug10', 'sep10','oct10','nov10','dec10', 'jan11', 'feb11'):
            ### This means the month is in semester 10B ##
            newflux[:,n] = flux[:,n] - model_j_flux[:,5]
        elif mon in ('jul11', 'aug11', 'sep11','oct11','nov11','dec11', 'jan12', 'feb12'):
            ### This means the month is in semester 11B ##
            newflux[:,n] = flux[:,n] - model_j_flux[:,6]
        elif mon in ('jul12', 'aug12', 'sep12','oct12','nov12','dec12', 'jan13', 'feb13'):
            ### This means the month is in semester 12B ##
            newflux[:,n] = flux[:,n] - model_j_flux[:,7]
            
    return newflux, model_j_flux

#%% Import data for variables selected in J and K + define constants ###
varydata = Table.read('variable_tables/J_and_K_variables_month_varystats_DR11data.fits')
#varydata = varydata[varydata['X-ray']==True]
#varydata = varydata[varydata['ID']==250411]

### Set constants for sample and bins ###
min_chi = 100 # chi cut required to create sample
    
tau_arr = np.arange(-24,25) # tau values to evaluate ccf at 
#tau_arr = np.arange(-10,10)

type = input('DCF or CCF? ')
type = type.lower()

save = input('Save plots? (Y or N) ')
if save == 'Y' or save == 'y':
    save = True
else:
    save = False

if save == True and type == 'ccf':
    filepath = 'plots/new_catalogue/JK_lightcurve_comp/cross-correlation/individual/prewhiten/'
elif save == True and type == 'dcf':
    filepath = 'plots/new_catalogue/JK_lightcurve_comp/dcf/individual/prewhiten/'

save_ccf = input('Create and save CCF plots? (Y or N) ')
if save_ccf == 'Y'or save_ccf == 'y':
    save_ccf = True
else:
    save_ccf = False


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

testmonths = ['sep05','sep06', 'sep07', 'sep08', 'sep09', 'sep10', 'sep11', 'sep12']
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

#%% Create sample and x/nox subsamples ###
mask = varydata['Chi_K'] > min_chi
varydata = varydata[mask]

#%% Run ccf on individual sources ###
num = len(varydata)

### Set up arrays for plots ###
max_lag = np.zeros(num)
mean_lag = np.zeros(num)
mean_lag_test = np.zeros(num)
median_lag = np.zeros(num)
max_ccf = np.zeros(num)
max_lag_fit = np.zeros(num)

#%% Get lightcurves with errors ###
j_flux, j_fluxerr, k_flux, k_fluxerr = getdata(varydata)

#%% Mean subract and normalise ###

test_j_flux = vari_funcs.cross_correlation.weighted_mean_subtract_normalise(
        j_flux, j_fluxerr)
test_k_flux = vari_funcs.cross_correlation.weighted_mean_subtract_normalise(
        k_flux, k_fluxerr)

#plt.figure(figsize=[9,6])
#plt.plot(Jx_months, test_j_flux[15,:], 'o', label='J-band')
#plt.plot(Kx_months, test_k_flux[15,:], 'o', label='K-band')

#%% Prewhiten using J semester data ###
#test_j_flux, model_j_flux = prewhiten(test_j_flux, varydata, jmonths)
#test_k_flux, model_j_flux = prewhiten(test_k_flux, varydata, kmonths)
#plt.plot(testx_months, model_j_flux[15,:], 'o', label='J-band Semester')
#plt.xlabel('Month')
#plt.ylabel('Normalised Flux')
#plt.xticks(inds, month_ticks, rotation = 'vertical')
#plt.legend()
#plt.tight_layout()
#
#plt.figure(figsize=[9,6])
#plt.plot(Jx_months, test_j_flux[15,:], 'o', label='J-band')
#plt.plot(Kx_months, test_k_flux[15,:], 'o', label='K-band')
#plt.xlabel('Month')
#plt.ylabel('Prewhitened Flux')
#plt.xticks(inds, month_ticks, rotation = 'vertical')
#plt.legend()
#plt.tight_layout()
#    break

for n in range(num):
    print(n)
    bindata = varydata[n]
#    
#    #%% Get lightcurves with errors ###
#    j_flux, j_fluxerr, k_flux, k_fluxerr = getdata(bindata)
#
#    #%% Mean subract and normalise ###
#    
#    test_j_flux = vari_funcs.cross_correlation.weighted_mean_subtract_normalise_single(
#            j_flux, j_fluxerr)
#    test_k_flux = vari_funcs.cross_correlation.weighted_mean_subtract_normalise_single(
#            k_flux, k_fluxerr)
#    
##    plt.figure()
##    plt.plot(Jx_months, test_j_flux, 'o')
##    plt.plot(Kx_months, test_k_flux, 'o')
#    
#    #%% Prewhiten using J semester data ###
#    test_j_flux, model_j_flux = prewhiten(test_j_flux, bindata, jmonths)
#    test_k_flux, model_j_flux = prewhiten(test_k_flux, bindata, kmonths)
##    plt.plot(testx_months, model_j_flux, 'o')
#    
##    plt.figure()
##    plt.plot(Jx_months, test_j_flux, 'o')
##    plt.plot(Kx_months, test_k_flux, 'o')
##    break
    #%% Create correlation arrays ###
    ''' Need arrays that have a space for every possible month so that the values 
    can be separated by the correct time periods'''

    ### Assign new arrays ###
    corr_test_j_flux = vari_funcs.cross_correlation.make_corr_arrays_single(test_j_flux[n,:], 
                                                                     jmask, 
                                                                     full_months)
    corr_test_k_flux = vari_funcs.cross_correlation.make_corr_arrays_single(test_k_flux[n,:], 
                                                                     kmask,
                                                                     full_months)

    #%% Calculate the CCF at various tau values ###
    out = np.array([vari_funcs.cross_correlation.cross_correlate(
            corr_test_k_flux, corr_test_j_flux, tau, type) for tau in tau_arr])

    ### Unpack values ###
    ccf = out[:,0]
    ccf_err = out[:,1]

    #%% Calculate the ACF at various tau values for J and K ###
        
    out_k = np.array([vari_funcs.cross_correlation.cross_correlate(
            corr_test_k_flux, corr_test_k_flux, tau, type) for tau in tau_arr])

    out_j = np.array([vari_funcs.cross_correlation.cross_correlate(
            corr_test_j_flux, corr_test_j_flux, tau, type) for tau in tau_arr])
    
    ### Unpack values ###
    acf_k = out_k[:,0]
    acf_k_err = out_k[:,1]
    acf_j = out_j[:,0]
    acf_j_err = out_j[:,1]
    

    #%% Fit a parabola for those points around the centre of the ccf function ###
    sub_tau = np.arange(-7,8)
    test_ccf = ccf[np.isin(tau_arr, sub_tau)]
    fit_params, pcov = scipy.optimize.curve_fit(parabola, sub_tau, test_ccf)
    plot_tau = np.linspace(-7,8,100)
    ccf_fit = parabola(plot_tau, *fit_params)
    max_lag_fit[n] = plot_tau[np.argmax(ccf_fit)]
    
    #%% Find mean lag (centroid) in the centre ###
    mean_lag[n] = ws.numpy_weighted_mean(sub_tau, weights=test_ccf)
    mean_lag_test[n] = np.average(sub_tau, weights=test_ccf)
    median_lag[n] = ws.numpy_weighted_median(sub_tau, weights=test_ccf)
    
    #%% Find lag where ccf is max in centre ###
    max_lag[n] = sub_tau[np.argmax(test_ccf)]
    max_ccf[n] = np.nanmax(test_ccf)
    
    #%% Create and save ccfs if prompted ###
    if save_ccf == True:
        plt.figure(figsize=[8,8])
        plt.title('DR11 ID: '+str(bindata['ID'])+' $\chi^{2}$ = '+str(bindata['Chi_K']))
        plt.errorbar(tau_arr, ccf, yerr=ccf_err, fmt='o', color='k', 
                     label = 'CCF')
        plt.errorbar(tau_arr, acf_k, yerr=acf_k_err, fmt='-', color='r', 
                     label = 'K-band ACF', alpha=0.5)
        plt.errorbar(tau_arr, acf_j, yerr=acf_j_err, fmt='-', color='b', 
                     label = 'J-band ACF', alpha=0.5)
        plt.vlines(max_lag[n], np.min(ccf), np.max(ccf), linestyle='dashed', 
                   label='Max Lag')
#        plt.vlines(max_lag_fit[n], np.min(ccf), np.max(ccf), linestyle='dotted', 
#                   label='Max Fitted Lag')
        plt.vlines(mean_lag[n], np.min(ccf), np.max(ccf), linestyle='dashdot', 
                   label='Mean Lag')
#        plt.plot(plot_tau, ccf_fit)
        plt.xlabel('Lag (months)')
        plt.ylabel('Cross-Correlation Function')
        plt.xlim(-25,25)
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.savefig(filepath+type+'s_with_acfs/'+str(bindata['ID'])+'.png', overwrite=True)
        plt.close('all')
        
#%% Plot max lag vs mass ###
mstar = varydata['Mstar_z_p']
mask = mstar != 0 
plt.figure()
plt.scatter(mstar[mask], max_lag[mask], c=max_ccf[mask])
plt.xlabel('Mstar_z_p')
plt.ylabel('Max CCF Lag (months)')
plt.xscale('log')
cbar = plt.colorbar()
cbar.set_label('Max CCF Value')
plt.tight_layout()
if save == True:
    plt.savefig(filepath+'max_lag_vs_mass.png')

plt.figure()
plt.scatter(mstar[mask], max_lag[mask], c=varydata['Chi_K'][mask],
            norm=colors.LogNorm(vmin=varydata['Chi_K'][mask].min(), 
                                vmax=varydata['Chi_K'][mask].max()))
plt.xlabel('Mstar_z_p')
plt.ylabel('Max CCF Lag (months)')
plt.xscale('log')
cbar = plt.colorbar()
cbar.set_label('Chi_K')
plt.tight_layout()
#%% Plot max lag vs z ###
z = varydata['z_use']
plt.figure()
plt.scatter(z, max_lag, c=max_ccf)
plt.xlabel('z')
plt.ylabel('Max CCF Lag (months)')
#plt.xscale('log')
cbar = plt.colorbar()
cbar.set_label('Max CCF Value')
plt.tight_layout()
if save == True:
    plt.savefig(filepath+'max_lag_vs_z.png')


plt.figure()
plt.scatter(z, max_lag, c=varydata['Chi_K'],
            norm=colors.LogNorm(vmin=varydata['Chi_K'].min(), 
                                vmax=varydata['Chi_K'].max()))
plt.xlabel('z')
plt.ylabel('Max CCF Lag (months)')
#plt.xscale('log')
cbar = plt.colorbar()
cbar.set_label('Chi_K')
plt.tight_layout()

#%% Plot max restricted vs mass and z ###
plt.figure()
plt.plot(z, max_lag,'o')
plt.xlabel('z')
plt.ylabel('Restricted Max Lag (months)')
plt.tight_layout()
if save == True:
    plt.savefig(filepath+'max_lag_res_vs_z.png')

plt.figure()
plt.plot(mstar[mask], max_lag[mask],'o')
plt.xlabel('Mstar_z_p')
plt.xscale('log')
plt.ylabel('Restricted Max Lag (months)')
plt.tight_layout()
if save == True:
    plt.savefig(filepath+'max_lag_res_vs_mass.png')

#%% Plot mean vs mass and z ###
plt.figure()
plt.plot(z, mean_lag,'o')
plt.xlabel('z')
plt.ylabel('Mean Lag (months)')
plt.ylim(-7,7)
plt.tight_layout()
if save == True:
    plt.savefig(filepath+'mean_lag_fit_vs_z.png')

plt.figure()
plt.plot(mstar[mask], mean_lag[mask],'o')
plt.xlabel('Mstar_z_p')
plt.ylabel('Mean Lag (months)')
plt.ylim(-7,7)
plt.tight_layout()
if save == True:
    plt.savefig(filepath+'mean_lag_fit_vs_mass.png')
    
#%% Plot median vs mass and z ###
plt.figure()
plt.plot(z, median_lag,'o')
plt.xlabel('z')
plt.ylabel('Median Lag (months)')
plt.ylim(-7,7)
plt.tight_layout()
if save == True:
    plt.savefig(filepath+'median_lag_fit_vs_z.png')

plt.figure()
plt.plot(mstar[mask], median_lag[mask],'o')
plt.xlabel('Mstar_z_p')
plt.ylabel('Median Lag (months)')
plt.ylim(-7,7)
plt.tight_layout()
if save == True:
    plt.savefig(filepath+'median_lag_fit_vs_mass.png')
    
end = time.time()
print(end-start)
    
    

    
    



