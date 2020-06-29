#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  1 13:30:42 2018

Code to compute maximum likelihood variability

@author: ppxee
"""
import time
start = time.time()
print(start)

### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
from astropy.table import Table, Column #for handling tables
import numpy as np #for handling arrays
#import math
#from astropy.stats import median_absolute_deviation
import vari_funcs #my module to help run code neatly
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
plt.close('all') #close any open plots
#from numpy.lib.recfunctions import append_fields

#%% Open the fits files and get data ###
### Import data for variables selected in J and K ###
varydata = Table.read('variable_tables/J_and_K_variables_month_varystats_DR11data.fits')
xvarydata = varydata[varydata['X-ray']==True]
noxvarydata = varydata[varydata['X-ray']==False]

### Import matched sample ###
#xvarydata = Table.read('variable_tables/matched_Xray_J_and_K_variables_varystats_DR11data_0.25_0.5_0.2.fits')
#noxvarydata = Table.read('variable_tables/matched_notXray_J_and_K_variables_varystats_DR11data_0.25_0.5_0.2.fits')

#%% Get lightcurves with errors ###
def getdata(data):
    ### Set up month arrays ###
    kmonths = ['sep05','oct05','nov05','dec05', 'jan06', #'dec06', 
              'jan07', 'aug07', 'sep07', 'oct07', 'sep08', 'oct08', 'nov08', 
              'jul09', 'aug09', 'sep09', 'oct09', 'nov09', 'dec09', 'jan10', 
              'feb10', 'aug10', 'sep10', 'oct10', 'nov10', 'dec10', 'jan11', #'feb11', 
              'aug11', 'sep11', 'oct11', 'nov11', 'dec11', 'jan12', 'feb12', 
              'jul12', 'aug12', 'sep12', 'oct12', 'nov12']
    
    jmonths = ['sep05', 'oct05', 'nov05', 'dec05', 'jan06', 'oct06', 'nov06',
              'dec06', 'aug07', 'sep07', 'oct07', 'oct08', 'nov08', 'aug09',
              'sep09', 'oct09', 'nov09', 'dec09', 'aug10', 'sep10', 'oct10',
              'nov10', 'dec10', 'jan11', 'aug11', 'sep11', 'oct11', 'nov11',
              'dec11', 'jan12', 'jul12', 'aug12', 'sep12', 'oct12', 'nov12']
    
    jmask = np.isin(jmonths, kmonths)
    kmask = np.isin(kmonths, jmonths)
    j_flux = data['Month_Flux_J'][:,jmask]
    j_fluxerr = data['Month_Fluxerr_J'][:,jmask]
    k_flux = data['Month_Flux_K'][:,kmask]
    k_fluxerr = data['Month_Fluxerr_K'][:,kmask]
    
    j_flux, j_fluxerr = vari_funcs.flux_funcs.normalise_flux_and_errors(j_flux, 
                                                                        j_fluxerr)
    k_flux, k_fluxerr = vari_funcs.flux_funcs.normalise_flux_and_errors(k_flux, 
                                                                        k_fluxerr)
    return j_flux, j_fluxerr, k_flux, k_fluxerr

j_flux, j_fluxerr, k_flux, k_fluxerr = getdata(varydata)
x_j_flux, x_j_fluxerr, x_k_flux, x_k_fluxerr = getdata(xvarydata)
nox_j_flux, nox_j_fluxerr, nox_k_flux, nox_k_fluxerr = getdata(noxvarydata)

#%% Find chi in relation to each other ###

def my_chisquare_exp_err(flux, fluxerr, expcurve, experr):
    '''Function that finds the chi square value of each light curve when passed
    an array of flux values, corresponding array of flux errors, and and expected
    lightcurve. 
    It assumes each row is a single light curves and tests against the null 
    hypothesis that the light curve is the same as the expected light curve given.
    The errors on the light curve being tested are taken into account.
    Inputs:
        flux = a 2D array of flux values where each row is a lightcurve for a 
               single object
        fluxerr = a 2D array of flux errors that correspond to the flux array
        expcurve = a 2D array of expected values for the light curves. 
        experr = a 2D array of expected error values for the light curves. 
    Outputs:
        chi = a 1D array of chi sq values for each row in the flux arrays
    '''
#    fluxn, fluxerrn = vari_funcs.normalise_flux_and_errors(flux, fluxerr)
    top = np.square(flux-expcurve)
    bot = np.square(fluxerr)+np.square(experr)
    chi = np.nansum(top/bot, axis=1)
    return chi

chi_JK = my_chisquare_exp_err(k_flux, k_fluxerr, j_flux, j_fluxerr)

x_chi_JK = my_chisquare_exp_err(x_k_flux, x_k_fluxerr, x_j_flux, x_j_fluxerr)

nox_chi_JK =  my_chisquare_exp_err(nox_k_flux, nox_k_fluxerr, nox_j_flux, nox_j_fluxerr)
    
#%% Find least squares in relation to each other ###

def least_squares(flux1, flux2):
    '''Function that finds the least squares between two lightcurves
    Inputs:
        flux1 = a 2D array of flux values where each row is a lightcurve for a 
               single object
        flux2 = a 2D array of flux values where each row is a lightcurve for a 
               single object
    Outputs:
        leastsq = a 1D array of least sq values for each row in the flux arrays
    '''
    leastsq = np.nansum(np.square(flux1 - flux2),axis=1)

    return leastsq

ls_J = least_squares(j_flux, k_flux)
ls_K = least_squares(k_flux, j_flux)

x_ls_J = least_squares(x_j_flux, x_k_flux)
x_ls_K = least_squares(x_k_flux, x_j_flux)

nox_ls_J = least_squares(nox_j_flux, nox_k_flux)
nox_ls_K = least_squares(nox_k_flux, nox_j_flux)
    
    
#%% plot chi_K hist ###
plt.figure()
#bins = np.logspace(np.log(np.min(chi_J)), np.log(np.max(chi_K)), 25)
bins = np.linspace(np.min(chi_JK), np.max(chi_JK), 25)
plt.hist([x_chi_JK, nox_chi_JK], bins=bins,  color=['r','b'], histtype='step',
         density=True)
plt.xlabel('$\chi^{2} = \sum((F_{K}-F_{J})^{2}/(\sigma_{K}^{2}+\sigma_{J}^{2})))$')
#plt.xscale('symlog')
#plt.yscale('log')
#x = np.linspace(min(chi_J), max(chi_J), 10)
#y = x
#plt.plot(x, y, 'k')
plt.tight_layout()
    
#%% Plot chi vs average flux ###
def get_avg_flux(data):
    j_flux = data['Flux_J']
    j_avgflux = np.nanmean(j_flux, axis=1)
    k_flux = data['Flux_K']
    k_avgflux = np.nanmean(k_flux, axis=1)
    
    return j_avgflux, k_avgflux

j_avgflux, k_avgflux = get_avg_flux(varydata)
x_j_avgflux, x_k_avgflux = get_avg_flux(xvarydata)
nox_j_avgflux, nox_k_avgflux = get_avg_flux(noxvarydata)


plt.figure()
plt.plot(nox_k_avgflux, nox_chi_JK,'bo')
plt.plot(x_k_avgflux, x_chi_JK,'rs')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Average K-band Flux')
plt.ylabel('$\chi^{2} = \sum((F_{K}-F_{J})^{2}/(\sigma_{K}^{2}+\sigma_{J}^{2})))$')
plt.hlines(30, 1e2, 5e5)
plt.tight_layout()
    
#%% Bin and plot average chi_JK vs flux ###

bins = np.arange(16,25,0.5)
#bins = np.append(bins, [26,27,28])

bins = 10**((30-bins)/2.5)
bins = np.flip(bins, axis=0)
#bins = bins[16:44] #because of flux limit

### Bin data ###
mean_chi_JK = np.zeros([len(bins)-1])
std_chi_JK = np.zeros([len(bins)-1])
stderr_chi_JK = np.zeros([len(bins)-1])
for n, binedge in enumerate(bins):
#    print(binedge)
    if n==np.size(bins)-1:
        break
    
    ### Define bin mask
    bina = k_avgflux >= binedge
    binb = k_avgflux < bins[n+1]
    bin = bina*binb
    
    ### find mean chi_JK ###
    bin_chi_JK = chi_JK[bin]
    mean_chi_JK[n] = np.nanmean(bin_chi_JK)
    std_chi_JK[n] = np.nanstd(bin_chi_JK)
    stderr_chi_JK[n] = np.nanstd(bin_chi_JK)/np.sqrt(len(bin_chi_JK))
    

plt.figure()
plt.errorbar(bins[:-1], mean_chi_JK, yerr=std_chi_JK, fmt='ko')
plt.xlabel('Mean K-Band Flux')
plt.ylabel('Mean $\chi^{2}_{JK}$')
plt.xscale('log')
plt.tight_layout()



##%% Get IDs of those with high chi ###
#
#high = varydata[chi_JK>30]
#chi_high = chi_JK[chi_JK>30]
#
#new_col = Column(chi_high, name='Chi_JK')
#### add to tbdata ###
#high.add_column(new_col)
#
#high.write('variable_tables/J_and_K_high_chi_JK_variables.fits')

##%% plot difference ###
#plt.figure()
##plt.plot(chi_J, chi_K, 'o')
#plt.plot(nox_ls_J, nox_ls_K, 'bo')
#plt.plot(x_ls_J, x_ls_K, 'rs')
#plt.xlabel('$Least Square = \sum(F_{J}-F_{K})^{2})$')
#plt.ylabel('$Least Square = \sum(F_{K}-F_{J})^{2})$')
##plt.xscale('log')
##plt.yscale('log')
#x = np.linspace(min(ls_J), max(ls_J), 10)
#y = x
#plt.plot(x, y, 'k')
#plt.tight_layout()
#    
##%% plot ls hist ###
#plt.figure()
##bins = np.logspace(np.log(np.min(chi_J)), np.log(np.max(chi_K)), 25)
##bins = np.linspace(np.min(ls_J), np.max(ls_K), 25)
#plt.hist([x_ls_J, nox_ls_J], bins=25,  color=['r','b'], histtype='step')#,
##         density=True)
#plt.xlabel('$Least Square = \sum(F_{J}-F_{K})^{2})$')
##plt.xscale('symlog')
##plt.yscale('log')
##x = np.linspace(min(chi_J), max(chi_J), 10)
##y = x
##plt.plot(x, y, 'k')
#plt.tight_layout()
#    
#    
##%% Plot ls vs average flux ###
##def get_avg_flux(data):
##    j_flux = data['Flux_J']
##    j_avgflux = np.nanmean(j_flux, axis=1)
##    k_flux = data['Flux_K']
##    k_avgflux = np.nanmean(k_flux, axis=1)
##    
##    return j_avgflux, k_avgflux
##
##j_avgflux, k_avgflux = get_avg_flux(varydata)
##x_j_avgflux, x_k_avgflux = get_avg_flux(xvarydata)
##nox_j_avgflux, nox_k_avgflux = get_avg_flux(noxvarydata)
#
#
#plt.figure()
#plt.plot(nox_k_avgflux, nox_ls_J,'bo')
#plt.plot(x_k_avgflux, x_ls_J,'rs')
#plt.xscale('log')
#plt.xlabel('Average K-band Flux')
#plt.ylabel('$Least Square = \sum(F_{J}-F_{K})^{2})$')
    
    
    
    
    
end = time.time()
print(end-start)












