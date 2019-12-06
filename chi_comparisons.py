#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 15:56:54 2019

Code to look at the different chisq distributions for different selection
populations (e.g. neg vs no neg and J vs K)

@author: ppxee
"""

### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
from astropy.table import Table #for handling tables
import numpy as np #for handling arrays
#import math
#from astropy.stats import median_absolute_deviation
from scipy import stats
import vari_funcs #my module to help run code neatly
#from scipy.stats import chisquare
plt.close('all') #close any open plots

def get_chi_and_flux_J(Jdata, Jsigtb, neg=True):
    
    ### Remove edges ###
    Jdata = vari_funcs.field_funcs.remove_edges(Jdata)
    
    ## Create arrays of flux values ###
    Jflux = vari_funcs.j_mag_flux.flux_stacks(Jdata, aper=4)
    
    if neg==False:
        ### remove values that are negative ###
        Jflux, Jdata = vari_funcs.flux_funcs.noneg(Jflux, Jdata)
        
#    flux, bindata = vari_funcs.fluxbin(binedge, bins[n+1], fluxn, qtbdata) #bindata
    
    ### Get error arrays ###
    Jflux, Jfluxerr, Jdata = vari_funcs.j_mag_flux.create_quad_error_array_J(Jsigtb, Jdata, aper=4)

    ### Get chi sq ###
    Jchisq = vari_funcs.vary_stats.my_chisquare_err(Jflux, Jfluxerr)  
    
    return Jflux, Jfluxerr, Jchisq, Jdata
    
def get_chi_and_flux_K(Kdata, Ksigtb, neg=True):
    
    ### Remove edges ###
    Kdata = vari_funcs.field_funcs.remove_edges(Kdata)
    
    ## Create arrays of flux values ###
    Kflux = vari_funcs.k_mag_flux.flux_stacks(Kdata, aper=4)
    
    if neg==False:
        ### remove values that are negative ###
        Kflux, Kdata = vari_funcs.flux_funcs.noneg(Kflux, Kdata)
    
    ### Get error arrays ###
    Kflux, Kfluxerr, Kdata = vari_funcs.k_mag_flux.create_quad_error_array(Ksigtb, Kdata, aper=4)

    ### Get chi sq ###
    Kchisq = vari_funcs.vary_stats.my_chisquare_err(Kflux, Kfluxerr)  
    
    return Kflux, Kfluxerr, Kchisq, Kdata
    
### Open the fits files and get data ###
Jdata = fits.open('mag_flux_tables/J/mag_flux_table_best_J_extra_clean.fits')[1].data
Jsdata = fits.open('mag_flux_tables/J/stars_mag_flux_table_J_extra_clean.fits')[1].data
Jsigtb = Table.read('sigma_tables/quad_epoch_sigma_table_extra_clean_2arcsec_J_noneg.fits')
Jsigtbneg = Table.read('sigma_tables/quad_epoch_sigma_table_extra_clean_2arcsec_J.fits')
Kdata = fits.open('mag_flux_tables/K/mag_flux_table_best_extra_clean_no06.fits')[1].data
Ksdata = fits.open('mag_flux_tables/K/stars_mag_flux_table_extra_clean_no06.fits')[1].data
Ksigtb = Table.read('sigma_tables/quad_epoch_sigma_table_extra_clean_no06_2arcsec.fits')
Ksigtbneg = Table.read('sigma_tables/quad_epoch_sigma_table_extra_clean_no06_2arcsec_neg.fits')

#%% Get chi sq and flux arrays ### 


### for no neg ###
Jflux, Jfluxerr, Jchisq, Jdatanoneg = get_chi_and_flux_J(Jdata, Jsigtb, 
                                                             neg=False)
Kflux, Kfluxerr, Kchisq, Kdatanoneg = get_chi_and_flux_K(Kdata, Ksigtb, 
                                                             neg=False)

### for negatives included
Jfluxneg, Jfluxerrneg, Jchisqneg, Jdataneg = get_chi_and_flux_J(Jdata, 
                                                             Jsigtbneg, 
                                                             neg=True)
Kfluxneg, Kfluxerrneg, Kchisqneg, Kdataneg = get_chi_and_flux_K(Kdata, 
                                                                Ksigtbneg, 
                                                                neg=True)

bins=np.logspace(np.log10(8e-2),np.log10(2e4),50)
#bins=np.logspace(np.log10(3e1),np.log10(2e4),50)
#bins=np.linspace(2e1,2e4,50)

x = np.logspace(np.log10(8e-2),np.log10(2e4),500)
#x=np.logspace(np.log10(3e1),np.log10(2e4),50)
#np.logspace(-2.3,4.4,500)
#x = np.linspace(3e-2,4e4,5000)
y = stats.chi2.pdf(x,6) #6 dof as 7 epochs
y2 = stats.chi2.pdf(x,7) #7 dof as 8 epochs
### Plot J and K ###
plt.figure()
plt.hist([Jchisq, Kchisq], bins, histtype='step', label=['J','K'], density=True)
plt.plot(x,y, label=r'Model $\chi^{2}$ with dof=6')
plt.plot(x,y2, label=r'Model $\chi^{2}$ with dof=7')
plt.xscale('log')
#plt.yscale('log')
plt.xlabel(r'$\chi^{2}$')
plt.legend()
plt.tight_layout()

### Plot J and K neg ###
plt.figure()
plt.hist([Jchisqneg, Kchisqneg], bins, histtype='step', label=['J','K'], density=True)
plt.plot(x,y, label=r'Model $\chi^{2}$ with dof=6')
plt.plot(x,y2, label=r'Model $\chi^{2}$ with dof=7')
plt.xscale('log')
#plt.yscale('log')
plt.xlabel(r'$\chi^{2}$')
plt.legend()
plt.tight_layout()

### Plot J no neg and neg ###
plt.figure()
plt.hist([Jchisq,Jchisqneg], bins, histtype='step', label=['J','J neg'], density=True)
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$\chi^{2}$')
plt.legend()
plt.tight_layout()

### Plot K no neg and neg ###
plt.figure()
plt.hist([Kchisq, Kchisqneg], bins, histtype='step', label=['K','K neg'], density=True)
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$\chi^{2}$')
plt.legend()
plt.tight_layout()