#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 29 12:13:52 2018

Code to look at ensemble variability of objects compared to X-ray luminosity,
with a split based on hardness ratio

@author: ppxee
"""

import time
start = time.time()
#print(start)

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
from scipy.optimize import curve_fit
#plt.close('all') #close any open plots
#from numpy.lib.recfunctions import append_fields

### Define cosmology ###
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
#def remove_edges(tbdata):
#    
#    ### Set X limits ###
#    x = tbdata['X_IMAGE_05B']
#    xmask1 = x > 1000
#    xmask2 = x < 24000
#    xmask = xmask1 * xmask2.astype(bool)
#    tbdata = tbdata[xmask]
#    
#    ### Set Y limits ###
#    y = tbdata['Y_IMAGE_05B']
#    ymask1 = y > 1000
#    ymask2 = y < 24000
#    ymask = ymask1 * ymask2.astype(bool)
#    tbdata = tbdata[ymask]
#    
#    return tbdata
    
#%% Open the fits files and get data ###
chandata = fits.open('variable_tables/no06_variables_chi30_2arcsec_chandata_DR11data_restframe.fits')[1].data
xmmdata = fits.open('variable_tables/no06_variables_chi30_2arcsec_xmmdata_DR11data_restframe.fits')[1].data
fullxray = fits.open('mag_flux_tables/novarys_chanDR11data_restframe_mag_flux_table_best_extra_clean_no06.fits')[1].data
#tbdata = fits.open('variable_tables/no06_variables_chi30_2arcsec_not_deviant_DR11data_restframe.fits')[1].data
tbdata = fits.open('variable_tables/no06_variables_chi30_2arcsec_noXray_DR11data_restframe.fits')[1].data
sndata = fits.open('variable_tables/no06_variables_chi30_2arcsec_DR11data_restframe_SN.fits')[1].data
radiodata = fits.open('variable_tables/no06_variables_chi30_2arcsec_DR11data_restframe_07B.fits')[1].data
#radiodata = fits.open('variable_tables/no06_variables_chi30_radiodata_DR11data_restframe.fits')[1].data
sigtb = Table.read('sigma_tables/quad_epoch_sigma_table_extra_clean_no06_2arcsec.fits')

posvar = np.linspace(0,2,5000)

#%% create function to do things 
def get_luminosity_and_flux(tbdata, xmm=False):
#    
    ### Extract magnitude table and error table ###
    flux = vari_funcs.flux4_stacks(tbdata)
    flux, tbdata = vari_funcs.noneg(flux, tbdata)
#    tbdata = tbdata[np.nanmean(flux,axis=1)>1e4]
    flux, fluxerr, tbdata = vari_funcs.create_quad_error_array(sigtb, tbdata, aper=4)
    
    ### Normalise ###
    fluxnorm, fluxerrnorm = vari_funcs.normalise_flux_and_errors(flux, fluxerr)
    
    ### Find luminosity distance ###
    z = tbdata['z_spec']#[mask]
    z[z==-1] = tbdata['z_p'][z==-1]
    DL = cosmo.luminosity_distance(z)
    DL = DL.to(u.m)
    
    ### Get AB magnitude ###
#    abmag = tbdata['KMAG_20']
    
    ### Convert to luminosity using formula worked out in lab book ###
#    L = 10 ** (-abmag/2.5) * 3631 * 1e-26 * (3e8/0.34e-6) * 4 * np.pi * (DL.value**2) # gives luminosity in W
    L = tbdata['M_K_z_p']
    
    ### remove any that have val of 99 ###
    L[L == 99] = np.nan # remove those with null values
#    L[L > -5] = np.nan # remove faintest as seem spurious
    mask = ~np.isnan(L)
#    [mask]
    return tbdata[mask], L[mask], fluxnorm[mask], fluxerrnorm[mask]


def run_max_likely(tbdata):
    posvar = np.linspace(0,2,5000)
    ### Remove edges ###
    tbdata = vari_funcs.remove_edges(tbdata)
    
    ### Get luminosity and flux ###
    tbdata, L, fluxnorm, fluxerrnorm= get_luminosity_and_flux(tbdata)
    
    ### Get sig values ###
    numobs = np.shape(fluxnorm)[0]
    meanflux = np.nanmean(fluxnorm, axis=1)
    out = np.array([vari_funcs.maximum_likelihood(fluxnorm[n,:], 
                                                  fluxerrnorm[n,:], meanflux[n], 
                                                  posvar, n=n, printn=100) for n in range(numobs)])
    return L, out, tbdata
    
#%% run function ###
def z_split(tbdata):
    z = tbdata['z_spec']#[mask]
    z[z==-1] = tbdata['z_p'][z==-1]
    tb1 = tbdata[z <= 0.3]
    mask = np.array(z>0.5)*np.array(z<=1).astype('bool')
    tb2 = tbdata[mask]
    mask = np.array(z>1)*np.array(z<=1.5).astype('bool')
    tb3 = tbdata[mask]
    mask = np.array(z>1.5)*np.array(z<=2).astype('bool')
    tb4 = tbdata[mask]
    mask = np.array(z>2)*np.array(z<=3).astype('bool')
    tb5 = tbdata[mask]
#    mask = np.array(z>3)*np.array(z<=6).astype('bool')
#    tb6 = tbdata[mask]
    tb6 = tbdata[z > 3]
#    tb2 = tbdata[z>2]
    return tb1, tb2, tb3, tb4, tb5, tb6

def func_pol(x,a,b):
    return a*+b*(x**-1)
def func_exp(x,a,b):
    return a*np.exp(b*x)

tbs = z_split(tbdata)
#L, out, tbdata = run_max_likely(tbdata)

#plt.figure(1, figsize=[12,7])
##z = tbdata['z_spec']#[mask]
##z[z==-1] = tbdata['z_p'][z==-1]
##z_colours = np.copy(z)
##z_colours[z <= 0.5] = 0
##mask = np.array(z>0.5)*np.array(z<=1).astype('bool')
##z_colours[mask] = 0.5
##mask = np.array(z>1)*np.array(z<=1.5).astype('bool')
##z_colours[mask] = 1
##mask = np.array(z>1.5)*np.array(z<=2).astype('bool')
##z_colours[mask] = 1.5
###mask = np.array(z>4)*np.array(z<=5).astype('bool')
###z_colours[mask] = 4
###mask = np.array(z>5)*np.array(z<=6).astype('bool')
###z_colours[mask] = 5
###mask = np.array(z>6)*np.array(z<=7).astype('bool')
###z_colours[mask] = 6
###z_colours[z > 7] = 7
##z_colours[z > 2] = 2
##plt.scatter(L, out[:,0], c=z_colours, zorder=1)
##plt.errorbar(L, out[:,0], yerr=out[:,1], fmt='ko', zorder=0, alpha=0.2)
#plt.xlabel('K Band Absolute Magnitude')
#plt.ylabel(r'$\sigma$')
#plt.title(r'$ax^{2} + bx + c$')
#plt.ylim(ymin=-0.05,ymax=1.5)
#plt.xlim(xmin=-29,xmax=-2)
#plt.gca().invert_xaxis()
##cbar = plt.colorbar()
#plt.legend()
##cbar.set_label('z')
#plt.tight_layout()
#
#plt.figure(2, figsize=[12,7])
##plt.scatter(L, out[:,0], c=z_colours, zorder=1)
##plt.errorbar(L, out[:,0], yerr=out[:,1], fmt='ko', zorder=0, alpha=0.2)
#plt.xlabel('K Band Absolute Magnitude')
#plt.ylabel(r'$\sigma$')
#plt.title('Exponential')
#plt.ylim(ymin=-0.05,ymax=1.5)
#plt.xlim(xmin=-29,xmax=-2)
#plt.gca().invert_xaxis()
##cbar = plt.colorbar()
#plt.legend()
##cbar.set_label('z')
#plt.tight_layout()

labels = ['z<=0.5', '0.5<z<=1','1<z<=1.5','1.5<z<=2','2<z<=3','z>3' ]
plt.figure(figsize=[12,7])
for n, tbdata in enumerate(tbs):
    L, out, tbdata = run_max_likely(tbdata)


    #%% Plot results ###
    
    ### run curve fit ###
    
#    popt, pcov = curve_fit(func_pol, L, out[:,0], sigma=out[:,1])
#    x = np.linspace(-28,-3)
#    y = func_pol(x, popt[0], popt[1],popt[2])
#    
#    perr = np.sqrt(np.diag(pcov))
#    yerr1 = func_pol(x, popt[0]+perr[0], popt[1]-perr[1], popt[2]+perr[2])
#    yerr2 = func_pol(x, popt[0]-perr[0], popt[1]+perr[1], popt[2]-perr[2])
    
#    plt.subplot(2,3,n+1)
#    plt.plot(L, out[:,0], 'o', zorder=1)
#    plt.errorbar(L, out[:,0], yerr=out[:,1], fmt='ko', zorder=0, alpha=0.2)
#    plt.plot(x,y, label=r'$ax^{2} + bx + c$')
#    plt.plot(x,yerr1,'r--')
#    plt.plot(x,yerr2,'r--')
#    plt.ylim(ymin=-0.05,ymax=1.5)
#    plt.xlim(xmin=-29,xmax=-2)
#    plt.xlabel('K Band Absolute Magnitude')
#    plt.ylabel(r'$\sigma$')
#    plt.title(r'$ax^{2} + bx + c$')
#    plt.gca().invert_xaxis()
#    cbar = plt.colorbar()
#    plt.legend()
#    cbar.set_label('z')
#    plt.tight_layout()
    
    
    popt, pcov = curve_fit(func_exp, L, out[:,0], sigma=out[:,1])
#    popt, pcov = curve_fit(func_pol, L, out[:,0], sigma=out[:,1])
    x = np.linspace(-28,-3)
    y = func_exp(x, popt[0], popt[1])
#    y = func_pol(x, popt[0], popt[1])
    
    perr = np.sqrt(np.diag(pcov))
    yerr1 = func_exp(x, popt[0]+perr[0], popt[1]-perr[1])
    yerr2 = func_exp(x, popt[0]-perr[0], popt[1]+perr[1])
#    yerr1 = func_pol(x, popt[0]+perr[0], popt[1]-perr[1])
#    yerr2 = func_pol(x, popt[0]-perr[0], popt[1]+perr[1])
    
#    plt.figure()
#    plt.plot(L, out[:,0], 'o', zorder=1)
#    plt.errorbar(L, out[:,0], yerr=out[:,1], fmt='ko', zorder=0, alpha=0.2)
    plt.plot(x,y, label=labels[n])
#    plt.plot(x,yerr1,'r--')
#    plt.plot(x,yerr2,'r--')
    plt.ylim(ymin=-0.05,ymax=1.5)
    plt.xlim(xmin=-29,xmax=-2)
    plt.xlabel('K Band Absolute Magnitude')
    plt.ylabel(r'$\sigma$')
#    plt.title(labels[n])
    plt.gca().invert_xaxis()
#    cbar = plt.colorbar()
    plt.legend()
#    cbar.set_label('z')

plt.tight_layout(h_pad=0.4, w_pad=0.5)
end = time.time()
print(end-start)

