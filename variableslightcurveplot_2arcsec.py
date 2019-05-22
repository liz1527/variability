#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  8 16:44:00 2018

@author: ppxee
"""

### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
from astropy.table import Table #for handling tables
import numpy as np #for handling arrays
from astropy.stats import median_absolute_deviation
import vari_funcs #my module to help run code neatly
plt.close('all') #close any open plots

font = {'family' : 'DejaVu Sans',
        'weight' : 'normal',
        'size'   : 14}

plt.rc('font', **font)

### Open the fits files and get data ###
varys = fits.open('variable_tables/no06_variables_chi30_2arcsec_deviant.fits')[1].data
sigtb = Table.read('sigma_tables/quad_epoch_sigma_table_extra_clean_no06_2arcsec.fits')

#varys = vari_funcs.chandra_only(varys)

flux = vari_funcs.flux4_stacks(varys)
flux, varys = vari_funcs.noneg(flux, varys)
flux, fluxerr, newvarys = vari_funcs.create_quad_error_array(sigtb, varys, aper=4)

mag = 30 - 2.5*np.log10(flux)
magerr = 1.086/(flux/fluxerr) # taken from http://faculty.virginia.edu/skrutskie/astr3130.s16/notes/astr3130_lec12.pdf?fbclid=IwAR0fe6lNYH8Azj1iVqusb5l-z3xeECx7JBv23ACDV0Xjdq04FHJPD3nPlxE

#flux,fluxerr = vari_funcs.normalise_flux_and_errors(flux, fluxerr)

#set up time variable for plot
t = np.linspace(1, 8, num=8)
years = ['05B', '06B', '07B', '08B', '09B', '10B', '11B', '12B']
x = [1,3,4,5,6,7,8]

chisq = vari_funcs.my_chisquare_err(flux, fluxerr)
mad = median_absolute_deviation(flux, axis=1)


mask = np.zeros(np.shape(mad))
for n in range(len(newvarys)):
#    plt.figure()
#    if newvarys['X-ray'][n] == True:
#        plt.errorbar(x, mag[n,:], yerr=magerr[n,:],fmt='o', color='r')
#    else:
#        plt.errorbar(x, mag[n,:], yerr=magerr[n,:],fmt='o', color='b')
#    plt.xlabel('Semester')
#    plt.ylabel('K-band magnitude')
#    plt.title('Lightcurve of Object '+str(newvarys['NUMBER_05B'][n])+' '+r' $\chi^{2} = $'+str(round(chisq[n], 2)))
#    plt.xticks(t, years)
#    plt.tight_layout()
#    plt.savefig('plots/Chi30Lightcurves/2arcsec/mag_'+str(n))#+str(varys['NUMBER_05B'][n])+'_lightcurve.png')
#    plt.close('all')
    
    plt.figure()
    if newvarys['X-ray'][n] == True:
        plt.errorbar(x, flux[n,:], yerr=fluxerr[n,:],fmt='o', color='r')
    else:
        plt.errorbar(x, flux[n,:], yerr=fluxerr[n,:],fmt='o', color='b')
    plt.xlabel('Semester')
    plt.ylabel('K-band flux')
    plt.title('Lightcurve of Object '+str(newvarys['NUMBER_05B'][n])+' '+r' $\chi^{2} = $'+str(round(chisq[n], 2)))
    plt.xticks(t, years)
    plt.tight_layout()
    plt.savefig('plots/Chi30Lightcurves/2arcsec/flux/deviant/flux_'+str(n))#+str(varys['NUMBER_05B'][n])+'_lightcurve.png')
    plt.close('all')
