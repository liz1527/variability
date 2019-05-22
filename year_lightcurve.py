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
#varys = fits.open('variable_tables/no06_variables_chi30_DR11data_restframe.fits')[1].data
tbdata = fits.open('mag_flux_tables/mag_flux_table_best_extra_clean_no06.fits')[1].data
sigtb = Table.read('sigma_tables/quad_epoch_sigma_table_extra_clean_no06.fits')

#varys = vari_funcs.chandra_only(varys)

flux = vari_funcs.flux5_stacks(tbdata)
flux, tbdata = vari_funcs.noneg(flux, tbdata)
flux, fluxerr, newtbdata = vari_funcs.create_quad_error_array(sigtb, tbdata)
#flux,fluxerr = vari_funcs.normalise_flux_and_errors(flux, fluxerr)

mag = 30 - 2.5*np.log10(flux)
magerr = 1.086/(flux/fluxerr) # taken from http://faculty.virginia.edu/skrutskie/astr3130.s16/notes/astr3130_lec12.pdf?fbclid=IwAR0fe6lNYH8Azj1iVqusb5l-z3xeECx7JBv23ACDV0Xjdq04FHJPD3nPlxE

#set up time variable for plot
t = np.linspace(1, 8, num=8)
years = ['05B', '06B', '07B', '08B', '09B', '10B', '11B', '12B']
x = [1,3,4,5,6,7,8]

#### reduce to just single lightcurve ###
#obnum = 94904
#mask = newvarys['NUMBER_05B'] == obnum
#flux = flux[mask].reshape(len(x))
#fluxerr = fluxerr[mask].reshape(len(x))
#
#### Plot ###
#plt.figure()
#if newvarys['X-ray'][mask] == 84:
#    plt.errorbar(x, flux, yerr=fluxerr, fmt='o', color='r')
#else:
#    plt.errorbar(x, flux, yerr=fluxerr, fmt='o', color='b')
#plt.xlabel('Semester')
#plt.ylabel('K-band flux')
#plt.title('Lightcurve of Object '+str(obnum))
#plt.xticks(t, years)
#plt.tight_layout()
##plt.savefig('Chi40Lightcurves/cleaned/no06/mag_'+str(n))#+str(varys['NUMBER_05B'][n])+'_lightcurve.png')
##plt.close('all')

nums = [94904,282832,250101,206455, 128670,175673]
#plt.figure()
### Plot multiple curves ###
for n, obnum in enumerate(nums):
    ### reduce to just single lightcurve ###
    mask = newtbdata['NUMBER_05B'] == obnum
    obflux = flux[mask].reshape(len(x))
    obfluxerr = fluxerr[mask].reshape(len(x))
    obmag = mag[mask].reshape(len(x))
    obmagerr = magerr[mask].reshape(len(x))
    
    ### Plot ###
#    plt.figure()
#    if newtbdata['X-ray'][mask] == 84:
#        plt.errorbar(x, obflux, yerr=obfluxerr, fmt='o', color='r')
#    else:
#        plt.errorbar(x, obflux, yerr=obfluxerr, fmt='o', color='b')
#    plt.xlabel('Semester')
#    plt.ylabel('K-band flux')
#    plt.title('Lightcurve of Object '+str(obnum))
#    plt.xticks(t, years)
#    plt.tight_layout()
#    #plt.savefig('Chi40Lightcurves/cleaned/no06/mag_'+str(n))#+str(varys['NUMBER_05B'][n])+'_lightcurve.png')
#    #plt.close('all')
    
    plt.figure()
#    plt.subplot(3,2,n+1)
    if newtbdata['X-ray'][mask] == 84:
        plt.errorbar(x, obmag, yerr=obmagerr, fmt='o', color='r')
    else:
        plt.errorbar(x, obmag, yerr=obmagerr, fmt='o', color='b')
    axes = plt.gca()
    ylims = axes.get_ylim()
    ymid = (ylims[1]+ylims[0])/2
    plt.ylim(ymin=ymid-0.2, ymax=ymid+0.2)
    plt.xlabel('Semester')
    plt.ylabel('3" K-band Magnitude')
#    plt.title('Lightcurve of Object '+str(obnum))
    plt.xticks(t, years)
    plt.tight_layout()
    plt.savefig(str(obnum)+'poster_lightcurve.pdf')
