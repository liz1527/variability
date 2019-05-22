#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  8 16:44:00 2018

Code to look at problem with 08B lightcurves

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
tbdata = fits.open('mag_flux_tables/mag_flux_table_best_1519match.fits')[1].data
varys = fits.open('variable_tables/variables_chi50.fits')[1].data
sigtb = Table.read('quad_epoch_sigma_table.fits')

#varys = vari_funcs.chandra_only(varys)

flux = vari_funcs.flux5_stacks(varys)
flux, varys = vari_funcs.noneg(flux, varys)
flux, fluxerr, newvarys = vari_funcs.create_quad_error_array(sigtb, varys)

flux,fluxerr = vari_funcs.normalise_flux_and_errors(flux, fluxerr)

#set up time variable for plot
t = np.linspace(1, 8, num=8)
years = ['05B', '06B', '07B', '08B', '09B', '10B', '11B', '12B']

chisq = vari_funcs.my_chisquare_err(flux, fluxerr)
mad = median_absolute_deviation(flux, axis=1)

mask = np.zeros(np.shape(mad))
for n in range(len(newvarys)):
    ### Need to find those that go wibbly in 08B ###
    diff08 = 1-flux[n,3]
    if diff08 > 3*mad[n] and diff08 > 0.1:
        mask[n] = 1
#        print(diff08)
#        print(mad[n])
#        plt.figure()
#        if varys['X-ray'][n] == True:
#            plt.errorbar(t, flux[n,:], yerr=fluxerr[n,:],fmt='o', color='r')
#        else:
#            plt.errorbar(t, flux[n,:], yerr=fluxerr[n,:],fmt='o', color='b')
#        plt.xlabel('Semester')
#        plt.ylabel('Flux')
#        plt.title('Lightcurve of Object '+str(newvarys['NUMBER_05B'][n])+' '+r' $\chi^{2} = $'+str(round(chisq[n], 2)))
#        plt.xticks(t, years)
#        plt.tight_layout()
#    #    axes = plt.gca()
#    #    ylims = axes.get_ylim()
#    #    ymid = (ylims[1]+ylims[0])/2
#    #    plt.ylim(ymin=ymid-0.26, ymax=ymid+0.26)
#        plt.savefig('Chi50Lightcurves/bad08B/'+str(n))#+str(varys['NUMBER_05B'][n])+'_lightcurve.png')
#        plt.close('all')

mask = mask.astype(bool)
bad = newvarys[mask]

badtb = Table(bad)
badtb.write('variable_tables/bad08_variables.fits')

### plot positions ###
fullra = tbdata['ALPHA_J2000_05B']
fulldec = tbdata['DELTA_J2000_05B']
varyra = bad['ALPHA_J2000_05B']
varydec = bad['DELTA_J2000_05B']

plt.scatter(fullra, fulldec, marker='.', c='b', s=1)
plt.scatter(varyra, varydec, marker='*', c='m')