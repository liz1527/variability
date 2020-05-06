#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 13:53:26 2019

code to create month lightcurves

@author: ppxee
"""

### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
from astropy.table import Table #for handling tables
import numpy as np #for handling arrays
#import math
#from astropy.stats import median_absolute_deviation
import vari_funcs #my module to help run code neatly
#from astropy.cosmology import FlatLambdaCDM
#from astropy import units as u
plt.close('all') #close any open plots


#months = ['sep05','oct05','nov05','dec05', 'jan06', #'dec06', 
#          'jan07', 'aug07', 'sep07', 'oct07', 'sep08', 'oct08', 'nov08', 
#          'jul09', 'aug09', 'sep09', 'oct09', 'nov09', 'dec09', 'jan10', 
#          'feb10', 'aug10', 'sep10', 'oct10', 'nov10', 'dec10', 'jan11', #'feb11', 
#          'aug11', 'sep11', 'oct11', 'nov11', 'dec11', 'jan12', 'feb12', 
#          'jul12', 'aug12', 'sep12', 'oct12', 'nov12']


months = ['sep05','oct05','nov05','dec05', 'jan06', 'dec06', 'jan07',  
          'aug07', 'sep07', 'oct07', 'sep08', 'oct08', 'nov08', 'jul09',  
          'aug09', 'sep09', 'oct09', 'nov09', 'dec09', 'jan10', 'feb10', 
          'aug10', 'sep10', 'oct10', 'nov10', 'dec10', 'jan11', 'feb11', 
          'aug11', 'sep11', 'oct11', 'nov11', 'dec11', 'jan12', 'feb12', 
          'jul12', 'aug12', 'sep12', 'oct12', 'nov12']

### set up month tick details ###
month_info = fits.open('Images/Convolving_Images/monthly_numbers.fits')[1].data #get month count data
full_months = month_info['Month'] #extract month nanes
tick_inds = np.load('Images/Convolving_Images/tick_inds_K.npy') #load tick locations
mask = np.zeros(len(full_months)) #set up mask
mask[tick_inds] = 1
mask = mask.astype(bool)
month_ticks = np.copy(full_months)
month_ticks = month_ticks[mask]#retrieve tick details

x = np.arange(0, len(month_info['Frames in v11']))
mask = np.isin(full_months, months)
x_months = x[mask]


tbdata = fits.open('mag_flux_tables/K/month/month_mag_flux_table_best_K_extra_quad_clean.fits')[1].data
sigtb = Table.read('sigma_tables/month_quad_epoch_sigma_table_K_extra_quad_clean_2arcsec_noneg_pvalue.fits')

print('Removing edges')
tbdata = vari_funcs.field_funcs.remove_edges(tbdata, 'sep05')
print('Removed edges')

### get flux data ###
print('getting flux')
#flux = vari_funcs.k_mag_flux.month_flux_stacks(tbdata)
flux, fluxerr, tbdata =  vari_funcs.k_mag_flux.create_quad_error_array_month(sigtb, tbdata)
print('got flux')
#flux[flux <= 0] = np.nan
flux, fluxerr, tbdata = vari_funcs.flux_funcs.nanneg(flux, fluxerr, tbdata)

del tbdata


#%% Plot ###
### average flux data ###
#fluxn = vari_funcs.flux_funcs.normalise_flux(flux)
fluxn, fluxerrn = vari_funcs.flux_funcs.normalise_flux_and_errors(flux, fluxerr)
avgflux = np.nanmean(fluxn, axis=0)
avgfluxn = np.nanmedian(fluxn, axis=0)

meanerr = np.sqrt(np.sum(np.square(fluxerrn), axis=0))/len(fluxerr)

plt.figure(figsize=[10,5])
plt.plot(x_months, avgfluxn, 'ro')
#plt.errorbar(x_months, flux[0,:], yerr=fluxerr[0,:], fmt='ro')
plt.xticks(tick_inds, month_ticks, rotation = 'vertical')
plt.ylabel('Flux')
plt.xlabel('Month')
plt.title('Median Light Curve')
plt.tight_layout()

plt.figure(figsize=[10,5])
#plt.plot(x_months, avgflux, 'ro')
plt.errorbar(x_months, avgflux, yerr=meanerr, fmt='ro')
plt.xticks(tick_inds, month_ticks, rotation = 'vertical')
plt.ylabel('Flux')
plt.xlabel('Month')
plt.title('Mean Light Curve')
plt.tight_layout()

#plt.ylim(ymin=40000)
#axes = plt.gca()
#ylims = axes.get_ylim()
#ymid = (ylims[1]+ylims[0])/2
#plt.ylim(ymin=ymid-0.25, ymax=ymid+0.25)
#plt.title('Light curve for object number %i' % obnum)
#plt.ylabel('K-band magnitude')
#plt.title('Lightcurve of Object '+str(obnum)+' '+r' $\chi^{2} = $'+str(round(chisq[0], 2)))
