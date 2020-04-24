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




tbdata = fits.open('mag_flux_tables/K/month/month_mag_flux_table_best_K_extra_quad_clean.fits')[1].data
sigtb = Table.read('sigma_tables/month_quad_epoch_sigma_table_K_extra_quad_clean_2arcsec_noneg.fits')
obnum = 173459

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

#mask = fitsdata['NUMBER_1'] == ob
obdata = tbdata[tbdata['NUMBER'] == obnum]

flux = vari_funcs.k_mag_flux.month_flux_stacks(obdata)
flux, fluxerr, obdata =  vari_funcs.k_mag_flux.create_quad_error_array_month(sigtb, obdata)

mag = 30 - 2.5*np.log10(flux)
mag += 1.9 #to get to AB mag
magerr = 1.086/(flux/fluxerr) # taken from http://faculty.virginia.edu/skrutskie/astr3130.s16/notes/astr3130_lec12.pdf?fbclid=IwAR0fe6lNYH8Azj1iVqusb5l-z3xeECx7JBv23ACDV0Xjdq04FHJPD3nPlxE

chisq = vari_funcs.vary_stats.my_chisquare_err(flux, fluxerr)

plt.figure(figsize=[10,5])
#plt.errorbar(x_months, flux[0,:], yerr=fluxerr[0,:], fmt='ro')
plt.errorbar(x_months, mag[0,:], yerr=magerr[0,:], fmt='ro')
plt.xticks(tick_inds, month_ticks, rotation = 'vertical')
plt.ylabel('Flux')
plt.xlabel('Month')
plt.title('Light curve for '+str(obnum)+r' $\chi^{2}$ = ' + str(chisq[0]))
plt.legend()
plt.tight_layout()


#plt.ylim(ymin=40000)
#axes = plt.gca()
#ylims = axes.get_ylim()
#ymid = (ylims[1]+ylims[0])/2
#plt.ylim(ymin=ymid-0.25, ymax=ymid+0.25)
#plt.title('Light curve for object number %i' % obnum)
plt.ylabel('K-band magnitude')
#plt.title('Lightcurve of Object '+str(obnum)+' '+r' $\chi^{2} = $'+str(round(chisq[0], 2)))
