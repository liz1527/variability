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




tbdata = fits.open('variable_tables/J_and_K_variables_varystats_DR11data_monthdata.fits')[1].data
month_sigtb = Table.read('sigma_tables/month_quad_epoch_sigma_table_K_extra_quad_clean_2arcsec_noneg.fits')
sem_sigtb = Table.read('sigma_tables/quad_epoch_sigma_table_K_extra_clean_2arcsec_noneg.fits')

months = ['sep05','oct05','nov05','dec05', 'jan06', 'dec06', 'jan07',  
          'aug07', 'sep07', 'oct07', 'sep08', 'oct08', 'nov08', 'jul09',  
          'aug09', 'sep09', 'oct09', 'nov09', 'dec09', 'jan10', 'feb10', 
          'aug10', 'sep10', 'oct10', 'nov10', 'dec10', 'jan11', 'feb11', 
          'aug11', 'sep11', 'oct11', 'nov11', 'dec11', 'jan12', 'feb12', 
          'jul12', 'aug12', 'sep12', 'oct12', 'nov12']
semesters = ['05B', '07B', '08B', '09B', '10B', '11B', '12B']

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

### set up sem tick details ###
t_year = np.linspace(1, 8, num=8)
years = ['05B', '06B', '07B', '08B', '09B', '10B', '11B', '12B']
x_sem = [1,3,4,5,6,7,8]

### get month data ###
month_flux, month_fluxerr, tbdata =  vari_funcs.k_mag_flux.create_quad_error_array_month(month_sigtb, tbdata)

month_mag = 30 - 2.5*np.log10(month_flux)
month_mag += 1.9 #to get to AB mag
month_magerr = 1.086/(month_flux/month_fluxerr) # taken from http://faculty.virginia.edu/skrutskie/astr3130.s16/notes/astr3130_lec12.pdf?fbclid=IwAR0fe6lNYH8Azj1iVqusb5l-z3xeECx7JBv23ACDV0Xjdq04FHJPD3nPlxE

month_mag[np.isinf(month_mag)] = np.nan# reset infs to nans
month_magerr[np.isinf(month_magerr)] = np.nan# reset infs to nans

month_chisq = vari_funcs.vary_stats.my_chisquare_err(month_flux, month_fluxerr)

### get semester data ###
sem_flux = tbdata['Flux_K']
sem_fluxerr = tbdata['Fluxerr_K']

sem_mag = 30 - 2.5*np.log10(sem_flux)
sem_mag += 1.9 #to get to AB mag
sem_magerr = 1.086/(sem_flux/sem_fluxerr) # taken from http://faculty.virginia.edu/skrutskie/astr3130.s16/notes/astr3130_lec12.pdf?fbclid=IwAR0fe6lNYH8Azj1iVqusb5l-z3xeECx7JBv23ACDV0Xjdq04FHJPD3nPlxE

sem_chisq = tbdata['Chi_K']

for n in range(len(tbdata)):
    ### find maxes and mins to define lims ###
    sem_max = np.nanmax(sem_mag[n,:])
    month_max = np.nanmax(month_mag[n,:])
    max = np.nanmax([sem_max, month_max])
    sem_min = np.nanmin(sem_mag[n,:])
    month_min = np.nanmin(month_mag[n,:])
    min = np.nanmin([sem_min, month_min])
    
    
    plt.figure(figsize=[10,10])
    ### plot year curve ###
    plt.subplot(211)
    if tbdata['X-ray_1'][n] == 84: #bool removed for some reason
        plt.errorbar(x_sem, sem_mag[n,:], yerr=sem_magerr[n,:],fmt='ro')
    else:
        plt.errorbar(x_sem, sem_mag[n,:], yerr=sem_magerr[n,:],fmt='bo')
    plt.ylim(ymin = min-0.025, ymax=max+0.025)
    plt.gca().invert_yaxis()
    plt.xlabel('Semester')
    plt.ylabel('K-band Magnitude')
    plt.title('Lightcurve of Object '+str(tbdata['ID'][n])+' '+r' $\chi^{2} = $'+str(round(sem_chisq[n], 2)))
    plt.xticks(t_year, years)
    
    plt.subplot(212)
    if tbdata['X-ray_1'][n] == 84:
        plt.errorbar(x_months, month_mag[n,:], yerr=month_magerr[n,:], fmt='ro')
    else:
        plt.errorbar(x_months, month_mag[n,:], yerr=month_magerr[n,:], fmt='bo')
    plt.xticks(tick_inds, month_ticks, rotation = 'vertical')
    plt.ylim(ymin=min-0.025, ymax=max+0.025)
    plt.gca().invert_yaxis()
    plt.ylabel('K-band Magnitude')
    plt.xlabel('Month')
    plt.title('Lightcurve of Object '+str(tbdata['ID'][n])+' '+r' $\chi^{2} = $'+str(round(month_chisq[n], 2)))
#    plt.legend()
    plt.tight_layout()
    plt.savefig('plots/new_catalogue/month_and_year_lightcurves/mag_'+str(n), 
                overwrite=True)#+str(varys['NUMBER_05B'][n])+'_lightcurve.png')
    plt.close('all')


#plt.ylim(ymin=40000)
#axes = plt.gca()
#ylims = axes.get_ylim()
#ymid = (ylims[1]+ylims[0])/2
#plt.ylim(ymin=ymid-0.25, ymax=ymid+0.25)
#plt.title('Light curve for object number %i' % obnum)
#plt.ylabel('K-band magnitude')
#plt.title('Lightcurve of Object '+str(obnum)+' '+r' $\chi^{2} = $'+str(round(chisq[0], 2)))
