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
import vari_funcs #my module to help run code neatly
plt.close('all') #close any open plots

font = {'family' : 'DejaVu Sans',
        'weight' : 'normal',
        'size'   : 14}

plt.rc('font', **font)

### Open the fits files and get data ###
varydata = Table.read('variable_tables/J_and_K_variables_month_varystats_DR11data.fits')

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

### Plot lightcurves
for obdata in varydata:
    obnum = obdata['ID']
    
    ### Get K band data ###
    Kflux = obdata['Month_Flux_K']
    Kfluxerr = obdata['Month_Fluxerr_K']
    
#    Kmag = 30 - 2.5*np.log10(Kflux)
#    Kmagerr = 1.086/(Kflux/Kfluxerr) # taken from http://faculty.virginia.edu/skrutskie/astr3130.s16/notes/astr3130_lec12.pdf?fbclid=IwAR0fe6lNYH8Azj1iVqusb5l-z3xeECx7JBv23ACDV0Xjdq04FHJPD3nPlxE
    
    Kflux, Kfluxerr = vari_funcs.flux_funcs.normalise_single_flux_and_errors(Kflux, Kfluxerr)
    
    ### Get J band data ###
    Jflux = obdata['Month_Flux_J']
    Jfluxerr = obdata['Month_Fluxerr_J']
    
#    Jmag = 30 - 2.5*np.log10(Jflux)
#    Jmagerr = 1.086/(Jflux/Jfluxerr) # taken from http://faculty.virginia.edu/skrutskie/astr3130.s16/notes/astr3130_lec12.pdf?fbclid=IwAR0fe6lNYH8Azj1iVqusb5l-z3xeECx7JBv23ACDV0Xjdq04FHJPD3nPlxE
    Jflux, Jfluxerr = vari_funcs.flux_funcs.normalise_single_flux_and_errors(Jflux, Jfluxerr)

    plt.figure(figsize=[11,9])
    
    if obdata['X-ray'] == True:
        plt.suptitle(str(obnum)+': X-ray detected')
    else:
        plt.suptitle(str(obnum)+': Not X-ray detected')

    ### Plot J band ###
#    plt.subplot(211)
#    if newKvarys['X-ray'][n] == True:
#        plt.errorbar(Jx, Jflux[Jmask,:].reshape(8), yerr=Jfluxerr[Jmask,:].reshape(8),fmt='o', color='r')
#    else:
#        plt.errorbar(Jx, Jflux[Jmask,:].reshape(8), yerr=Jfluxerr[Jmask,:].reshape(8),fmt='o', color='b')
    plt.errorbar(Jx_months, Jflux, yerr=Jfluxerr,fmt='.', color='b', 
                 zorder=0, alpha=0.2)
    plt.plot(Jx_months, Jflux, 'bo', label='J-band')
#    plt.xlabel('Month')
#    plt.ylabel('J-band flux')
##    plt.xlim(xmin=0.5,xmax=8.5)
#    plt.title('$\chi^{2}_{J} = $'+str(round(obdata['Month_Chi_J'], 2)))
#    plt.xticks(tick_inds, month_ticks, rotation = 'vertical')
#    plt.xticks(inds, month_ticks, rotation = 'vertical')
#    plt.grid(True, axis='x')
   
    ### Plot K band ###
#    plt.subplot(212)
#    if newKvarys['X-ray'][n] == True:
#        plt.errorbar(Kx, Kflux[n,:], yerr=Kfluxerr[n,:],fmt='o', color='r')
#    else:
#        plt.errorbar(Kx, Kflux[n,:], yerr=Kfluxerr[n,:],fmt='o', color='b')
    plt.errorbar(Kx_months, Kflux, yerr=Kfluxerr, fmt='.', color='r', 
                 zorder=0, alpha=0.2)
    plt.plot(Kx_months, Kflux, 'ro', label='K-band')
    plt.xlabel('Month')
    plt.ylabel('Normalised Flux')
#    plt.ylabel('K-band flux')
#    plt.xlim(xmin=0.5,xmax=8.5)
#    plt.title('$\chi^{2}_{K} = $'+str(round(obdata['Month_Chi_K'], 2)))
#    plt.title('Lightcurve of Object '+str(newvarys['NUMBER_05B'][n])+' '+r' $\chi^{2} = $'+str(round(chisq[n], 2)))
#    plt.xticks(tick_inds, month_ticks, rotation = 'vertical')
    plt.xticks(inds, month_ticks, rotation = 'vertical')
    plt.legend()
#    plt.grid(True, axis='x')
    plt.tight_layout()
    plt.subplots_adjust(top=0.89)
#    break
    plt.savefig('plots/new_catalogue/Chi30Lightcurves/month_JK_lightcurves/normflux_'+str(obnum))#+str(varys['NUMBER_05B'][n])+'_lightcurve.png')
    plt.close('all')
