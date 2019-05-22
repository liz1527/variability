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
#plt.close('all') #close any open plots

font = {'family' : 'DejaVu Sans',
        'weight' : 'normal',
        'size'   : 14}

plt.rc('font', **font)

### Open the fits files and get data ###
tbdata = fits.open('mag_flux_tables/mag_flux_table_best_extra_clean_no06_DR11data.fits')[1].data
sigtb = Table.read('sigma_tables/quad_epoch_sigma_table_cleaned_no06.fits')

obnum = 173520 #252446

tbdata = tbdata[tbdata['NUMBER_1']==obnum]

flux = vari_funcs.flux5_stacks(tbdata)
flux, tbdata = vari_funcs.noneg(flux, tbdata)
#fluxerr = vari_funcs.fluxerr5_stacks(tbdata)
flux, fluxerr, tbdata = vari_funcs.create_quad_error_array(sigtb, tbdata)

mag = 30 - 2.5*np.log10(flux)
magerr = 1.086/(flux/fluxerr) # taken from http://faculty.virginia.edu/skrutskie/astr3130.s16/notes/astr3130_lec12.pdf?fbclid=IwAR0fe6lNYH8Azj1iVqusb5l-z3xeECx7JBv23ACDV0Xjdq04FHJPD3nPlxE

#flux,fluxerr = vari_funcs.normalise_flux_and_errors(flux, fluxerr)

#set up time variable for plot
t = np.linspace(1, 8, num=8)
years = ['05B', '06B', '07B', '08B', '09B', '10B', '11B', '12B']
x = [1,3,4,5,6,7,8]

chisq = vari_funcs.my_chisquare_err(flux, fluxerr)

plt.figure()
if tbdata['X-ray'] == True:
    plt.errorbar(x, mag[0,:], yerr=magerr[0,:],fmt='o', color='r')
else:
    plt.errorbar(x, mag[0,:], yerr=magerr[0,:],fmt='o', color='b')
#plt.errorbar(x, mag[0,:], yerr=magerr[0,:],fmt='o', color='b')
plt.xlabel('Semester')
plt.ylabel('K-band magnitude')
#plt.title('Lightcurve of Object '+str(obnum)+' '+r' $\chi^{2} = $'+str(round(chisq[0], 2)))
plt.xticks(t, years)
#plt.ylim(ymin=19.03,ymax=19.69)
plt.tight_layout()
#plt.savefig('plots/Chi30Lightcurves/extra_clean_no06/mag_'+str(obnum)+'_notitle.pdf')
#    plt.savefig('Chi40Lightcurves/cleaned/no06/mag_'+str(n))#+str(varys['NUMBER_05B'][n])+'_lightcurve.png')
#plt.close('all')
