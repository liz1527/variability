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
varys_Kex = fits.open('variable_tables/J/K_extraction/J_variables_chi32_noneg_DR11data.fits')[1].data
sigtb_Kex = Table.read('sigma_tables/quad_epoch_sigma_table_J_extra_clean_2arcsec_noneg.fits')
varys_Jex = fits.open('mag_flux_tables/J/J_extraction/mag_flux_table_J_extra_clean_DR11data.fits')[1].data
sigtb_Jex = Table.read('sigma_tables/quad_epoch_sigma_table_extra_clean_2arcsec_J_noneg_J_extracted.fits')

mask = np.isin(varys_Jex['ID'], varys_Kex['ID'])
varys_Jex = varys_Jex[mask]
#varys = vari_funcs.chandra_only(varys)

flux_Kex = vari_funcs.j_mag_flux.flux4_stacks(varys_Kex)
flux_Kex, varys_Kex = vari_funcs.flux_funcs.noneg(flux_Kex, varys_Kex)
flux_Kex, fluxerr_Kex, newvarys_Kex = vari_funcs.j_mag_flux.create_quad_error_array_J(sigtb_Kex, varys_Kex, aper=4)

flux_Jex = vari_funcs.j_mag_flux.flux4_stacks(varys_Jex)
flux_Jex, varys_Jex = vari_funcs.flux_funcs.noneg(flux_Jex, varys_Jex)
flux_Jex, fluxerr_Jex, newvarys_Jex = vari_funcs.j_mag_flux.create_quad_error_array_J(sigtb_Jex, varys_Jex, aper=4)

mag_Kex = 30 - 2.5*np.log10(flux_Kex)
mag_Kex += 1.9 #to get to AB mag
magerr_Kex = 1.086/(flux_Kex/fluxerr_Kex) # taken from http://faculty.virginia.edu/skrutskie/astr3130.s16/notes/astr3130_lec12.pdf?fbclid=IwAR0fe6lNYH8Azj1iVqusb5l-z3xeECx7JBv23ACDV0Xjdq04FHJPD3nPlxE

mag_Jex = 30 - 2.5*np.log10(flux_Jex)
mag_Jex += 1.9 #to get to AB mag
magerr_Jex = 1.086/(flux_Jex/fluxerr_Jex) # taken from http://faculty.virginia.edu/skrutskie/astr3130.s16/notes/astr3130_lec12.pdf?fbclid=IwAR0fe6lNYH8Azj1iVqusb5l-z3xeECx7JBv23ACDV0Xjdq04FHJPD3nPlxE

#flux,fluxerr = vari_funcs.normalise_flux_and_errors(flux, fluxerr)

#set up time variable for plot
t = np.linspace(1, 8, num=8)
years = ['05B', '06B', '07B', '08B', '09B', '10B', '11B', '12B']
x = [1,2,3,4,5,6,7,8]

chisq_Kex = vari_funcs.vary_stats.my_chisquare_err(flux_Kex, fluxerr_Kex)
chisq_Jex = vari_funcs.vary_stats.my_chisquare_err(flux_Jex, fluxerr_Jex)
#mad = median_absolute_deviation(flux, axis=1)


#mask = np.zeros(np.shape(mad))
for n in range(len(newvarys_Kex)-1):#2):
    plt.figure()
#    if newvarys_Kex['X-ray'][n] == True:
#        plt.errorbar(x, mag_Kex[n,:], yerr=magerr_Kex[n,:],fmt='o', color='r')
#        plt.errorbar(x, mag_Jex[n,:], yerr=magerr_Jex[n,:],fmt='o', color='r')
#    else:
#        plt.errorbar(x, mag_Kex[n,:], yerr=magerr_Kex[n,:],fmt='o', color='b')
#        plt.errorbar(x, mag_Jex[n,:], yerr=magerr_Jex[n,:],fmt='o', color='b')
    
    plt.errorbar(x, mag_Kex[n,:], yerr=magerr_Kex[n,:],fmt='o', label='K extraction')
    plt.errorbar(x, mag_Jex[n,:], yerr=magerr_Jex[n,:],fmt='o', label='J extraction')
    
    plt.xlabel('Semester')
    plt.ylabel('J-band magnitude')
    plt.title('Lightcurve of Object '+str(newvarys_Kex['NUMBER_05B'][n])+' '+
              r' $\chi^{2}_{K} = $'+str(round(chisq_Kex[n], 2))+
              r' $\chi^{2}_{J} = $'+str(round(chisq_Jex[n], 2)))
    plt.xticks(t, years)
    plt.legend()
    plt.tight_layout()
    plt.savefig('plots/new_catalogue/Chi30Lightcurves/J_lightcurves/mag_'+str(varys_Kex['ID'][n]))#+'_lightcurve.png')+str(n)
    plt.close('all')
    
#    plt.figure()
#    if newvarys['X-ray'][n] == True:
#        plt.errorbar(x, flux[n,:], yerr=fluxerr[n,:],fmt='o', color='r')
#    else:
#        plt.errorbar(x, flux[n,:], yerr=fluxerr[n,:],fmt='o', color='b')
#    plt.xlabel('Semester')
#    plt.ylabel('K-band flux')
#    plt.title('Lightcurve of Object '+str(newvarys['NUMBER_05B'][n])+' '+r' $\chi^{2} = $'+str(round(chisq[n], 2)))
#    plt.xticks(t, years)
#    plt.tight_layout()
#    plt.savefig('plots/new_catalogue/Chi30Lightcurves/neg_only/not_deviant/flux_'+str(n))#+str(varys['NUMBER_05B'][n])+'_lightcurve.png')
#    plt.close('all')
