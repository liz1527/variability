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
Kvarys = fits.open('variable_tables/no06_variables_chi30_2arcsec_DR11data_restframe.fits')[1].data
Jvarys = fits.open('variable_tables/no06_variables_chi30_2arcsec_DR11data_J.fits')[1].data
Hvarys = fits.open('variable_tables/no06_variables_chi30_2arcsec_DR11data_H.fits')[1].data
Ksigtb = Table.read('sigma_tables/quad_epoch_sigma_table_extra_clean_no06_2arcsec.fits')
Jsigtb = Table.read('sigma_tables/quad_epoch_sigma_table_extra_clean_2arcsec_J.fits')

#Kvarys = Kvarys[Kvarys['NUMBER_1']==62253]
#varys = vari_funcs.chandra_only(varys)

### Get K band data ###
Kflux = vari_funcs.flux4_stacks(Kvarys)
#Kflux, Kvarys = vari_funcs.noneg(Kflux, Kvarys)
Kflux, Kfluxerr, newKvarys = vari_funcs.create_quad_error_array(Ksigtb, Kvarys, aper=4)

Kmag = 30 - 2.5*np.log10(Kflux)
Kmagerr = 1.086/(Kflux/Kfluxerr) # taken from http://faculty.virginia.edu/skrutskie/astr3130.s16/notes/astr3130_lec12.pdf?fbclid=IwAR0fe6lNYH8Azj1iVqusb5l-z3xeECx7JBv23ACDV0Xjdq04FHJPD3nPlxE

#Kflux,Kfluxerr = vari_funcs.normalise_flux_and_errors(Kflux, Kfluxerr)

### Get J band data ###
#Jflux = vari_funcs.jflux4_stacks(Jvarys)
#Jflux, Jvarys = vari_funcs.noneg(Jflux, Jvarys)
Jflux, Jfluxerr, newJvarys = vari_funcs.create_quad_error_array_J(Jsigtb, Jvarys, aper=4)
#Jfluxerr = vari_funcs.jfluxerr4_stacks(Jvarys)

Jmag = 30 - 2.5*np.log10(Jflux)
Jmagerr = 1.086/(Jflux/Jfluxerr) # taken from http://faculty.virginia.edu/skrutskie/astr3130.s16/notes/astr3130_lec12.pdf?fbclid=IwAR0fe6lNYH8Azj1iVqusb5l-z3xeECx7JBv23ACDV0Xjdq04FHJPD3nPlxE

#Jflux,Jfluxerr = vari_funcs.normalise_flux_and_errors(Jflux, Jfluxerr)

### Get H band data ###
Hflux = vari_funcs.hflux4_stacks(Hvarys)
Hflux, Hvarys = vari_funcs.noneg(Hflux, Hvarys)
#Hflux, Hfluxerr, newHvarys = vari_funcs.create_quad_error_array(Hsigtb, Hvarys, aper=4)
Hfluxerr = vari_funcs.hfluxerr4_stacks(Hvarys)

mag = 30 - 2.5*np.log10(Hflux)
magerr = 1.086/(Hflux/Hfluxerr) # taken from http://faculty.virginia.edu/skrutskie/astr3130.s16/notes/astr3130_lec12.pdf?fbclid=IwAR0fe6lNYH8Azj1iVqusb5l-z3xeECx7JBv23ACDV0Xjdq04FHJPD3nPlxE

#flux,fluxerr = vari_funcs.normalise_flux_and_errors(flux, fluxerr)

#set up time variable for plot
t = np.linspace(1, 8, num=8)
years = ['05B', '06B', '07B', '08B', '09B', '10B', '11B', '12B']
Kx = [1,3,4,5,6,7,8]
Jx = [1,2,3,4,5,6,7,8]
Hx = [2,3,4,5,6,7,8]
#nums = range(len(Jvarys))


for n in range(len(newKvarys)): #Will plot all with blank plots of those with no J/H
    ### Find corresponding objects in J and H ###
    obnum = newKvarys['NUMBER_1'][n]
    Jmask = np.isin(newJvarys['NUMBER_1'], obnum)
    Hmask = np.isin(Hvarys['NUMBER_1'], obnum)
    if ~np.any(Jmask) or ~np.any(Hmask):
        continue
#    Jn = 
#    Hmask = Hvarys[np.isin(obnum, Hvarys['NUMBER_1'])]
    
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
    
    plt.figure(figsize=[7,8])
    
        
    ### Plot J band ###
    plt.subplot(311)
    if newKvarys['X-ray'][n] == True:
        plt.title(str(obnum)+': X-ray detected')
    else:
        plt.title(str(obnum)+': Not X-ray detected')
#    if newKvarys['X-ray'][n] == True:
#        plt.errorbar(Jx, Jflux[Jmask,:].reshape(8), yerr=Jfluxerr[Jmask,:].reshape(8),fmt='o', color='r')
#    else:
#        plt.errorbar(Jx, Jflux[Jmask,:].reshape(8), yerr=Jfluxerr[Jmask,:].reshape(8),fmt='o', color='b')
    plt.errorbar(Jx, Jflux[Jmask,:].reshape(8), yerr=Jfluxerr[Jmask,:].reshape(8),fmt='o', color='b')
    plt.xlabel('Semester')
    plt.ylabel('J-band flux')
    plt.xlim(xmin=0.5,xmax=8.5)
#    plt.title('Lightcurve of Object '+str(newvarys['NUMBER_05B'][n])+' '+r' $\chi^{2} = $'+str(round(chisq[n], 2)))
    plt.xticks(t, years)
    
    ### Plot H band ###
    plt.subplot(312)
#    if newKvarys['X-ray'][n] == True:
#        plt.errorbar(Hx, Hflux[Hmask,:].reshape(7), yerr=Hfluxerr[Hmask,:].reshape(7),fmt='o', color='r')
#    else:
#        plt.errorbar(Hx, Hflux[Hmask,:].reshape(7), yerr=Hfluxerr[Hmask,:].reshape(7),fmt='o', color='b')
    plt.errorbar(Hx, Hflux[Hmask,:].reshape(7), yerr=Hfluxerr[Hmask,:].reshape(7),fmt='o', color='k')
    plt.xlabel('Semester')
    plt.ylabel('H-band flux')
    plt.xlim(xmin=0.5,xmax=8.5)
#    plt.title('Lightcurve of Object '+str(newvarys['NUMBER_05B'][n])+' '+r' $\chi^{2} = $'+str(round(chisq[n], 2)))
    plt.xticks(t, years)
    
    ### Plot K band ###
    plt.subplot(313)
#    if newKvarys['X-ray'][n] == True:
#        plt.errorbar(Kx, Kflux[n,:], yerr=Kfluxerr[n,:],fmt='o', color='r')
#    else:
#        plt.errorbar(Kx, Kflux[n,:], yerr=Kfluxerr[n,:],fmt='o', color='b')
    plt.errorbar(Kx, Kflux[n,:], yerr=Kfluxerr[n,:],fmt='o', color='r')
    plt.xlabel('Semester')
    plt.ylabel('K-band flux')
    plt.xlim(xmin=0.5,xmax=8.5)
#    plt.title('Lightcurve of Object '+str(newvarys['NUMBER_05B'][n])+' '+r' $\chi^{2} = $'+str(round(chisq[n], 2)))
    plt.xticks(t, years)
    
    
    plt.tight_layout()
    plt.savefig('plots/Chi30Lightcurves/2arcsec/JHK_lightcurves/flux_'+str(n))#+str(varys['NUMBER_05B'][n])+'_lightcurve.png')
    plt.close('all')
