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
varydata = Table.read('variable_tables/J_and_K_high_chi_JK_variables.fits')


#Kvarys = Kvarys[Kvarys['NUMBER_1']==62253]
#varys = vari_funcs.chandra_only(varys)

#Jflux,Jfluxerr = vari_funcs.normalise_flux_and_errors(Jflux, Jfluxerr)

#flux,fluxerr = vari_funcs.normalise_flux_and_errors(flux, fluxerr)

#set up time variable for plot
t = np.linspace(1, 8, num=8)
years = ['05B', '06B', '07B', '08B', '09B', '10B', '11B', '12B']
Kx = [1,3,4,5,6,7,8]
Jx = [1,2,3,4,5,6,7,8]
#nums = range(len(Jvarys))


for obdata in varydata:
    obnum = obdata['ID']
    
    ### Get K band data ###
    Kflux = obdata['Flux_K']
    Kfluxerr = obdata['Fluxerr_K']
    
#    Kmag = 30 - 2.5*np.log10(Kflux)
#    Kmagerr = 1.086/(Kflux/Kfluxerr) # taken from http://faculty.virginia.edu/skrutskie/astr3130.s16/notes/astr3130_lec12.pdf?fbclid=IwAR0fe6lNYH8Azj1iVqusb5l-z3xeECx7JBv23ACDV0Xjdq04FHJPD3nPlxE
    
    #Kflux,Kfluxerr = vari_funcs.normalise_flux_and_errors(Kflux, Kfluxerr)
    
    ### Get J band data ###
    Jflux = obdata['Flux_J']
    Jfluxerr = obdata['Fluxerr_J']
    
#    Jmag = 30 - 2.5*np.log10(Jflux)
#    Jmagerr = 1.086/(Jflux/Jfluxerr) # taken from http://faculty.virginia.edu/skrutskie/astr3130.s16/notes/astr3130_lec12.pdf?fbclid=IwAR0fe6lNYH8Azj1iVqusb5l-z3xeECx7JBv23ACDV0Xjdq04FHJPD3nPlxE

    plt.figure(figsize=[7,8])
    
    if obdata['X-ray'] == True:
        plt.suptitle(str(obnum)+': X-ray detected'+' '+r' $\chi^{2}_{JK} = $'+str(round(obdata['Chi_JK'], 2)))
    else:
        plt.suptitle(str(obnum)+': Not X-ray detected'+' '+r' $\chi^{2}_{JK} = $'+str(round(obdata['Chi_JK'], 2)))

    ### Plot J band ###
    plt.subplot(211)
#    if newKvarys['X-ray'][n] == True:
#        plt.errorbar(Jx, Jflux[Jmask,:].reshape(8), yerr=Jfluxerr[Jmask,:].reshape(8),fmt='o', color='r')
#    else:
#        plt.errorbar(Jx, Jflux[Jmask,:].reshape(8), yerr=Jfluxerr[Jmask,:].reshape(8),fmt='o', color='b')
    plt.errorbar(Jx, Jflux, yerr=Jfluxerr,fmt='o', color='b')
    plt.xlabel('Semester')
    plt.ylabel('J-band flux')
    plt.xlim(xmin=0.5,xmax=8.5)
    plt.title('$\chi^{2}_{J} = $'+str(round(obdata['Chi_J'], 2)))
    plt.xticks(t, years)
    
    ### Plot K band ###
    plt.subplot(212)
#    if newKvarys['X-ray'][n] == True:
#        plt.errorbar(Kx, Kflux[n,:], yerr=Kfluxerr[n,:],fmt='o', color='r')
#    else:
#        plt.errorbar(Kx, Kflux[n,:], yerr=Kfluxerr[n,:],fmt='o', color='b')
    plt.errorbar(Kx, Kflux, yerr=Kfluxerr,fmt='o', color='r')
    plt.xlabel('Semester')
    plt.ylabel('K-band flux')
    plt.xlim(xmin=0.5,xmax=8.5)
    plt.title('$\chi^{2}_{K} = $'+str(round(obdata['Chi_K'], 2)))
#    plt.title('Lightcurve of Object '+str(newvarys['NUMBER_05B'][n])+' '+r' $\chi^{2} = $'+str(round(chisq[n], 2)))
    plt.xticks(t, years)
    
    
    plt.tight_layout()
    plt.subplots_adjust(top=0.89)
#    break
#    plt.savefig('plots/new_catalogue/Chi30Lightcurves/JHK_lightcurves/high_chi_JK/flux_'+str(obnum))#+str(varys['NUMBER_05B'][n])+'_lightcurve.png')
#    plt.close('all')
