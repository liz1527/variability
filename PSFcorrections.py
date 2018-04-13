#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 11:55:33 2017

@author: ppxee
"""

### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
#from astropy.table import Table #for handling tables
import numpy as np #for handling arrays
#import math
#from astropy.stats import median_absolute_deviation
import vari_funcs #my module to help run code neatly
plt.close('all') #close any open plots

### Open the fits files and get data ###
combined = fits.open('mag_flux_table_best.fits')
tbdata = combined[1].data
chandra = fits.open('chandra_mag_flux_table_best.fits')
chandata = chandra[1].data
stars = fits.open('stars_mag_flux_table.fits')
sdata = stars[1].data
#chanposstars = fits.open('chandra_possible_stars.fits')
#pstardata = chanposstars[1].data

### Restrict objects to those in the Chandra field ###
tbdata = vari_funcs.chandra_only(tbdata)

## Create arrays of flux values ###
flux = vari_funcs.flux5_stacks(tbdata)
fluxchan = vari_funcs.flux5_stacks(chandata) # for chandra non-stellar objects
#fluxpstar = vari_funcs.flux5_stacks(pstardata) # for chandra objects identified as stars
sflux = vari_funcs.flux5_stacks(sdata)

### Calculate the average flux for each year within each flux bin ###
avgfluxperepoch = np.mean(flux, axis=0)#for UDS
avgfluxchanperepoch = np.mean(fluxchan, axis=0) #for non-stellar chandra
savgfluxperepoch = np.median(sflux, axis=0) #median so randomly bright star shouldn't skew it

### Multiply all flux values in a yearstack by the correct constant ###
fluxcorr = vari_funcs.psf_correct(flux, flux, 'mean') 
fluxchancorr = vari_funcs.psf_correct(flux, fluxchan, 'mean') 
avgfluxcorr = vari_funcs.psf_correct(flux, avgfluxperepoch, 'mean') 

### Multiply all flux values in a yearstack by the star constant ###
fluxscorr = vari_funcs.psf_correct(sflux, flux, 'median') 
fluxchanscorr = vari_funcs.psf_correct(sflux, fluxchan, 'median') 
avgfluxscorr = vari_funcs.psf_correct(sflux, avgfluxperepoch, 'median') 

### get error arrays and correct them ###
fluxerr = vari_funcs.fluxerr5_stacks(tbdata)
fluxerrcorr = vari_funcs.err_correct(flux, fluxerr, fluxcorr)
fluxerrscorr = vari_funcs.err_correct(flux, fluxerr, fluxscorr)

## Create arrays of flux values but without 06B ###
fluxn = vari_funcs.flux5_stacks_no06(tbdata)
fluxchann = vari_funcs.flux5_stacks_no06(chandata) # for chandra non-stellar objects
#fluxpstar = vari_funcs.flux5_stacks(pstardata) # for chandra objects identified as stars
sfluxn = vari_funcs.flux5_stacks_no06(sdata)

### Calculate the average flux for each year within each flux bin ###
avgfluxperepochn = np.mean(fluxn, axis=0)#for UDS
avgfluxchanperepochn = np.mean(fluxchann, axis=0) #for non-stellar chandra
savgfluxperepochn = np.median(sfluxn, axis=0) #median so randomly bright star shouldn't skew it

### Multiply all flux values in a yearstack by the correct constant ###
fluxcorrn = vari_funcs.psf_correct(fluxn, fluxn, 'mean') 
fluxchancorrn = vari_funcs.psf_correct(fluxn, fluxchann, 'mean') 
avgfluxcorrn = vari_funcs.psf_correct(fluxn, avgfluxperepochn, 'mean') 

### Multiply all flux values in a yearstack by the star constant ###
fluxscorrn = vari_funcs.psf_correct(sfluxn, fluxn, 'median') 
fluxchanscorrn = vari_funcs.psf_correct(sfluxn, fluxchann, 'median') 
avgfluxscorrn = vari_funcs.psf_correct(sfluxn, avgfluxperepochn, 'median') 

### get error arrays and correct them ###
fluxerrn = vari_funcs.fluxerr5_stacks_no06(tbdata)
fluxerrcorrn = vari_funcs.err_correct(fluxn, fluxerrn, fluxcorrn)
fluxerrscorr = vari_funcs.err_correct(fluxn, fluxerrn, fluxscorrn)

fig = vari_funcs.flux_variability_plot(fluxcorrn, fluxchancorrn, 'mad', normalised = True)

fig.canvas.mpl_connect('pick_event', vari_funcs.onpickchanonly)

