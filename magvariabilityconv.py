#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 19 10:56:37 2017

Variability using magnitudes on the output from the convolved images

@author: ppxee
"""

### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
#from astropy.table import Table #for handling tables
import numpy as np #for handling arrays
#import math
from astropy.stats import median_absolute_deviation
import vari_funcs_no06 #my module to help run code neatly
plt.close('all') #close any open plots

### Open the fits files and get data ###
combined = fits.open('mag_flux_table_best_qS.fits')
tbdata = combined[1].data
chandra = fits.open('xray_mag_flux_table_best_qS.fits')
chandata = chandra[1].data
stars = fits.open('stars_mag_flux_table_qS.fits')
sdata = stars[1].data

#### Restrict objects to those in the Chandra field ###
#tbdata = vari_funcs_no06.chandra_only(tbdata)
#sdata = vari_funcs_no06.chandra_only(sdata)

## Create arrays of flux values but without 06B ###
fluxn = vari_funcs_no06.mag5_stacks(tbdata)
fluxchann = vari_funcs_no06.mag5_stacks(chandata) # for chandra non-stellar objects
sfluxn = vari_funcs_no06.mag5_stacks(sdata)

### remove values that are +/-99 ###
fluxn, tbdata = vari_funcs_no06.no99(fluxn, tbdata)
fluxchann, chandata = vari_funcs_no06.no99(fluxchann, chandata)
sfluxn, sdata = vari_funcs_no06.no99(sfluxn, sdata)

### Calculate the average flux for each year within each flux bin ###
#avgfluxperepochn = np.mean(fluxn, axis=0)#for UDS
#avgfluxchanperepochn = np.mean(fluxchann, axis=0) #for non-stellar chandra
#avgfluxperepochn = np.mean(fluxchann, axis=0) #for non-stellar chandra

### Multiply all flux values in a yearstack by the correct constant ###
#fluxcorrn = vari_funcs_no06.psf_correct(fluxn, fluxn, 'median') 
#fluxchancorrn = vari_funcs_no06.psf_correct(fluxn, fluxchann, 'median') 
#sfluxcorrn = vari_funcs_no06.psf_correct(fluxn, sfluxn, 'median') 

### get error arrays and correct them ###
#fluxerrn = vari_funcs_no06.magerr5_stacks(tbdata)
#fluxerrn = fluxerrn[mask]
#fluxerrcorrn = vari_funcs_no06.err_correct(fluxn, fluxerrn, fluxcorrn)

### Create Plot ###
#fig = vari_funcs_no06.flux_variability_plot(fluxn, fluxchann, 'mad',
#                                            starflux=sfluxn, stars=True)
#vari_funcs_no06.flux_variability_plot(fluxcorrn, fluxchancorrn, 'mad', 
#                                            starflux=sfluxcorrn, stars=True)
#plt.title('Using Convolved Images')
#fig.canvas.mpl_connect('pick_event', vari_funcs_no06.onpickchanonly)
#fig = vari_funcs_no06.flux_variability_plot(fluxn, fluxchann, 'mad',
#                                            starflux=sfluxn, stars=True,
#                                            normalised=True, psfcorrect=True)
#plt.title('Convolved Normalised and PSF-Corrected')
#fig = vari_funcs_no06.flux_variability_plot(fluxn, fluxchann, 'mad',
#                                            starflux=sfluxn, stars=True,
#                                            normalised=True)
#plt.title('Convolved Normalised')
#fig = vari_funcs_no06.flux_variability_plot(fluxn, fluxchann, 'mad',
#                                            starflux=sfluxn, stars=True,
#                                            psfcorrect=True)
#plt.title('Convolved PSF-Corrected')
fig = vari_funcs_no06.flux_variability_plot(fluxn, fluxchann, 'mad',
                                            starflux=sfluxn, stars=True)
plt.title('Convolved No Changes')

### close fits files ###
combined.close()
chandra.close()
stars.close()