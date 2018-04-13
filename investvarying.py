#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 15:11:24 2017

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
combined = fits.open('mag_flux_table_variables.fits')
tbdata = combined[1].data
chandra = fits.open('xray_visible_variables.fits')
chandata = chandra[1].data
nonchan = fits.open('non_xray_visible_variables.fits')
nonchandata = nonchan[1].data

### Create arrays of flux values but without 06B ###
flux = vari_funcs_no06.flux5_stacks(tbdata)
fluxchan = vari_funcs_no06.flux5_stacks(chandata) # for chandra non-stellar objects
fluxnonchan = vari_funcs_no06.flux5_stacks(nonchandata)

### Calculate average fluxes ###
avgfluxall = np.mean(flux, axis=0)
avgfluxchan = np.mean(fluxchan, axis=0)
avgfluxnon = np.mean(fluxnonchan, axis=0)

### Plot average light curves ###
vari_funcs_no06.avg_lightcurve(avgfluxall)
vari_funcs_no06.avg_lightcurve(avgfluxchan)
vari_funcs_no06.avg_lightcurve(avgfluxnon)
