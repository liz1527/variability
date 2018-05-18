#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 20 12:16:04 2017

@author: ppxee
"""
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
from astropy.table import Table #for handling tables
import numpy as np #for handling arrays
#import math
#from astropy.stats import median_absolute_deviation
import vari_funcs #my module to help run code neatly
plt.close('all') #close any open plots


### Open the fits files and get data ###
tbdata = fits.open('mag_flux_tables/mag_flux_table_best.fits')[1].data
chandata = fits.open('mag_flux_tables/xray_mag_flux_table_best.fits')[1].data
sdata = fits.open('mag_flux_tables/stars_mag_flux_table.fits')[1].data

### extract magnitude arrays ###
allmag = vari_funcs.flux5_stacks(tbdata)
allchanmag = vari_funcs.flux5_stacks(chandata) # for chandra non-stellar objects
allsmag = vari_funcs.flux5_stacks(sdata)

### remove 99s ###
allmag, tbdata = vari_funcs.noneg(allmag, tbdata)
allchanmag, chandata = vari_funcs.noneg(allchanmag, chandata)
allsmag, sdata = vari_funcs.noneg(allsmag, sdata)

# get error data
magerr = vari_funcs.fluxerr5_stacks(tbdata)

#normalise
normmag = vari_funcs.normalise_flux(allmag)
normmagerr = vari_funcs.err_correct_flux(allmag, magerr)

vari_funcs.avg_lightcurve(allmag[100,:], magerr[100,:])
vari_funcs.avg_lightcurve(normmag[100,:], normmagerr[100,:])