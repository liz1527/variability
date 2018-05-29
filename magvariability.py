#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 15:04:52 2017

Variability using magnitudes instead of fluxes

@author: ppxee
"""

### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
#from astropy.table import Table #for handling tables
import numpy as np #for handling arrays
#import math
from astropy.stats import median_absolute_deviation
import vari_funcs #my module to help run code neatly
#plt.close('all') #close any open plots
#from numpy.lib.recfunctions import append_fields

### Open the fits files and get data ###
combined = fits.open('mag_flux_tables/mag_flux_table_best.fits')
tbdata = combined[1].data
chandra = fits.open('mag_flux_tables/xray_mag_flux_table_best.fits')
chandata = chandra[1].data
stars = fits.open('mag_flux_tables/stars_mag_flux_table.fits')
sdata = stars[1].data

### Restrict objects to those in the Chandra field ###
#tbdata = vari_funcs.chandra_only(tbdata)
#chandata = vari_funcs.chandra_only(chandata)
#sdata = vari_funcs.chandra_only(sdata)

#### Split chandra and xmm ###
#mask = np.isnan(chandata['RA'])
#xmmdata = chandata[mask]
#chandata = chandata[~mask]

## Create arrays of flux values ###
fluxn = vari_funcs.mag5_stacks(tbdata)
fluxchann = vari_funcs.mag5_stacks(chandata) # for chandra non-stellar objects
sfluxn = vari_funcs.mag5_stacks(sdata)
#fluxxmm = vari_funcs.mag5_stacks(xmmdata)

### remove values that are +/-99 ###
fluxn, tbdata = vari_funcs.no99(fluxn, tbdata)
fluxchann, chandata = vari_funcs.no99(fluxchann, chandata)
sfluxn, sdata = vari_funcs.no99(sfluxn, sdata)
#fluxxmm, xmmdata = vari_funcs.no99(fluxxmm, xmmdata)

#### normalise to 1 ###
#fluxn = vari_funcs.normalise(fluxn)
#fluxchann = vari_funcs.normalise(fluxchann)
#sfluxn = vari_funcs.normalise(sfluxn)

#
### Multiply all flux values in a yearstack by the correct constant ###
#fluxcorrn = vari_funcs.psf_correct_mag(fluxn, fluxn, 'median') 
#fluxchancorrn = vari_funcs.psf_correct(fluxn, fluxchann, 'median') 
#sfluxcorrn = vari_funcs.psf_correct(fluxn, sfluxn, 'median') 

### get error arrays and correct them ###
fluxerrn = vari_funcs.magerr5_stacks(tbdata)
fluxerrchan = vari_funcs.magerr5_stacks(chandata)
sfluxerr = vari_funcs.magerr5_stacks(sdata)
#fluxerrcorrn = vari_funcs.err_correct(fluxn, fluxerrn, fluxcorrn)


fig = vari_funcs.flux_variability_plot(fluxn, fluxchann, 'mad',
                                            starflux=sfluxn, stars=True)

#fig = vari_funcs.flux_variability_plot(fluxn, fluxchann, 'excess',
#                                       fluxerr = fluxerrn, starfluxerr = sfluxerr,
#                                            starflux=sfluxn, stars=True, 
#                                            chanerr = fluxerrchan,
#                                            normalised=True)

#fig.canvas.mpl_connect('pick_event', vari_funcs.onpick)
#

combined.close()
chandra.close()
stars.close()