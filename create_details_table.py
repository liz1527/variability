#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 11:55:57 2020

Create tables that save all the important info so it doesn't need to be 
calculated many times

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
import sys

### make filenames from parts (not sure if useful but trying it) ###
folder = 'variable_tables'
band = 'K'
neg = 'noneg'
filename = 'variables_no06_chi30_DR11data_DR11_photometry'
fullpath = folder+'/'+band+'/'+filename+'.fits'

aper = 4 #set to 2arcsec for now

### open orginal table ###
tb = fits.open(fullpath)[1].data

### Get flux and mag arrays ###
sigtb = Table.read('sigma_tables/quad_epoch_sigma_table_'+band+
                   '_extra_clean_2arcsec_'+neg+'.fits')
if band == 'K':    
    flux, fluxerr, tb = vari_funcs.k_mag_flux.create_quad_error_array(sigtb, 
                                                                      tb, 
                                                                      aper)
    mag = vari_funcs.k_mag_flux.mag_stacks(tb, aper)
elif band == 'J':
    flux, fluxerr, tb = vari_funcs.j_mag_flux.create_quad_error_array_J(sigtb, 
                                                                        tb, 
                                                                        aper)
    mag = vari_funcs.j_mag_flux.mag_stacks(tb, aper)
elif band == 'H':
    flux = vari_funcs.h_mag_flux.flux_stacks(tb, aper)
    fluxerr = vari_funcs.h_mag_flux.fluxerr_stacks(tb, aper)
    mag = vari_funcs.h_mag_flux.mag_stacks(tb, aper)
else:
    print('Band not recognised. Quitting script')
    sys.exit()
    
### Calc mean flux and mag ###
meanflux = np.nanmean(flux, axis=1)
meanmag = np.nanmean(mag, axis=1)

### find chi squared ###
chisq = vari_funcs.vary_stats.my_chisquare_err(flux, fluxerr)

### get x-ray lum and mixed z column ###
xflux, xL, z = vari_funcs.xray_funcs.get_xray_L_mixedtable(tb) # this gives upp limit for not detected sourcs

### get stellarity ###
