#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  8 16:44:00 2018

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
plt.close('all') #close any open plots

font = {'family' : 'DejaVu Sans',
        'weight' : 'normal',
        'size'   : 14}

plt.rc('font', **font)

### Open the fits files and get data ###
tbdata = fits.open('mag_flux_tables/mag_flux_table_best.fits')[1].data
varys = fits.open('variable_tables/non_xray_variables_with_spectra_for_Ismael.fits')[1].data

#obnum = input('Object ID: ')#252446 158984
obnums = varys['NUMBER_05B']
dr11 = varys['ID_DR11']

for n, obnum in enumerate(obnums):
    vari_funcs.lightcurve5(int(obnum), tbdata)
    plt.tight_layout()
    axes = plt.gca()
    ylims = axes.get_ylim()
    ymid = (ylims[1]+ylims[0])/2
#    plt.ylim(ymin=ymid-0.26, ymax=ymid+0.26)
    plt.savefig('IsmealLightcurves/DR11_'+str(dr11[n])+'_lightcurve.png')
