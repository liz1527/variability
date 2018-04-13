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
combined = fits.open('mag_flux_tables/mag_flux_table_best_qS.fits')
tbdata = combined[1].data

#obnum = input('Object ID: ')#252446 158984
obnums = [67178, 227399, 203391, 187967,193325,243208]
for obnum in obnums:
    vari_funcs.lightcurve5(int(obnum), tbdata)
    plt.tight_layout()
    axes = plt.gca()
    ylims = axes.get_ylim()
    ymid = (ylims[1]+ylims[0])/2
    plt.ylim(ymin=ymid-0.26, ymax=ymid+0.26)
    #plt.savefig(str(obnum)+'.png')
