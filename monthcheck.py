#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 13 14:46:27 2018

Code to check the output of the initial month stacks

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

tbdata = fits.open('mag_flux_tables/month_mag_flux_table.fits')[1].data

mag = np.stack(([tbdata['MAG_APER_sep05'][:,5], tbdata['MAG_APER_oct05'][:,5],
            tbdata['MAG_APER_nov05'][:,5], tbdata['MAG_APER_dec05'][:,5],
            tbdata['MAG_APER_jan06'][:,5], tbdata['MAG_APER_jan07'][:,5]]), axis=1)

mag, tbdata = vari_funcs.no99(mag, tbdata)

avgflux = np.mean(mag, axis=0)

plt.plot(avgflux, 'o')
plt.hlines(np.mean(mag), 0, 5)