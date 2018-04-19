#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 19 12:07:06 2018

@author: ppxee
"""

### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
from astropy.table import Table #for handling tables
import numpy as np #for handling arrays
#import math
import vari_funcs #my module to help run code neatly
plt.close('all') #close any open plots


tbdata = fits.open('mag_flux_tables/month_mag_flux_table_best.fits')[1].data

vari_funcs.lightcurve5months(203391, tbdata)