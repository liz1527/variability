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
plt.close('all') #close any open plots
#from numpy.lib.recfunctions import append_fields

### Open arrays ###
fluxn = np.load('mag_arrays/monthmagarray.npy')
fluxerrn = np.load('mag_arrays/monthmagerrdiffarray.npy')
fluxchann = np.load('mag_arrays/monthxraymagarray.npy')
fluxerrchan = np.load('mag_arrays/monthxraymagerrdiffarray.npy')
sfluxn = np.load('mag_arrays/monthstarmagarray.npy')
sfluxerr = np.load('mag_arrays/monthstarmagerrdiffarray.npy')


#%%
## Create Plot for non corrected ###
fig1 = vari_funcs.flux_variability_plot(fluxn, fluxchann, 'mad',
                                            starflux=sfluxn, stars=True)
fig2 = vari_funcs.flux_variability_plot(fluxn, fluxchann, 'excess',
                                       fluxerr = fluxerrn, 
                                       starfluxerr = sfluxerr,
                                            starflux=sfluxn, stars=True,
                                            chanerr = fluxerrchan,
                                            normalised=True)
fig2.canvas.mpl_connect('pick_event', vari_funcs.onpickmonth)
fig1.canvas.mpl_connect('pick_event', vari_funcs.onpickmonth)

