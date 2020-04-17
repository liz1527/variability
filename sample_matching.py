#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 12:18:07 2020

Code to create matched samples of X-ray and non-X-ray variables

@author: ppxee
"""


import time
start = time.time()
print(start)

### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
from astropy.table import Table #for handling tables
import numpy as np #for handling arrays
#import math
#from astropy.stats import median_absolute_deviation
import vari_funcs #my module to help run code neatly
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
plt.close('all') #close any open plots
#from numpy.lib.recfunctions import append_fields

#%% Open the fits files and split by X-ray detection ###
varydata = Table.read('variable_tables/J_and_K_variables_varystats_DR11data.fits')
xvarydata = varydata[varydata['X-ray']==True]
noxvarydata = varydata[varydata['X-ray']==False]

#%% get data arrays so don't extract on every loop ###
noxsize = noxvarydata['FWHM_WORLD'] # FWHM good estimate of size
noxlum = noxvarydata['M_K_z_p'] # Absolute mag = luminosity
noxz = noxvarydata['z_p'] # best photometric z so consistant with lum calcs

#%% Iterate over objects ###
for ob in xvarydata:
    ### Extract arrays ###
    obsize = ob['FWHM_WORLD'] # FWHM good estimate of size
    oblum = ob['M_K_z_p'] # Absolute mag = luminosity
    obz = ob['z_p'] # best photometric z so consistant with lum calcs
    
    ### Find differences ###
    sizediff = noxsize - obsize
    lumdiff = noxlum - oblum
    zdiff = noxz - obz
    
    ### Start with size match ###