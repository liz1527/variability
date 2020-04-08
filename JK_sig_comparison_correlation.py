#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  1 13:30:42 2018

Code to compute maximum likelihood variability

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

#%% Open the fits files and get data ###
### Import data for variables selected in K ###
varydata = Table.read('variable_tables/J_and_K_variables_varystats_DR11data.fits')
#varydata = varydata[varydata['KMAG_20']<=25]
xvarydata = varydata[varydata['X-ray']==True]
noxvarydata = varydata[varydata['X-ray']==False]

### Extract flux table and error tables ###
Kflux = varydata['Flux_K'] 
Kfluxerr = varydata['Fluxerr_K']
Jflux = varydata['Flux_J'] 
Jfluxerr = varydata['Fluxerr_J']

### Cross correlate the J and K band data ###
test = np.correlate(Kflux[0], Jflux[0], mode="full")

plt.plot(test)

end = time.time()
print(end-start)












