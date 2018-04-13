#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 16:36:21 2017

Code to plot variability against luminosity

@author: ppxee
"""

### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
#from astropy.table import Table #for handling tables
import numpy as np #for handling arrays
#import math
#from astropy.stats import median_absolute_deviation
import vari_funcs_no06 #my module to help run code neatly
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
plt.close('all') #close any open plots

### Open the fits files and get data ###
combined = fits.open('variable_sources_data.fits')
tbdata = combined[1].data
chandra = fits.open('chandra_variable_sources_data.fits')
chandata = chandra[1].data

### Define cosmology ###
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

### Find luminosity distance both for all and just for chandra sources ###
z = tbdata['z']
DL = cosmo.luminosity_distance(z)
DL = DL.to(u.cm)
chanz = chandata['z']
chanDL = cosmo.luminosity_distance(chanz)
chanDL = chanDL.to(u.cm)

### Calculate the luminosity ###
xrayF = chandata['Full_flux']
xrayL = xrayF*4*np.pi*chanDL.value**2

### Plot against MAD ###
mad = chandata['MAD']
plt.figure()
plt.scatter(xrayL, mad, c=chanz)
#plt.plot(xrayL, mad, 'o')
clb = plt.colorbar()
clb.set_label('redshift')
plt.xlabel('X-ray Luminosity (ergs/s)')
plt.ylabel('MAD')

plt.title('Amplitude of variability against x-ray luminosity')

### Plot redshift against MAD ###
mad = tbdata['MAD']
plt.figure()
plt.plot(z, mad, 'o')
plt.xlabel('redshift')
plt.ylabel('MAD')
plt.title('Amplitude of variability against redshift')

combined.close()
chandra.close()
