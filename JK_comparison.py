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
KKdata = fits.open('variable_tables/no06_variables_chi30_2arcsec_DR11data_restframe.fits')[1].data
KJdata = fits.open('variable_tables/no06_variables_chi30_2arcsec_DR11data_J_best.fits')[1].data

### Import data for variables selected in J ###
JKdata = fits.open('variable_tables/J_variables_chi40_noneg_DR11data_Kdata.fits')[1].data
JJdata = fits.open('variable_tables/J_variables_chi40_noneg_DR11data.fits')[1].data

### Import sig data ###
Jsigtb = Table.read('sigma_tables/quad_epoch_sigma_table_extra_clean_2arcsec_J.fits')
Ksigtb = Table.read('sigma_tables/quad_epoch_sigma_table_extra_clean_no06_2arcsec_neg.fits')

#Jxraydata = Jdata[Jdata['X-ray']==True]

#### Limit to Chandra region for simplicity ###
#KKdata = vari_funcs.chandra_only(KKdata)
#KJdata = vari_funcs.chandra_only(KJdata)
#JKdata = vari_funcs.chandra_only(JKdata)
#JJdata = vari_funcs.chandra_only(JJdata)

#%% plot z distributions of the sets ###
Jz = vari_funcs.get_z(JJdata)
Kz = vari_funcs.get_z(KKdata)

plt.figure()
plt.hist([Jz, Kz], bins=25, histtype='step', label=['J Variables', 'K Variables'])
plt.xlabel('z')
plt.legend(loc='upper center')


#%% plot M_star distributions of the sets ###
Jmstar = JJdata['Mstar_z_p'][Jz<6]
Kmstar = KKdata['Mstar_z_p'][Kz<6]

plt.figure()
plt.hist([Jmstar, Kmstar], bins=np.logspace(6,12,25), histtype='step', label=['J Variables', 'K Variables'])
plt.xlabel('Stellar Mass')
plt.xscale('log')
plt.legend(loc='upper left')

end = time.time()
print(end-start)















