#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 10:34:49 2019

Code to create a venn diagram describing the selection methods and their 
overlaps.

@author: ppxee
"""

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
plt.close('all')

from matplotlib_venn import venn2

varydata = Table.read('variable_tables/J_and_K_variables_varystats_DR11data.fits')

### Split table into J and K selected ###
jdata = varydata[varydata['Chi_J'] > 32.08]
kdata = varydata[varydata['Chi_K'] > 30]

### Find overlap ###
jIDs = jdata['ID']
kIDs = kdata['ID']
jmask = np.isin(jIDs, kIDs)
kmask = np.isin(kIDs, jIDs)

### Mask tables ###
bothdata = jdata[jmask]
jdata = jdata[~jmask]
kdata = kdata[~kmask]

### get lengths ###
jnum = len(jdata)
knum = len(kdata)
bothnum = len(bothdata)
tot = jnum + knum + bothnum
print(tot)

### Make Venn diagram ###
plt.figure()
subsets = {
        '10': jnum, # just J
        '01': knum, # just K
        '11': bothnum # J + K
           }
values = venn2(subsets = subsets, set_labels = ('J', 'K'), set_colors=('green', 'purple'))
plt.show()
plt.tight_layout()