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

from matplotlib_venn import venn3
plt.figure()
subsets = {
        '100': 185, # just NIR
        '010': 725, # just X-ray
        '001': 5874, # just stern
        '110': 20, # NIR + X-ray
        '011': 355, # X-ray + stern
        '101': 41, # NIR + Stern
        '111': 147 # all
           }
values = venn3(subsets = subsets, 
      set_labels = ('NIR Variability', 'X-ray', 'Mid-IR'), normalize_to=10)
plt.show()
plt.tight_layout()