#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 20 12:16:04 2017

@author: ppxee
"""
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
from astropy.table import Table #for handling tables
import numpy as np #for handling arrays
#import math
#from astropy.stats import median_absolute_deviation
import vari_funcs #my module to help run code neatly
plt.close('all') #close any open plots


values = [-7, -6, -5, -4, -3, -2, -1,  0,  1,  2,  3,  4,  5,  6,  7]
weights = [-0.01168533, -0.01891393, -0.00103619,  0.00439063,  0.00318069,
        0.00221236,  0.00222293,  0.00202136,  0.00244354,  0.00285592,
        0.00347878,  0.00448091,  0.00568668,  0.00546587, -0.00630108]

s = 0
for x, y in zip(values, weights):
    s += x * y

average = s / sum(weights)
print(average)