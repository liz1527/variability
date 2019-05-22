#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 27 11:49:02 2018

@author: ppxee
"""

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
plt.close('all')

xray = fits.open('variable_tables/xray_variable_table_details_report.fits')[1].data
nonxray = fits.open('variable_tables/non_xray_variable_table_details_report.fits')[1].data

magx = xray['KMAG_20']
magn = nonxray['KMAG_20']
bins = np.arange(17,23,0.2)

plt.hist(magn, bins, color='b', zorder=2,alpha=0.5, label='Non X-ray')
plt.hist(magx, bins, color='r',alpha=0.5, label='X-ray')
plt.xlabel('K band magnitude')
plt.ylabel('Number')
plt.legend()