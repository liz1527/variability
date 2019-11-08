#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 14:53:15 2019

Code to show number of observations taken in each semester

@author: ppxee
"""

### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
from astropy.table import Table #for handling tables
import numpy as np #for handling arrays
import vari_funcs #my module to help run code neatly
plt.close('all')

semesters = ['05B', '06B', '07B', '08B', '09B', '10B', '11B', '12B']
framenum = np.zeros(len(semesters))
for n, sem in enumerate(semesters):
    list = open('semester_lists/K/'+sem+'_DR11_K.lst','r')
    lines = list.readlines()
    framenum[n] = len(lines)

x = [1,2,3,4,5,6,7,8]

plt.figure()
plt.bar(x, framenum)
plt.xticks(x, semesters)

### number of observations is number of frames/4 as 4 cameras ###
obsnum = framenum/4

plt.figure()
plt.bar(x, obsnum, color='k')
plt.xticks(x, semesters)
plt.xlabel('Semester')
plt.ylabel('Number of Observations')
plt.tight_layout()