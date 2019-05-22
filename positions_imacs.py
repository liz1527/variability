#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 11:42:35 2019

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
plt.close('all')
plt.style.use('default')

### Get fits files ###
#xdata = fits.open('variable_tables/no06_variables_chi30_nospectra_xray.fits')[1].data
#tbdata = fits.open('variable_tables/no06_variables_chi30_nospectra.fits')[1].data
#xdata = fits.open('variable_tables/no06_variables_chi30_noUDSz_xray.fits')[1].data
#tbdata = fits.open('variable_tables/no06_variables_chi30_noUDSz.fits')[1].data
chandata = fits.open('variable_tables/no06_variables_chi30_chandata_DR11data_restframe.fits')[1].data
xmmdata = fits.open('variable_tables/no06_variables_chi30_xmmdata_DR11data_restframe.fits')[1].data
tbdata = fits.open('variable_tables/no06_variables_chi30_DR11data_restframe.fits')[1].data

### Remove edges ###
tbdata = vari_funcs.remove_edges(tbdata)
#xdata = vari_funcs.remove_edges(xdata)

chandata = vari_funcs.remove_edges(chandata)
xmmdata = vari_funcs.remove_edges(xmmdata)

nontb = tbdata[~tbdata['X-ray']]
#xtb = tbdata[tbdata['X-ray']]

### Apply a magnitude cut ###
hightbdata = tbdata[tbdata['iMAG_20']<24]
highnontb = nontb[nontb['iMAG_20']<24]
#xdata = xdata[xdata['iMAG_20']<24]
#
highchandata = chandata[chandata['iMAG_20']<24]
highxmmdata = xmmdata[xmmdata['iMAG_20']<24]

#tbdata = tbdata[tbdata['RMAG_20']<24.5]
#nontb = nontb[nontb['RMAG_20']<24.5]
#xdata = xdata[xdata['RMAG_20']<24.5]



### Get positions ###
nonRA = nontb['ALPHA_J2000']
nonDec = nontb['DELTA_J2000']
#xRA = xdata['ALPHA_J2000']
#xDec = xdata['DELTA_J2000']

chanRA = chandata['ALPHA_J2000']
chanDec = chandata['DELTA_J2000']
xmmRA = xmmdata['ALPHA_J2000']
xmmDec = xmmdata['DELTA_J2000']
xRA = np.append(chanRA, xmmRA)
xDec = np.append(chanDec, xmmDec)

highnonRA = highnontb['ALPHA_J2000']
highnonDec = highnontb['DELTA_J2000']
#xRA = xdata['ALPHA_J2000']
#xDec = xdata['DELTA_J2000']

highchanRA = highchandata['ALPHA_J2000']
highchanDec = highchandata['DELTA_J2000']
highxmmRA = highxmmdata['ALPHA_J2000']
highxmmDec = highxmmdata['DELTA_J2000']
highxRA = np.append(highchanRA, highxmmRA)
highxDec = np.append(highchanDec, highxmmDec)

### Plot positions ###
plt.figure(figsize=[8,8])
plt.plot(highnonRA, highnonDec, 'bo',markersize=10)
plt.plot(highxRA, highxDec, 'ro',markersize=10)
plt.plot(nonRA, nonDec, 'b+',alpha=0.5,zorder=0)
plt.plot(xRA, xDec, 'r+',alpha=0.5,zorder=0)
plt.xlim(xmin=34, xmax=34.9)
plt.ylim(ymin=-5.6, ymax=-4.65)
#plt.xlim(xmin=34.06, xmax=34.87)
#plt.ylim(ymin=-5.54, ymax=-4.7)
plt.gca().invert_xaxis()
plt.xlabel('RA')
plt.ylabel('Dec')

### Plot on the IMACS FoV ###
diam = 27.4/60 # 27.4 arcmin
circle1=  plt.Circle((34.645,-4.88), radius=diam/2, fill=False, color='k')
plt.gca().add_artist(circle1)
circle2=  plt.Circle((34.23,-4.98), radius=diam/2, fill=False, color='k')
plt.gca().add_artist(circle2)
circle3=  plt.Circle((34.65,-5.25), radius=diam/2, fill=False, color='k')
plt.gca().add_artist(circle3)
circle4=  plt.Circle((34.34,-5.33), radius=diam/2, fill=False, color='k')
plt.gca().add_artist(circle4)
plt.tight_layout()

print('Full sample')
print('Total = '+str(len(tbdata)))
print('X-ray = '+str(len(xRA)))

### plot r band hist of sample ###
#mask1 = tbdata['RMAG_20']!=99
#rmag = tbdata['RMAG_20']#[mask1]
#mask1 = xdata['RMAG_20']!=99
#xrmag = xdata['RMAG_20'][mask1]

#xrayrmag = np.append(chanrmag,xmmrmag)
#med = np.median(nonrmag)
#plt.figure()
#plt.hist(nonrmag, color='b', label = 'Non X-ray',histtype='step')
#plt.hist(xrayrmag, color='r', label='X-ray',histtype='step')
#plt.xlabel('R Band Magnitude')
#plt.ylabel('Counts')
#plt.legend()
#plt.tight_layout()
###%% Restrict to >1e4 Average Flux ###
##
#### Create arrays of flux values ###
##allflux = vari_funcs.flux5_stacks(tbdata)
##flux = vari_funcs.flux5_stacks(nontb)
##chanflux = vari_funcs.flux5_stacks(chandata) 
##xmmflux = vari_funcs.flux5_stacks(xmmdata)
##
##### Get average flux ###
##avgallflux = np.nanmean(allflux, axis=1)
##avgflux = np.nanmean(flux, axis=1)
##avgchanflux = np.nanmean(chanflux, axis=1)
##avgxmmflux = np.nanmean(xmmflux, axis=1)
##
##hightbdata = tbdata[avgallflux > 1e4]
##highnontb = nontb[avgflux > 1e4]
##highchandata = chandata[avgchanflux > 1e4]
##highxmmdata = xmmdata[avgxmmflux > 1e4]
##
##### Get positions ###
##nonRA = highnontb['ALPHA_J2000']
##nonDec = highnontb['DELTA_J2000']
##chanRA = highchandata['ALPHA_J2000']
##chanDec = highchandata['DELTA_J2000']
##xmmRA = highxmmdata['ALPHA_J2000']
##xmmDec = highxmmdata['DELTA_J2000']
#
#### Plot positions ###
#plt.figure(figsize=[8,8])
#plt.plot(nonRA, nonDec, 'bo')
#plt.plot(chanRA, chanDec, 'ro')
#plt.plot(xmmRA, xmmDec, 'ro')
#plt.xlim(xmin=34, xmax=34.9)
#plt.ylim(ymin=-5.6, ymax=-4.65)
#plt.gca().invert_xaxis()
#plt.xlabel('RA')
##plt.savefig('plots/Chi_variables/no06_extra_clean/positions_highflux.pdf')
#
#### Plot on the IMACS FoV ###
#diam = 27.4/60 # 27.4 arcmin
#circle1=  plt.Circle((34.64,-4.93), radius=diam/2, fill=False, color='k')
#plt.gca().add_artist(circle1)
#circle2=  plt.Circle((34.23,-4.93), radius=diam/2, fill=False, color='k')
#plt.gca().add_artist(circle2)
#circle3=  plt.Circle((34.64,-5.32), radius=diam/2, fill=False, color='k')
#plt.gca().add_artist(circle3)
#circle4=  plt.Circle((34.23,-5.32), radius=diam/2, fill=False, color='k')
#plt.gca().add_artist(circle4)
##print('High Flux')
#print('Total = '+str(len(hightbdata)))
#print('X-ray = '+str(len(chanRA)+len(xmmRA)))

#### plot r band hist of high flux sample ###
#mask1 = hightbdata['RMAG_20']!=99
#rmag = hightbdata['RMAG_20'][mask1]
#mask1 = highnontb['RMAG_20']!=99
#nonrmag = highnontb['RMAG_20'][mask1]
#mask1 = highchandata['RMAG_20']!=99
#chanrmag = highchandata['RMAG_20'][mask1]
#mask1 = highxmmdata['RMAG_20']!=99
#xmmrmag = highxmmdata['RMAG_20'][mask1]
#
#xrayrmag = np.append(chanrmag,xmmrmag)
#plt.figure()
#plt.hist(nonrmag, color='b', label = 'Non X-ray',histtype='step')
#plt.hist(xrayrmag, color='r', label='X-ray',histtype='step')
#plt.xlabel('R Band Magnitude')
#plt.ylabel('Counts')
#plt.legend()
#plt.tight_layout()