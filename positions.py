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
plt.style.use('dark_background')

### Get fits files ###
chandata = fits.open('variable_tables/no06_variables_chi30_2arcsec_chandata_DR11data_restframe.fits')[1].data
xmmdata = fits.open('variable_tables/no06_variables_chi30_2arcsec_xmmdata_DR11data_restframe.fits')[1].data
tbdata = fits.open('variable_tables/no06_variables_chi30_2arcsec_DR11data_restframe.fits')[1].data

### Remove edges ###
tbdata = vari_funcs.remove_edges(tbdata)
chandata = vari_funcs.remove_edges(chandata)
xmmdata = vari_funcs.remove_edges(xmmdata)

nontb = tbdata[~tbdata['X-ray']]
#xtb = tbdata[tbdata['X-ray']]

### Get positions ###
nonRA = nontb['ALPHA_J2000']
nonDec = nontb['DELTA_J2000']
#xRA = xtb['ALPHA_J2000']
#xDec = xtb['DELTA_J2000']
chanRA = chandata['ALPHA_J2000']
chanDec = chandata['DELTA_J2000']
xmmRA = xmmdata['ALPHA_J2000']
xmmDec = xmmdata['DELTA_J2000']

### Plot positions ###
plt.figure(figsize=[8,8])
plt.plot(nonRA, nonDec, 'bo')
#plt.plot(xRA, xDec, 'ro', mfc='None', markersize=10)
plt.plot(chanRA, chanDec, 'ro')
plt.plot(xmmRA, xmmDec, 'ro')
plt.xlim(xmin=34, xmax=34.9)
plt.ylim(ymin=-5.55, ymax=-4.65)
#plt.xlim(xmin=34.06, xmax=34.87)
#plt.ylim(ymin=-5.54, ymax=-4.7)
plt.gca().invert_xaxis()
plt.xlabel('RA')
plt.ylabel('Dec')
#plt.savefig('plots/Chi_variables/no06_extra_clean/positions_all.png')

print('Full sample')
print('Total = '+str(len(tbdata)))
print('X-ray = '+str(len(chanRA)+len(xmmRA)))
##%% Restrict to >1e4 Average Flux ###
#
### Create arrays of flux values ###
#allflux = vari_funcs.flux5_stacks(tbdata)
#flux = vari_funcs.flux5_stacks(nontb)
#chanflux = vari_funcs.flux5_stacks(chandata) 
#xmmflux = vari_funcs.flux5_stacks(xmmdata)
#
#### Get average flux ###
#avgallflux = np.nanmean(allflux, axis=1)
#avgflux = np.nanmean(flux, axis=1)
#avgchanflux = np.nanmean(chanflux, axis=1)
#avgxmmflux = np.nanmean(xmmflux, axis=1)
#
#hightbdata = tbdata[avgallflux > 1e4]
#highnontb = nontb[avgflux > 1e4]
#highchandata = chandata[avgchanflux > 1e4]
#highxmmdata = xmmdata[avgxmmflux > 1e4]
#
#### Get positions ###
#nonRA = highnontb['ALPHA_J2000']
#nonDec = highnontb['DELTA_J2000']
#chanRA = highchandata['ALPHA_J2000']
#chanDec = highchandata['DELTA_J2000']
#xmmRA = highxmmdata['ALPHA_J2000']
#xmmDec = highxmmdata['DELTA_J2000']
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
#print('High Flux')
#print('Total = '+str(len(hightbdata)))
#print('X-ray = '+str(len(chanRA)+len(xmmRA)))
#%% Restrict to Chandra region ###
tbdata = vari_funcs.chandra_only(tbdata)
#chandata = vari_funcs.chandra_only(chandata)
xmmdata = vari_funcs.chandra_only(xmmdata)

nontb = tbdata#[~tbdata['X-ray']]

### Get positions ###
nonRA = nontb['ALPHA_J2000']
nonDec = nontb['DELTA_J2000']
chanRA = chandata['ALPHA_J2000']
chanDec = chandata['DELTA_J2000']
xmmRA = xmmdata['ALPHA_J2000']
xmmDec = xmmdata['DELTA_J2000']

### Plot positions ###
plt.figure(figsize=[8,8])
plt.plot(nonRA, nonDec, 'bo')
plt.plot(chanRA, chanDec, 'ro')
#plt.plot(xmmRA, xmmDec, 'ro')
plt.xlim(xmin=34, xmax=34.9)
plt.ylim(ymin=-5.55, ymax=-4.65)
plt.gca().invert_xaxis()
plt.xlabel('RA')
plt.ylabel('Dec')
#plt.savefig('plots/Chi_variables/no06_extra_clean/positions_chandra.png')
#
#print('In Chandra')
#print('Total = '+str(len(tbdata)))
#print('X-ray = '+str(len(chanRA)+len(xmmRA)))
#
### with only high flux ###
#hightbdata = vari_funcs.chandra_only(hightbdata)
#highnontb = vari_funcs.chandra_only(highnontb)
##highchandata = vari_funcs.chandra_only(highchandata)
#highxmmdata = vari_funcs.chandra_only(highxmmdata)
#
#### Get positions ###
#nonRA = highnontb['ALPHA_J2000']
#nonDec = highnontb['DELTA_J2000']
#chanRA = highchandata['ALPHA_J2000']
#chanDec = highchandata['DELTA_J2000']
#xmmRA = highxmmdata['ALPHA_J2000']
#xmmDec = highxmmdata['DELTA_J2000']
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
#plt.ylabel('Dec')
##plt.savefig('plots/Chi_variables/no06_extra_clean/positions_chandra_highflux.pdf')
#
#print('In Chandra and high flux')
#print('Total = '+str(len(hightbdata)))
#print('X-ray = '+str(len(chanRA)+len(xmmRA)))