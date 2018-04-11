#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 11 10:10:23 2018

Investigate properties of the variable sources

@author: ppxee
"""

### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
from astropy.table import Table #for handling tables
import numpy as np #for handling arrays
#import math
from astropy.stats import median_absolute_deviation
import vari_funcs #my module to help run code neatly
# Imports for catalogue matching
from astropy.coordinates import match_coordinates_sky
from astropy.coordinates import SkyCoord
from astropy import units as u
# Imports for luminosity distance
from astropy.cosmology import FlatLambdaCDM
# Set up cosmology
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
plt.close('all') #close any open plots

# Import fits table
varydata = fits.open('variable_tables/variable_with06B_mag21_DR11details.fits')[1].data

# Extract magnitude table and error table
mag = vari_funcs.mag5_stacks(varydata)
magerr = vari_funcs.magerr5_stacks(varydata)

# Calculate excess variance
excess = vari_funcs.sigmasq(mag, magerr)

# Get X-ray data
varycoord = SkyCoord(varydata['ALPHA_J2000_05B']*u.degree, varydata['DELTA_J2000_05B']*u.degree)

#match with xmm
print('Matching XMM')
xmm = Table.read('UDS_catalogues/XMM_not_Chandra.fits')
xmmcoord = SkyCoord(xmm['RAJ2000'], xmm['DEJ2000'])
idx, d2d , _ = match_coordinates_sky(varycoord, xmmcoord)
mask = d2d<=5*u.arcsec #make sure match is within 5 arcsec (like in topcat)
xmmidx = idx[mask]
xmmdata = xmm[xmmidx] #create table containing xmm details for the variable sources
idx, d2d , _ = match_coordinates_sky(xmmcoord, varycoord)
mask = d2d<=5*u.arcsec #make sure match is within 5 arcsec (like in topcat)
xmmvaryidx = idx[mask]

#match with chandra
print('Matching Chandra')
chan = Table.read('UDS_catalogues/chandra_catalogue.fits')
chan['RA'].unit = u.deg
chan['Dec'].unit = u.deg
chancoord = SkyCoord(chan['RA'], chan['Dec'])
idx, d2d , _ = match_coordinates_sky(varycoord, chancoord)
mask = d2d<=1*u.arcsec #make sure match is within 1 arcsec (like in topcat)
chanidx = idx[mask]
chandata = chan[chanidx] #create table containing chandra details for the variable sources
idx, d2d , _ = match_coordinates_sky(chancoord, varycoord)
mask = d2d<=1*u.arcsec #make sure match is within 5 arcsec (like in topcat)
chanvaryidx = idx[mask]

xrayvaryidx = np.append(chanvaryidx, xmmvaryidx)

# Calculate luminosity distances with both photometric and spectroscopic data
DLspec = cosmo.luminosity_distance(varydata['z_spec'][xrayvaryidx])
DLphot = cosmo.luminosity_distance(varydata['z_m2'][xrayvaryidx])

# Calculate luminosity










