#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 12 10:48:02 2018

code to create star stamp table to detemine psf

@author: ppxee
"""

### Import Modules Required ###
from astropy.coordinates import match_coordinates_sky
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table
import sys

sem = str(sys.argv[1]) #1st arguement from call is the semester

### read in fits_LDAC table and stars table ###
data = Table.read(sem+'_stamp_output.fits',hdu=2)
stars = Table.read('UDS_catalogues/DR8-secure-stars.fits')

### Define coordinates ###
stampcoord = SkyCoord(data['ALPHA_J2000'], data['DELTA_J2000'])
starscoord = SkyCoord(stars['RA']*u.degree, stars['DEC']*u.degree)

### Match catalogues and create new table ###
idx, d2d , _ = match_coordinates_sky(starscoord, stampcoord)
starstamps = data[idx]

starstamps.write(sem+'_star_stamps_table.fits')














