#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 12:36:55 2020

add x-ray boolean to UDS cat

@author: ppxee
"""

from astropy.table import Table, Column
import numpy as np

uds = Table.read('UDS_catalogues/DR11-2arcsec-Jun-30-2019.fits')
xray = Table.read('mag_flux_tables/K/xray_mag_flux_table_best_K_extra_clean_DR11data.fits')

mask = np.isin(uds['ID'], xray['ID'])

xraycol = np.zeros(len(uds))
xraycol[mask] = 1
xraycol = xraycol.astype('bool')

new_col = Column(xraycol, name='X-ray')
uds.add_column(new_col)

uds.write('UDS_catalogues/DR11-2arcsec-Jun-30-2019_database.fits')