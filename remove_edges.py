#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  5 14:55:10 2020

@author: ppxee
"""
from astropy.io import fits #for handling fits
from astropy.table import Table
import vari_funcs

tbdata = Table.read('mag_flux_tables/K/month/month_mag_flux_table_best_K_extra_quad_clean.fits')#[1].data

tbdata = vari_funcs.field_funcs.remove_edges(tbdata, 'sep05')

tbdata.write('mag_flux_tables/K/month/month_mag_flux_table_best_K_extra_quad_clean_edged.fits')