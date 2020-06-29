#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  1 13:30:42 2018

Code to create lightcurve tables for our summer student to work on

@author: ppxee
"""
import time
start = time.time()
print(start)

### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
from astropy.table import Table, Column, vstack #for handling tables
import numpy as np #for handling arrays
#import math
#from astropy.stats import median_absolute_deviation
import vari_funcs #my module to help run code neatly
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
plt.close('all') #close any open plots
#from numpy.lib.recfunctions import append_fields

def add_col(arr, tbdata, name):
    ### create column ###
    new_col = Column(arr, name=name)
    ### add to tbdata ###
    tbdata.add_column(new_col)
    return tbdata

#%% Open the fits files and get data ###
### Import sem data for all ###
full = Table.read('UDS_catalogues/full_varystats_neg.fits')

### Import month data (doing separately as need fits format) ###
Kdata = fits.open('mag_flux_tables/K/month/month_flux_table_best_K_extra_quad_clean_38_DR11_ID.fits')[1].data
Jdata = fits.open('mag_flux_tables/J/month/month_flux_table_best_J_extra_quad_clean_DR11_ID.fits')[1].data

### Import sig data ###
Ksigtb = Table.read('sigma_tables/month_quad_epoch_sigma_table_K_extra_quad_clean_38_2arcsec_neg.fits')
Jsigtb = Table.read('sigma_tables/month_quad_epoch_sigma_table_J_extra_quad_clean_2arcsec_neg.fits')


#%% Extract magnitude table and error tables ###
Kfluxmon, Kfluxerrmon, Kdata = vari_funcs.k_mag_flux.create_quad_error_array_month(Ksigtb, Kdata, aper=4)
Jfluxmon, Jfluxerrmon, Jdata = vari_funcs.j_mag_flux.create_quad_error_array_month_J(Jsigtb, Jdata, aper=4)

#%% Match ID Columns ###
### Get matching IDs in J and K to full first first ###
kmask1 = np.isin(Kdata['ID'], full['ID'])
jmask1 = np.isin(Jdata['ID'], full['ID'])

### Now match J and K tables ###
kmask2 = np.isin(Kdata['ID'], Jdata['ID'])
jmask2 = np.isin(Jdata['ID'], Kdata['ID'])

### Join and apply masks ###
kmask = kmask1*kmask2.astype(bool)
Kdata = Kdata[kmask]
Kfluxmon = Kfluxmon[kmask]
Kfluxerrmon = Kfluxerrmon[kmask]
jmask = jmask1*jmask2.astype(bool)
Jdata = Jdata[jmask]
Jfluxmon = Jfluxmon[jmask]
Jfluxerrmon = Jfluxerrmon[jmask]

### Match full to final J and K length ###
mask = np.isin(full['ID'], Jdata['ID'])
full = full[mask]
#%% Set up base table ###
basetable = Table()

### add ID, RA and Dec columns ###
basetable = add_col(full['ID'], basetable, 'DR11_ID')
basetable = add_col(full['RA'], basetable, 'RA')
basetable = add_col(full['Dec'], basetable, 'Dec')

#%% Set up semester table ###
semtable = Table(np.copy(basetable))

### set up semester arrays ###
ksemesters = ['05B', '07B', '08B', '09B', '10B', '11B', '12B']
jsemesters = ['05B', '06B', '07B', '08B', '09B', '10B', '11B', '12B']

### iterate over semesters to separate flux columns ###
for n, sem in enumerate(ksemesters):
    ksemflux = full['Flux_K'][:,n]
    ksemfluxerr = full['Fluxerr_K'][:,n]
    semtable = add_col(ksemflux, semtable, 'K_flux_'+sem)
    semtable = add_col(ksemfluxerr, semtable, 'K_fluxerr_'+sem)


for n, sem in enumerate(jsemesters):
    jsemflux = full['Flux_J'][:,n]
    jsemfluxerr = full['Fluxerr_J'][:,n]
    semtable = add_col(jsemflux, semtable, 'J_flux_'+sem)
    semtable = add_col(jsemfluxerr, semtable, 'J_fluxerr_'+sem)

semtable.write('UDS_catalogues/semester_lightcurves_J_and_K.fits', overwrite=True)


#%% Set up month table ###
montable = Table(np.copy(basetable))

### Set up month arrays ###
kmonths = ['sep05','oct05','nov05','dec05', 'jan06', 'jan07', 'aug07', 'sep07', 
           'oct07', 'sep08', 'oct08', 'nov08', 'jul09', 'aug09', 'sep09', 
           'oct09', 'nov09', 'dec09', 'jan10', 'feb10', 'aug10', 'sep10', 
           'oct10', 'nov10', 'dec10', 'jan11', 'aug11', 'sep11', 'oct11', 
           'nov11', 'dec11', 'jan12', 'feb12', 'jul12', 'aug12', 'sep12', 
           'oct12', 'nov12']

jmonths = ['sep05', 'oct05', 'nov05', 'dec05', 'jan06', 'oct06', 'nov06',
          'dec06', 'aug07', 'sep07', 'oct07', 'oct08', 'nov08', 'aug09',
          'sep09', 'oct09', 'nov09', 'dec09', 'aug10', 'sep10', 'oct10',
          'nov10', 'dec10', 'jan11', 'aug11', 'sep11', 'oct11', 'nov11',
          'dec11', 'jan12', 'jul12', 'aug12', 'sep12', 'oct12', 'nov12']

### iterate over months to separate flux columns ###
for n, mon in enumerate(kmonths):
    kmonflux = Kfluxmon[:,n]
    kmonfluxerr = Kfluxerrmon[:,n]
    montable = add_col(kmonflux, montable, 'K_flux_'+mon)
    montable = add_col(kmonfluxerr, montable, 'K_fluxerr_'+mon)


for n, mon in enumerate(jmonths):
    jmonflux = Jfluxmon[:,n]
    jmonfluxerr = Jfluxerrmon[:,n]
    montable = add_col(jmonflux, montable, 'J_flux_'+mon)
    montable = add_col(jmonfluxerr, montable, 'J_fluxerr_'+mon)
    
montable.write('UDS_catalogues/month_lightcurves_J_and_K.fits', overwrite=True)

##finaltable.write('variable_tables/NIR_variables_J_and_K.fits', overwrite=True)
#finaltable.write('variable_tables/J_and_K_variables_month_varystats_DR11data.fits', 
#                 overwrite=True)

end = time.time()
print(end-start)












