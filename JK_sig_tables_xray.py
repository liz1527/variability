#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  1 13:30:42 2018

Code to compute maximum likelihood variability and save table of joint J and K 
detected variables with variability amplitudes in both bands

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

#%% Open the fits files and get data ###
### Import data for xray sources ###
Kdata = fits.open('mag_flux_tables/K/mag_flux_table_best_K_extra_clean_DR11data.fits')[1].data
Jdata = fits.open('mag_flux_tables/J/mag_flux_table_best_J_extra_clean_DR11data.fits')[1].data

### Import sig data ###
Jsigtb = Table.read('sigma_tables/quad_epoch_sigma_table_J_extra_clean_2arcsec_noneg.fits')
Ksigtb = Table.read('sigma_tables/quad_epoch_sigma_table_K_extra_clean_2arcsec_noneg.fits')

#Jxraydata = Jdata[Jdata['X-ray']==True]

#### Limit to Chandra region for simplicity ###
#KKdata = vari_funcs.field_funcs.chandra_only(KKdata)
#KJdata = vari_funcs.field_funcs.chandra_only(KJdata)
#JKdata = vari_funcs.field_funcs.chandra_only(JKdata)
#JJdata = vari_funcs.field_funcs.chandra_only(JJdata)

### Remove negatives ###
Kflux = vari_funcs.k_mag_flux.flux_stacks(Kdata, aper=4)
Jflux = vari_funcs.j_mag_flux.flux_stacks(Jdata, aper=4)
Kflux, Kdata = vari_funcs.flux_funcs.noneg(Kflux, Kdata)
Jflux, Jdata = vari_funcs.flux_funcs.noneg(Jflux, Jdata)

### Make sure tables match ###
Jmask = np.isin(Jdata['ID'], Kdata['ID'])
Kmask = np.isin(Kdata['ID'], Jdata['ID'])
Jdata = Jdata[Jmask]
Kdata = Kdata[Kmask]

### Extract magnitude table and error tables ###
Kflux, Kfluxerr, Kdata = vari_funcs.k_mag_flux.create_quad_error_array(Ksigtb, Kdata, aper=4)
Jflux, Jfluxerr, Jdata = vari_funcs.j_mag_flux.create_quad_error_array_J(Jsigtb, Jdata, aper=4)

### Make sure tables still match ###
Kmask = np.isin(Kdata['ID'], Jdata['ID'])
Jmask = np.isin(Jdata['ID'], Kdata['ID'])
Kdata = Kdata[Kmask]
Jdata = Jdata[Jmask]
Kflux = Kflux[Kmask]
Jflux = Jflux[Jmask]
Kfluxerr = Kfluxerr[Kmask]
Jfluxerr = Jfluxerr[Jmask]

### Normalise ###
Kfluxnorm, Kfluxerrnorm = vari_funcs.flux_funcs.normalise_flux_and_errors(Kflux, Kfluxerr)
Jfluxnorm, Jfluxerrnorm = vari_funcs.flux_funcs.normalise_flux_and_errors(Jflux, Jfluxerr)

### Find z for full tbdata ###
z = vari_funcs.get_z(Kdata)
#%% All points
posvar = np.linspace(0,2,5000)
Kmeanflux = np.nanmean(Kfluxnorm, axis=1)
Jmeanflux = np.nanmean(Jfluxnorm, axis=1)

### Set up arrays for K selected ###
allKout = np.array([])
allKouterr = np.array([])
allJout = np.array([])
allJouterr = np.array([])

Kbad = [] # array for IDs that do not have a match in the other bad
KbadIDs = [] # array for IDs that do not have a match in the other bad
for n in range(len(Kdata)): #loop over the selection band
    obnum = Kdata['ID'][n] #get DR11 number
    print(obnum)
    Jmask = np.isin(Jdata['ID'], obnum) #find equivilant J
    if ~np.any(Jmask):
        Kbad.append(n)
        KbadIDs.append(obnum)
        continue
    
    ### Get maximum likelihoods in J and K for that object ###
    Kout = vari_funcs.vary_stats.maximum_likelihood(Kfluxnorm[n,:], 
                                                     Kfluxerrnorm[n,:], 
                                                     Kmeanflux[n], posvar)
    
    Jout = vari_funcs.vary_stats.maximum_likelihood(Jfluxnorm[Jmask,:].reshape(8), 
                                                     Jfluxerrnorm[Jmask,:].reshape(8), 
                                                     Jmeanflux[Jmask], posvar)
    

    ### save output into x-ray and band specfic arrays ###
    allKout = np.append(allKout, Kout[0])
    allJout = np.append(allJout, Jout[0])
    allKouterr = np.append(allKouterr, Kout[1])
    allJouterr = np.append(allJouterr, Jout[1])
        
#### remove rows with no J match ###
#Ktable = Table(np.copy(Kdata))
#Ktable.remove_rows(Kbad)
#Krowmask = ~np.isin(Kdata['ID'], KbadIDs)
#Jrowmask = ~np.isin(Jdata['ID'], KbadIDs)
#Kflux = Kflux[Krowmask]
#Jflux = Jflux[Jrowmask]
#Kfluxerr = Kfluxerr[Krowmask]
#Jfluxerr = Jfluxerr[Jrowmask]
#z = z[Krowmask]

### Sort out columns wanted for tables ###
def get_cols(tbdata):
    ID = tbdata['ID']
    ra = tbdata['ALPHA_J2000']
    dec = tbdata['DELTA_J2000']
    return ID, ra, dec

Ktable = Table(np.copy(Kdata))
KID, Kra, Kdec = get_cols(Ktable)

### Make ID, ra, dec tables for K and J separately ###
Ktable = Table([KID, Kra, Kdec],names=['ID','RA','Dec'])

### add sig and err columns for J and K to both tables ###
def add_col(arr, tbdata, name):
    ### create column ###
    new_col = Column(arr, name=name)
    ### add to tbdata ###
    tbdata.add_column(new_col)
    return tbdata
Ktable = add_col(z, Ktable, 'z')
Ktable = add_col(Kflux, Ktable, 'Flux_K')
Ktable = add_col(Kfluxerr, Ktable, 'Fluxerr_K')
Ktable = add_col(Jflux, Ktable, 'Flux_J')
Ktable = add_col(Jfluxerr, Ktable, 'Fluxerr_J')
Ktable = add_col(allKout, Ktable, 'sig_K')
Ktable = add_col(allKouterr, Ktable, 'sig_K_err')
Ktable = add_col(allJout, Ktable, 'sig_J')
Ktable = add_col(allJouterr, Ktable, 'sig_J_err')

#
#### Match ID columns so that can make one megatable ###
#mask = np.isin(Ktable['ID'], Jtable['ID'])
#bothtable = Ktable[mask]
#Jtable = Jtable[~np.isin(Jtable['ID'], Ktable['ID'])]
#Ktable = Ktable[~mask]

#finaltable = vstack([bothtable, Ktable, Jtable])

### Calculate and add a chi_sq values for J and K ###
Kchi = vari_funcs.vary_stats.my_chisquare_err(Ktable['Flux_K'], Ktable['Fluxerr_K'])
Jchi = vari_funcs.vary_stats.my_chisquare_err(Ktable['Flux_J'], Ktable['Fluxerr_J'])
Ktable = add_col(Kchi, Ktable, 'Chi_K')
Ktable = add_col(Jchi, Ktable, 'Chi_J')

#finaltable.write('variable_tables/NIR_variables_J_and_K.fits', overwrite=True)
Ktable.write('UDS_catalogues/full_varystats_neg.fits', overwrite=True)

end = time.time()
print(end-start)












