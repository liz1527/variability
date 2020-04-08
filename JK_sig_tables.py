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
### Import data for variables selected in K ###
KKdata = fits.open('variable_tables/K/variables_no06_chi30_DR11data_noxtalkcontam.fits')[1].data
KJdata = fits.open('variable_tables/K/variables_no06_chi30_DR11data_noxtalkcontam_J.fits')[1].data

### Import data for variables selected in J ###
JKdata = fits.open('variable_tables/J/K_extraction/J_variables_chi32_noneg_DR11data_noxtalkcontam_K.fits')[1].data
JJdata = fits.open('variable_tables/J/K_extraction/J_variables_chi32_noneg_DR11data_noxtalkcontam.fits')[1].data

### Import sig data ###
Jsigtb = Table.read('sigma_tables/quad_epoch_sigma_table_J_extra_clean_2arcsec_noneg.fits')
Ksigtb = Table.read('sigma_tables/quad_epoch_sigma_table_K_extra_clean_2arcsec_noneg.fits')

#Jxraydata = Jdata[Jdata['X-ray']==True]

#### Limit to Chandra region for simplicity ###
#KKdata = vari_funcs.field_funcs.chandra_only(KKdata)
#KJdata = vari_funcs.field_funcs.chandra_only(KJdata)
#JKdata = vari_funcs.field_funcs.chandra_only(JKdata)
#JJdata = vari_funcs.field_funcs.chandra_only(JJdata)

### Extract magnitude table and error tables ###
KKflux, KKfluxerr, KKdata = vari_funcs.k_mag_flux.create_quad_error_array(Ksigtb, KKdata, aper=4)
KJflux, KJfluxerr, KJdata = vari_funcs.j_mag_flux.create_quad_error_array_J(Jsigtb, KJdata, aper=4)
JKflux, JKfluxerr, JKdata = vari_funcs.k_mag_flux.create_quad_error_array(Ksigtb, JKdata, aper=4)
JJflux, JJfluxerr, JJdata = vari_funcs.j_mag_flux.create_quad_error_array_J(Jsigtb, JJdata, aper=4)

### Normalise ###
KKfluxnorm, KKfluxerrnorm = vari_funcs.flux_funcs.normalise_flux_and_errors(KKflux, KKfluxerr)
KJfluxnorm, KJfluxerrnorm = vari_funcs.flux_funcs.normalise_flux_and_errors(KJflux, KJfluxerr)
JKfluxnorm, JKfluxerrnorm = vari_funcs.flux_funcs.normalise_flux_and_errors(JKflux, JKfluxerr)
JJfluxnorm, JJfluxerrnorm = vari_funcs.flux_funcs.normalise_flux_and_errors(JJflux, JJfluxerr)

#%% All points
posvar = np.linspace(0,2,5000)
KKmeanflux = np.nanmean(KKfluxnorm, axis=1)
KJmeanflux = np.nanmean(KJfluxnorm, axis=1)
JKmeanflux = np.nanmean(JKfluxnorm, axis=1)
JJmeanflux = np.nanmean(JJfluxnorm, axis=1)

### Set up arrays for K selected ###
allKKout = np.array([])
allKKouterr = np.array([])
allKJout = np.array([])
allKJouterr = np.array([])

### Find z for full tbdata ###
Kz = vari_funcs.get_z(KKdata)

Kbad = [] # array for IDs that do not have a match in the other bad
KbadIDs = [] # array for IDs that do not have a match in the other bad
for n in range(len(KKdata)): #loop over the selection band
    obnum = KKdata['ID'][n] #get DR11 number
    KJmask = np.isin(KJdata['ID'], obnum) #find equivilant J
    if ~np.any(KJmask):
        Kbad.append(n)
        KbadIDs.append(obnum)
        continue
    
    ### Get maximum likelihoods in J and K for that object ###
    KKout = vari_funcs.vary_stats.maximum_likelihood(KKfluxnorm[n,:], 
                                                     KKfluxerrnorm[n,:], 
                                                     KKmeanflux[n], posvar)
    
    KJout = vari_funcs.vary_stats.maximum_likelihood(KJfluxnorm[KJmask,:].reshape(8), 
                                                     KJfluxerrnorm[KJmask,:].reshape(8), 
                                                     KJmeanflux[KJmask], posvar)
    

    ### save output into x-ray and band specfic arrays ###
    allKKout = np.append(allKKout, KKout[0])
    allKJout = np.append(allKJout, KJout[0])
    allKKouterr = np.append(allKKouterr, KKout[1])
    allKJouterr = np.append(allKJouterr, KJout[1])
        
### remove rows with no J match ###
KKtable = Table(np.copy(KKdata))
KKtable.remove_rows(Kbad)
rowmask = ~np.isin(KKdata['ID'], KbadIDs)
KKflux = KKflux[rowmask]
KJflux = KJflux
KKfluxerr = KKfluxerr[rowmask]
KJfluxerr = KJfluxerr
Kz = Kz[rowmask]
K_xray = KKtable['X-ray']
               
### Set up arrays for J selected ###
allJKout = np.array([])
allJKouterr = np.array([])
allJJout = np.array([])
allJJouterr = np.array([])

### Find z for full tbdata ###
Jz = vari_funcs.get_z(JJdata)

Jbad = [] # array for IDs that do not have a match in the other bad
JbadIDs = [] # array for IDs that do not have a match in the other bad
for n in range(len(JJdata)): #loop over the selection band
    obnum = JJdata['ID'][n] #find DR11 number
    JKmask = np.isin(JKdata['ID'], obnum) #find that object in K array
    if ~np.any(JKmask):
        Jbad.append(n)
        JbadIDs.append(obnum)
        continue
    
    ### Get maximum likelihoods in J and K for that object ###
    JJout = vari_funcs.vary_stats.maximum_likelihood(JJfluxnorm[n,:], 
                                                     JJfluxerrnorm[n,:], 
                                                     JJmeanflux[n], posvar)
    
    JKout = vari_funcs.vary_stats.maximum_likelihood(JKfluxnorm[JKmask,:].reshape(7), 
                                                     JKfluxerrnorm[JKmask,:].reshape(7), 
                                                     JKmeanflux[JKmask], posvar)


    ### save output into x-ray and band specfic arrays ###
    allJKout = np.append(allJKout, JKout[0])
    allJJout = np.append(allJJout, JJout[0])
    allJKouterr = np.append(allJKouterr, JKout[1])
    allJJouterr = np.append(allJJouterr, JJout[1])

#### remove rows with no K match ###
JJtable = Table(np.copy(JJdata))
J_xray = JJtable['X-ray']
#JJtable.remove_rows(Jbad)
#rowmask = ~np.isin(JJdata['ID'], JbadIDs)
#JKflux = JKflux[rowmask]
#JJflux = JJflux[rowmask]
#JKfluxerr = JKfluxerr[rowmask]
#JJfluxerr = JJfluxerr[rowmask]
#Jz = Jz[rowmask]

### Sort out columns wanted for tables ###
def get_cols(tbdata):
    ID = tbdata['ID']
    ra = tbdata['ALPHA_J2000']
    dec = tbdata['DELTA_J2000']
    return ID, ra, dec

KID, Kra, Kdec = get_cols(KKtable)
JID, Jra, Jdec = get_cols(JJtable)

### Make ID, ra, dec tables for K and J separately ###
Ktable = Table([KID, Kra, Kdec],names=['ID','RA','Dec'])
Jtable = Table([JID, Jra, Jdec],names=['ID','RA','Dec'])

### add sig and err columns for J and K to both tables ###
def add_col(arr, tbdata, name):
    ### create column ###
    new_col = Column(arr, name=name)
    ### add to tbdata ###
    tbdata.add_column(new_col)
    return tbdata
Ktable = add_col(Kz, Ktable, 'z')
Ktable = add_col(KKflux, Ktable, 'Flux_K')
Ktable = add_col(KKfluxerr, Ktable, 'Fluxerr_K')
Ktable = add_col(KJflux, Ktable, 'Flux_J')
Ktable = add_col(KJfluxerr, Ktable, 'Fluxerr_J')
Ktable = add_col(allKKout, Ktable, 'sig_K')
Ktable = add_col(allKKouterr, Ktable, 'sig_K_err')
Ktable = add_col(allKJout, Ktable, 'sig_J')
Ktable = add_col(allKJouterr, Ktable, 'sig_J_err')
Ktable = add_col(K_xray, Ktable, 'X-ray')
Jtable = add_col(Jz, Jtable, 'z')
Jtable = add_col(JKflux, Jtable, 'Flux_K')
Jtable = add_col(JKfluxerr, Jtable, 'Fluxerr_K')
Jtable = add_col(JJflux, Jtable, 'Flux_J')
Jtable = add_col(JJfluxerr, Jtable, 'Fluxerr_J')
Jtable = add_col(allJKout, Jtable, 'sig_K')
Jtable = add_col(allJKouterr, Jtable, 'sig_K_err')
Jtable = add_col(allJJout, Jtable, 'sig_J')
Jtable = add_col(allJJouterr, Jtable, 'sig_J_err')
Jtable = add_col(J_xray, Jtable, 'X-ray')


### Match ID columns so that can make one megatable ###
mask = np.isin(Ktable['ID'], Jtable['ID'])
bothtable = Ktable[mask]
Jtable = Jtable[~np.isin(Jtable['ID'], Ktable['ID'])]
Ktable = Ktable[~mask]

finaltable = vstack([bothtable, Ktable, Jtable])

### Calculate and add a chi_sq values for J and K ###
Kchi = vari_funcs.vary_stats.my_chisquare_err(finaltable['Flux_K'], finaltable['Fluxerr_K'])
Jchi = vari_funcs.vary_stats.my_chisquare_err(finaltable['Flux_J'], finaltable['Fluxerr_J'])
finaltable = add_col(Kchi, finaltable, 'Chi_K')
finaltable = add_col(Jchi, finaltable, 'Chi_J')

#finaltable.write('variable_tables/NIR_variables_J_and_K.fits', overwrite=True)
finaltable.write('variable_tables/J_and_K_variables_varystats.fits', overwrite=True)

end = time.time()
print(end-start)












