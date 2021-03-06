#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  1 13:30:42 2018

Code to compute maximum likelihood variability and save table of joint J and K 
detected variables with month lightcures and variability amplitudes in the K 
bands - should be expanded to do J band once these lcs have been processed

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
### Import data for variables selected in J and K Yearstacks ###
varys = Table.read('variable_tables/J_and_K_variables_varystats_DR11data.fits')

### Import month data (doing separately as need fits format) ###
Kdata = fits.open('variable_tables/J_and_K_variables_varystats_DR11data_monthdata.fits')[1].data
Jdata = fits.open('variable_tables/J_and_K_variables_varystats_DR11data_J_monthdata.fits')[1].data

### Import sig data ###
Ksigtb = Table.read('sigma_tables/month_quad_epoch_sigma_table_K_extra_quad_clean_38_2arcsec_neg.fits')
Jsigtb = Table.read('sigma_tables/month_quad_epoch_sigma_table_J_extra_quad_clean_2arcsec_neg.fits')

#Jxraydata = Jdata[Jdata['X-ray']==True]

#### Limit to Chandra region for simplicity ###
#KKdata = vari_funcs.field_funcs.chandra_only(KKdata)
#KJdata = vari_funcs.field_funcs.chandra_only(KJdata)
#JKdata = vari_funcs.field_funcs.chandra_only(JKdata)
#JJdata = vari_funcs.field_funcs.chandra_only(JJdata)

### Extract magnitude table and error tables ###
Kflux, Kfluxerr, Kdata = vari_funcs.k_mag_flux.create_quad_error_array_month(Ksigtb, Kdata, aper=4)
Jflux, Jfluxerr, Jdata = vari_funcs.j_mag_flux.create_quad_error_array_month_J(Jsigtb, Jdata, aper=4)
#Kflux, Kfluxerr, Kdata = vari_funcs.flux_funcs.nanneg(Kflux, Kfluxerr, Kdata)

### Normalise ###
Kfluxnorm, Kfluxerrnorm = vari_funcs.flux_funcs.normalise_flux_and_errors(Kflux, Kfluxerr)
Jfluxnorm, Jfluxerrnorm = vari_funcs.flux_funcs.normalise_flux_and_errors(Jflux, Jfluxerr)

#%% All points
posvar = np.linspace(0,2,5000)
Kmeanflux = np.nanmean(Kfluxnorm, axis=1)
Jmeanflux = np.nanmean(Jfluxnorm, axis=1)

### Set up arrays for K selected ###
allKout = np.array([])
allKouterr = np.array([])

### Find z for full tbdata ###
Kz = vari_funcs.get_z(Kdata)

Kbad = [] # array for IDs that do not have a match in the other bad
KbadIDs = [] # array for IDs that do not have a match in the other bad
for n in range(len(Kdata)): #loop over the selection band
    obnum = Kdata['ID'][n] #get DR11 number
    print(obnum) # give an idea of progress
    
    ### Get maximum likelihoods in J and K for that object ###
    Kout = vari_funcs.vary_stats.maximum_likelihood(Kfluxnorm[n,:], 
                                                     Kfluxerrnorm[n,:], 
                                                     Kmeanflux[n], posvar)

    ### save output into x-ray and band specfic arrays ###
    allKout = np.append(allKout, Kout[0])
    allKouterr = np.append(allKouterr, Kout[1])
        
               
## Set up arrays for J selected ###
allJout = np.array([])
allJouterr = np.array([])

### Find z for full tbdata ###
Jz = vari_funcs.get_z(Jdata)

Jbad = [] # array for IDs that do not have a match in the other bad
JbadIDs = [] # array for IDs that do not have a match in the other bad
for n in range(len(Jdata)): #loop over the selection band
    obnum = Jdata['ID'][n] #find DR11 number
    
    ### Get maximum likelihoods in J and K for that object ###
    Jout = vari_funcs.vary_stats.maximum_likelihood(Jfluxnorm[n,:], 
                                                    Jfluxerrnorm[n,:], 
                                                    Jmeanflux[n], posvar)
    

    ### save output into x-ray and band specfic arrays ###
    allJout = np.append(allJout, Jout[0])
    allJouterr = np.append(allJouterr, Jout[1])

##### remove rows with no K match ###
#Jtable = Table(np.copy(Jdata))
#J_xray = Jtable['X-ray']
##JJtable.remove_rows(Jbad)
##rowmask = ~np.isin(JJdata['ID'], JbadIDs)
##JKflux = JKflux[rowmask]
##JJflux = JJflux[rowmask]
##JKfluxerr = JKfluxerr[rowmask]
##JJfluxerr = JJfluxerr[rowmask]
##Jz = Jz[rowmask]

#%% Set up table with original varys as base ###
### Check all variables made it to the end ###
#mask = np.isin(varys['ID'], Kdata['ID'])
#Ktable = varys[mask]
Ktable = Table(Kdata)
Ktable['X-ray'] = Ktable['X-ray'].astype('bool')

### remove monthdata columns ###
Ktable.keep_columns(varys.colnames)

### add sig and err columns for J and K to both tables ###
def add_col(arr, tbdata, name):
    ### create column ###
    new_col = Column(arr, name=name)
    ### add to tbdata ###
    tbdata.add_column(new_col)
    return tbdata
Ktable = add_col(Kflux, Ktable, 'Month_Flux_K')
Ktable = add_col(Kfluxerr, Ktable, 'Month_Fluxerr_K')
Ktable = add_col(allKout, Ktable, 'Month_sig_K')
Ktable = add_col(allKouterr, Ktable, 'Month_sig_K_err')
Ktable = add_col(Jflux, Ktable, 'Month_Flux_J')
Ktable = add_col(Jfluxerr, Ktable, 'Month_Fluxerr_J')
Ktable = add_col(allJout, Ktable, 'Month_sig_J')
Ktable = add_col(allJouterr, Ktable, 'Month_sig_J_err')


### Match ID columns so that can make one megatable ###
#mask = np.isin(Ktable['ID'], Jtable['ID'])
#bothtable = Ktable[mask]
#Jtable = Jtable[~np.isin(Jtable['ID'], Ktable['ID'])]
#Ktable = Ktable[~mask]

finaltable = Ktable#vstack([bothtable, Ktable, Jtable])

### Calculate and add a chi_sq values for J and K ###
Kchi = vari_funcs.vary_stats.my_chisquare_err(finaltable['Month_Flux_K'], finaltable['Month_Fluxerr_K'])
Jchi = vari_funcs.vary_stats.my_chisquare_err(finaltable['Month_Flux_J'], finaltable['Month_Fluxerr_J'])
finaltable = add_col(np.array(Kchi), finaltable, 'Month_Chi_K')
finaltable = add_col(np.array(Jchi), finaltable, 'Month_Chi_J')

#finaltable.write('variable_tables/NIR_variables_J_and_K.fits', overwrite=True)
finaltable.write('variable_tables/J_and_K_variables_month_varystats_DR11data.fits', 
                 overwrite=True)

end = time.time()
print(end-start)












