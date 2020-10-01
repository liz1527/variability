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
### Import data for J and K Yearstacks ###
varys = Table.read('UDS_catalogues/chandra_varystats_noneg.fits')

### Import month data ###
monthdata = Table.read('UDS_catalogues/month_lightcurves_J_and_K.fits')

### Make sure IDs match up ###
mask = np.isin(monthdata['DR11_ID'], varys['ID'])
monthdata = monthdata[mask]
mask = np.isin(varys['ID'], monthdata['DR11_ID'])
varys = varys[mask]

#### Limit to Chandra region for simplicity ###
#KKdata = vari_funcs.field_funcs.chandra_only(KKdata)
#KJdata = vari_funcs.field_funcs.chandra_only(KJdata)
#JKdata = vari_funcs.field_funcs.chandra_only(JKdata)
#JJdata = vari_funcs.field_funcs.chandra_only(JJdata)

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
Kflux = np.empty([len(monthdata), len(kmonths)])
Kfluxerr = np.empty([len(monthdata), len(kmonths)])
for n, mon in enumerate(kmonths):
    Kflux[:,n] = monthdata['K_flux_'+mon]
    Kfluxerr[:,n] = monthdata['K_fluxerr_'+mon]


Jflux = np.empty([len(monthdata), len(jmonths)])
Jfluxerr = np.empty([len(monthdata), len(jmonths)])
for n, mon in enumerate(jmonths):
    Jflux[:,n] = monthdata['J_flux_'+mon]
    Jfluxerr[:,n] = monthdata['J_fluxerr_'+mon]

### Normalise ###
Kfluxnorm, Kfluxerrnorm = vari_funcs.flux_funcs.normalise_flux_and_errors(Kflux, Kfluxerr)
Jfluxnorm, Jfluxerrnorm = vari_funcs.flux_funcs.normalise_flux_and_errors(Jflux, Jfluxerr)

#%% All points
posvar = np.linspace(0,2,5000)
Kmeanflux = np.nanmean(Kfluxnorm, axis=1)
Jmeanflux = np.nanmean(Jfluxnorm, axis=1)

### Set up arrays for K ###
allKout = np.array([])
allKouterr = np.array([])

## Set up arrays for J ###
allJout = np.array([])
allJouterr = np.array([])

for n in range(len(monthdata)): #loop over the selection band
    obnum = monthdata['DR11_ID'][n] #get DR11 number
    print(obnum) # give an idea of progress
    
    ### Get maximum likelihoods in K for that object ###
    Kout = vari_funcs.vary_stats.maximum_likelihood(Kfluxnorm[n,:], 
                                                     Kfluxerrnorm[n,:], 
                                                     Kmeanflux[n], posvar)

    ### save output into band specfic arrays ###
    allKout = np.append(allKout, Kout[0])
    allKouterr = np.append(allKouterr, Kout[1])
        
    ### Get maximum likelihoods in J for that object ###
    Jout = vari_funcs.vary_stats.maximum_likelihood(Jfluxnorm[n,:], 
                                                    Jfluxerrnorm[n,:], 
                                                    Jmeanflux[n], posvar)
    

    ### save output into  band specfic arrays ###
    allJout = np.append(allJout, Jout[0])
    allJouterr = np.append(allJouterr, Jout[1])


#%% Set up table with original varys as base ###
Ktable = Table(varys)
#Ktable['X-ray'] = Ktable['X-ray'].astype('bool')

#### remove monthdata columns ###
#Ktable.keep_columns(varys.colnames)

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
finaltable.write('UDS_catalogues/chandra_month_varystats_noneg.fits', 
                 overwrite=True)

end = time.time()
print(end-start)












