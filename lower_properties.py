#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 24 15:10:37 2019

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
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
plt.close('all') #close any open plots

#%% Check position on chi-flux plot first ###

### Open the fits files and get data ###
tbdata = fits.open('mag_flux_tables/mag_flux_table_best_extra_clean_no06.fits')[1].data
chandata = fits.open('mag_flux_tables/xray_mag_flux_table_best_extra_clean_no06.fits')[1].data
sdata = fits.open('mag_flux_tables/stars_mag_flux_table_extra_clean_no06.fits')[1].data
sigtb = Table.read('sigma_tables/quad_epoch_sigma_table_extra_clean_no06_2arcsec.fits')

def prep_data(tbdata):
    ### Remove edges ###
    tbdata = vari_funcs.remove_edges(tbdata)
    
    ## Create arrays of flux values ###
    flux = vari_funcs.flux4_stacks(tbdata)
    
    ### remove values that are negative ###
    flux, tbdata = vari_funcs.noneg(flux, tbdata)
    
    ### Get error arrays ###
    flux, fluxerr, tbdata = vari_funcs.create_quad_error_array(sigtb, tbdata, aper=4)
    
    return flux, fluxerr, tbdata

### Prep data ###
flux, fluxerr, tbdata = prep_data(tbdata)
fluxchan, chanerr, chandata = prep_data(chandata)
sflux, serr, sdata = prep_data(sdata)


### reset X-ray column as messed up by stacking ###
tbdata['X-ray'][tbdata['X-ray']==70] = False 
tbdata['X-ray'][tbdata['X-ray']==84] = True

### Check chisq plot looks correct ###
fig,_ = vari_funcs.flux_variability_plot(flux, fluxchan, 'chisq', 
                                       fluxerr=fluxerr, chanerr=chanerr,
                                       starflux=sflux, starfluxerr=serr,
                                       #normalised=True, 
                                       stars=True, scale='log')
fig.canvas.mpl_connect('pick_event', vari_funcs.onpickflux_2arcsec)

#varydata = fits.open('variable_tables/no06_variables_chi30_2arcsec_DR11data_restframe.fits')[1].data
#varydatalow = vari_funcs.flux_split(varydata, 'lower')
varydatalow = fits.open('variable_tables/no06_variables_chi30_2arcsec_spec_DR11.fits')[1].data
#varydatalow = vari_funcs.flux_split(varydata, 'lower')

varyfluxlow, varyfluxerrlow, varydatalow = prep_data(varydatalow)
varymeanlow = np.nanmean(varyfluxlow, axis=1)
varychilow = vari_funcs.my_chisquare_err(varyfluxlow, varyfluxerrlow)
plt.plot(varymeanlow, varychilow, 'kd', mfc='None')

