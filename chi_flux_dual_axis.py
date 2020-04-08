#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 11:35:11 2018

Code to creat chi-flux plot using the bootstraped error bars and 2 arcsec K-band

This version has dual x-axis showing magnitudes and fluxes, and is the version
used in my first paper.

@author: ppxee
"""

### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
from astropy.table import Table #for handling tables
import numpy as np #for handling arrays
#import math
#from astropy.stats import median_absolute_deviation
from scipy import stats
import vari_funcs #my module to help run code neatly
#from scipy.stats import chisquare
plt.close('all') #close any open plots
font = {'family' : 'DejaVu Sans',
        'weight' : 'normal',
        'size'   : 16}

plt.rc('font', **font)

### Open the fits files and get data ###
tbdata = fits.open('mag_flux_tables/K/mag_flux_table_best_K_extra_clean.fits')[1].data
chandata = fits.open('mag_flux_tables/K/xray_mag_flux_table_best_K_extra_clean.fits')[1].data
sdata = fits.open('mag_flux_tables/K/stars_mag_flux_table_K_extra_clean.fits')[1].data
#xtalkdata = fits.open('mag_flux_tables/K/mag_flux_table_xtalk_K_extra_clean.fits')[1].data
xtalkdata = fits.open('variable_tables/K/variables_no06_chi30_neg_DR11data_xtalkcontam.fits')[1].data
sigtb = Table.read('sigma_tables/quad_epoch_sigma_table_K_extra_clean_2arcsec_neg.fits')

def prep_data(tbdata):
    ### Remove edges ###
    tbdata = vari_funcs.field_funcs.remove_edges(tbdata)
    
    ### Create arrays of flux values ###
    flux = vari_funcs.k_mag_flux.flux4_stacks(tbdata)
    
    ### remove values that are negative ###
#    flux, tbdata = vari_funcs.flux_funcs.noneg(flux, tbdata)
    
    ### Get error arrays ###
    flux, fluxerr, tbdata = vari_funcs.k_mag_flux.create_quad_error_array(sigtb, tbdata, aper=4)
    
    
    ### get average flux for each object ###
    avgfluxperob = np.nanmean(flux, axis=1)
    
    ### find chi squared ###
    vary = vari_funcs.vary_stats.my_chisquare_err(flux, fluxerr)
    
    return tbdata, avgfluxperob, vary

tbdata, avgfluxperob, vary = prep_data(tbdata)
chandata, avgfluxchanperob, varychan = prep_data(chandata)
sdata, savgfluxperob, varystar = prep_data(sdata)
xtalkdata, xtalkavgfluxperob, varyxtalk = prep_data(xtalkdata)

### reset X-ray column as messed up by stacking ###
tbdata['X-ray'][tbdata['X-ray']==70] = False 
tbdata['X-ray'][tbdata['X-ray']==84] = True

### Create plot with two axes ##
fig = plt.figure(figsize=[8,8])
ax1 = fig.add_subplot(111) # Subplot covering whole plot
ax2 = ax1.twiny() # Twin of the first subplot

### Plot the chi sq v mean ###
plt.plot(savgfluxperob, varystar, 'm*', mfc = 'none', markersize = 10,
             label='DR11 Star') 
line, = plt.plot(avgfluxperob, vary, 'b+', label='Galaxy', picker=2)
plt.plot(avgfluxchanperob, varychan, 'ro', mfc = 'none', markersize = 10,
         label='X-ray detected') #no picker as will be selected in the UDS point
plt.plot(xtalkavgfluxperob, varyxtalk, 'kd', mfc = 'none', markersize = 10,
         label='X-talk') #no picker as will be selected in the UDS point

    
### Apply required plot charateristics ###
ax1.set_xscale('log')
ax1.set_yscale('log')
    
ax1.set_ylim(3e-2,3e4)
#ax1.set_xlim(8e1, 1e7) #if flux noneg
ax1.set_xlim(4e0, 1e7) #if flux neg
#    plt.xlim(13,26) #if mag

ax1.set_xlabel('Mean Flux')
ax1.set_ylabel(r'$\chi^{2}$')

ticks = np.array([1e2, 1e3, 1e4, 1e5, 1e6, 1e7]) #set flux ticks
new_ticks = 30 - 2.5*np.log10(ticks) #find equivilant mag ticks

ax2.set_xlim(ax1.get_xlim()) #set twin to same limits
ax2.set_xscale('log') #set twin to same scale
ax2.set_xticks(ticks) #set flux ticks
ax2.set_xticklabels(new_ticks) #set mag ticks
ax2.minorticks_off() #switch off log ticks on the upper axis
ax1.minorticks_on() #keep them on on the lower axis
ax2.set_xlabel('Mean $K$-band Magnitude (AB)')


plt.legend()

### Activate on click properties
fig.canvas.mpl_connect('pick_event', vari_funcs.selection_plot_funcs.onpickflux_2arcsec)

#
### Select Variables as those with chisq > 22.458 and >50 ###
varydata24 = tbdata[vary>22.458]
varydata30 = tbdata[vary>30]
varydata40 = tbdata[vary>40]
varydata50 = tbdata[vary>50]

xtalkvarydata30 = xtalkdata[varyxtalk>30]

plt.hlines(30, 4e0, 1e7,'g', zorder=4,label='Chi=30')

plt.legend()
plt.tight_layout()
