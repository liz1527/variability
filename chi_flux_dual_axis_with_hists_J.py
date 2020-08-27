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
tbdata = fits.open('mag_flux_tables/J/mag_flux_table_best_J_extra_clean.fits')[1].data
chandata = fits.open('mag_flux_tables/J/xray_mag_flux_table_best_J_extra_clean.fits')[1].data
sdata = fits.open('mag_flux_tables/J/stars_mag_flux_table_J_extra_clean.fits')[1].data
sigtb = Table.read('sigma_tables/quad_epoch_sigma_table_J_extra_clean_2arcsec_noneg.fits')

def prep_data(tbdata):
    ### Remove edges ###
    tbdata = vari_funcs.field_funcs.remove_edges(tbdata)
    
    ### Create arrays of flux values ###
    flux = vari_funcs.j_mag_flux.flux4_stacks(tbdata)
    
    ### remove values that are negative ###
    flux, tbdata = vari_funcs.flux_funcs.noneg(flux, tbdata)
    
    ### Get error arrays ###
    flux, fluxerr, tbdata = vari_funcs.j_mag_flux.create_quad_error_array_J(sigtb, 
                                                                            tbdata, 
                                                                            aper=4)
    
    ### get average flux for each object ###
    avgfluxperob = np.nanmean(flux, axis=1)
    
    ### find chi squared ###
    vary = vari_funcs.vary_stats.my_chisquare_err(flux, fluxerr)
    
    return tbdata, avgfluxperob, vary

tbdata, avgfluxperob, vary = prep_data(tbdata)
chandata, avgfluxchanperob, varychan = prep_data(chandata)
sdata, savgfluxperob, varystar = prep_data(sdata)

### reset X-ray column as messed up by stacking ###
tbdata['X-ray'][tbdata['X-ray']==70] = False 
tbdata['X-ray'][tbdata['X-ray']==84] = True

### Create plot with two axes ##
fig = plt.figure(figsize=[9,8])
ax1 = plt.subplot2grid((5,5), (0,0), colspan=4, rowspan=5)# Subplot covering majority of plot
ax4 = plt.subplot2grid((5,5), (0,4), rowspan=5, sharey=ax1)

ax2 = ax1.twiny() # Twin of the first subplot

### set tick values ###
low_ticks = np.array([1e1, 1e2, 1e3, 1e4, 1e5, 1e6])#, 1e7]) #set flux ticks
#new_ticks = 30 - 2.5*np.log10(ticks) #find equivilant mag ticks
#new_ticks = new_ticks + 1.9

new_ticks =np.array([29, 27, 25, 23, 21, 19, 17, 15])
upp_ticks = new_ticks - 1.9
upp_ticks = 10**((upp_ticks-30)/(-2.5))

### Plot the chi sq v mean ###
ax1.plot(savgfluxperob, varystar, 'm*', mfc = 'none', markersize = 10,
             label='DR11 Star', rasterized=True) 
line, = ax1.plot(avgfluxperob, vary, 'b+', label='Galaxy', picker=2, rasterized=True)
ax1.plot(avgfluxchanperob, varychan, 'ro', mfc = 'none', markersize = 10,
         label='X-ray detected', rasterized=True) #no picker as will be selected in the UDS point

### Plot chi sq hist ###
bins = np.logspace(np.log10(3e-2), np.log10(3e4),100)
ax4.hist(varystar, bins, color='m', histtype='step', linestyle='--', 
         linewidth=1.5, density=True, orientation=u'horizontal')
ax4.hist(vary, bins, color='b', histtype='step', linestyle='-.', 
         linewidth=1.5, density=True, orientation=u'horizontal')
ax4.hist(varychan, bins, color='r', histtype='step', linestyle='-', 
         linewidth=1.5, density=True, orientation=u'horizontal')

### Apply plot characteristics ###
ax1.set_xscale('log')
ax1.set_yscale('log')
    
ax1.set_ylim(3e-2,3e4)
ax1.set_xlim(1e1, 1e7) #if flux
#    plt.xlim(13,26) #if mag

ax1.set_xlabel('Mean Flux')
ax1.set_ylabel(r'$\chi^{2}$')


ax2.set_xlim(ax1.get_xlim()) #set twin to same limits
ax2.set_xscale('log') #set twin to same scale
ax2.set_xticks(upp_ticks) #set flux ticks
ax2.set_xticklabels(new_ticks) #set mag ticks
ax2.minorticks_off() #switch off log ticks on the upper axis
ax2.set_xlabel('Mean $J$-band Magnitude (AB)')

ax1.set_xticks(low_ticks) #set flux ticks
ax1.minorticks_on() #keep them on on the lower axis

plt.setp(ax4.get_yticklabels(), visible=False)
ax4.set_xticks(ax4.get_xticks()[1:]) 
#ax4.tick_params(axis='y', which='both', left=False)
#ax4.majorticks_off() #switch off log ticks on the upper axis
ax4.set_xlabel('Norm. Counts')

### Plot flux hist ###
#bins2 = np.logspace(np.log10(8e1), np.log10(1e7),100)
#ax3.hist(savgfluxperob, bins2, color='m', histtype='step', linestyle='--', 
#         linewidth=1.5, density=True)
#ax3.hist(avgfluxperob, bins2, color='b', histtype='step', linestyle='-.', 
#         linewidth=1.5, density=True)
#ax3.hist(avgfluxchanperob, bins2, color='r', histtype='step', linestyle='-', 
#         linewidth=1.5, density=True)
#
#### Apply plot characteristics ###
#plt.setp(ax3.get_xticklabels(), visible=False)

### Activate on click properties
fig.canvas.mpl_connect('pick_event', vari_funcs.selection_plot_funcs.onpickflux_2arcsec)


### Select Variables as those with chisq > 32.08 and >50 ###
varydata24 = tbdata[vary>32.08]
varydata30 = tbdata[vary>30]
varydata40 = tbdata[vary>40]
varydata50 = tbdata[vary>50]
#
#svarydata30 = sdata[varystar>30]
#save30 = Table(svarydata30)
#save30.write('variable_tables/K/variables_stars.fits')

ax1.hlines(32.08, 8e1, 1e7,'g', zorder=4,label='Chi=32.08')

ax1.legend(fontsize=14)
#plt.subplots_adjust(hspace=0)
plt.tight_layout(w_pad=0, h_pad=0.6)
