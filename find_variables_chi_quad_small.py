#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 11:35:11 2018

Code to isolate variable sources using the bootstraped error bars and chi^2

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

### Open the fits files and get data ###
tbdata = fits.open('mag_flux_tables/K/mag_flux_table_best_K_extra_clean.fits')[1].data
chandata = fits.open('mag_flux_tables/K/xray_mag_flux_table_best_K_extra_clean.fits')[1].data
sdata = fits.open('mag_flux_tables/K/stars_mag_flux_table_K_extra_clean.fits')[1].data
sigtb = Table.read('sigma_tables/quad_epoch_sigma_table_K_extra_clean_1arcsec_noneg.fits')

ap = 2 # 1arcsec aper

### Remove edges ###
tbdata = vari_funcs.field_funcs.remove_edges(tbdata)
chandata = vari_funcs.field_funcs.remove_edges(chandata)
sdata = vari_funcs.field_funcs.remove_edges(sdata)

## Create arrays of flux values ###
flux = vari_funcs.k_mag_flux.flux_stacks(tbdata, aper=ap)
fluxchan = vari_funcs.k_mag_flux.flux_stacks(chandata, aper=ap) 
sflux = vari_funcs.k_mag_flux.flux_stacks(sdata, aper=ap)

### remove values that are negative ###
flux, tbdata = vari_funcs.flux_funcs.noneg(flux, tbdata)
fluxchan, chandata = vari_funcs.flux_funcs.noneg(fluxchan, chandata)
sflux, sdata = vari_funcs.flux_funcs.noneg(sflux, sdata)

#flux, tbdata = vari_funcs.flux_funcs.nremove_low_flux(flux, tbdata)
#fluxchan, chandata = vari_funcs.flux_funcs.nremove_low_flux(fluxchan, chandata)
#sflux, sdata = vari_funcs.flux_funcs.nremove_low_flux(sflux, sdata)

### Get error arrays ###
flux, fluxerr, tbdata = vari_funcs.k_mag_flux.create_quad_error_array(sigtb, tbdata, aper=ap)
fluxchan, chanerr, chandata = vari_funcs.k_mag_flux.create_quad_error_array(sigtb, chandata,aper=ap)
sflux, serr, sdata = vari_funcs.k_mag_flux.create_quad_error_array(sigtb, sdata, aper=ap)

### reset X-ray column as messed up by stacking ###
tbdata['X-ray'][tbdata['X-ray']==70] = False 
tbdata['X-ray'][tbdata['X-ray']==84] = True

### Check chisq plot looks correct ###
fig,_ = vari_funcs.selection_plot_funcs.flux_variability_plot(flux, fluxchan, 'chisq', 
                                       fluxerr=fluxerr, chanerr=chanerr,
                                       starflux=sflux, starfluxerr=serr,
                                       #normalised=True, 
                                       stars=True, scale='log')
fig.canvas.mpl_connect('pick_event', vari_funcs.selection_plot_funcs.onpickflux_2arcsec)

plt.ylim(3e-2,3e4)
#plt.xlim(8e1, 1e7) #if flux noneg
plt.xlim(4e0, 1e7) #if flux neg
#devdata = fits.open('variable_tables/no06_variables_chi30_2arcsec_deviant_neg.fits')[1].data
#devflux, devfluxerr, devdata = vari_funcs.k_mag_flux.create_quad_error_array(sigtb, devdata, aper=ap)
#devmean = np.nanmean(devflux, axis=1)

### Calculate chi^2 values ###
chisq = vari_funcs.vary_stats.my_chisquare_err(flux, fluxerr)
chanchisq = vari_funcs.vary_stats.my_chisquare_err(fluxchan, chanerr)
#devchisq = vari_funcs.vary_stats.my_chisquare_err(devflux, devfluxerr)

#plt.plot(devmean, devchisq, 'kd', mfc='None', markersize=10)
### Select Variables as those with chisq > 22.458 and >50 ###
varydata24 = tbdata[chisq>22.458]
varydata30 = tbdata[chisq>30]
varydata40 = tbdata[chisq>40]
varydata50 = tbdata[chisq>50]

#plt.hlines(22.458, 8e1, 1e7,zorder=4,label='Chi>22.5')
plt.hlines(30, 4e0, 1e7, 'g', zorder=4,label='Chi=30')
#plt.hlines(40, 8e1, 1e7,'y', zorder=4,label='Chi>40')
#plt.hlines(50, 8e1, 1e7,'c', zorder=4,label='Chi>50')

plt.legend()
plt.tight_layout()

#### Check if variables vanish if the most deviant semester is removed ###
#varyflux = vari_funcs.k_mag_flux.flux5_stacks(varydata40)
#varyflux, varyfluxerr, varydata40 = vari_funcs.k_mag_flux.create_quad_error_array(sigtb, varydata40)
#fluxn, fluxnerr = vari_funcs.flux_funcs.normalise_flux_and_errors(varyflux, varyfluxerr)
#diff = abs(fluxn - 1)
#mask = np.argmax(diff, axis=1)
#varyfluxnew = np.copy(varyflux)
#varyfluxerrnew = np.copy(varyfluxerr)
#varyfluxnew[np.arange(len(varyflux)),mask] = np.nan
#varyfluxerrnew[np.arange(len(varyflux)),mask] = np.nan
#
#varychi = vari_funcs.vary_stats.my_chisquare_err(varyflux, varyfluxerr)
#varychinew = vari_funcs.vary_stats.my_chisquare_err(varyfluxnew, varyfluxerrnew)
#chidiff = varychi - varychinew
### Save new tables ###
#save24 = Table(varydata24)
#save24.write('variable_tables/no06_variables_chi22.fits')
#save30 = Table(varydata30)
#save30.write('variable_tables/K/variables_no06_chi30_1-5arcsec_neg.fits', overwrite=True)
#save40 = Table(varydata40)
#save40.write('variable_tables/no06_variables_chi40.fits')
#save50 = Table(varydata50)
#save50.write('variable_tables/no06_variables_chi50.fits')
#
#%% checks for false positives ###

### Create binedge array ###
bins = np.array(sigtb.colnames)
binarr = np.empty(int(len(bins)/4))
for k, bin in enumerate(bins):
    if bin[0] == '1':
        binarr[k] = int(bin[2:])
        k+=1
binarr = binarr.astype(int)

binsizes =  np.empty(int(len(binarr)))
binmean =  np.empty(int(len(binarr)))
for n, binedge in enumerate(binarr):
    if binedge==binarr[-1]:
        binupp = np.nanmax(flux)
    else:
        binupp = binarr[n+1]
    binflux, bindata = vari_funcs.flux_funcs.fluxbin(binedge, binupp, flux, tbdata) #bindata
#    binsflux, binsdata = vari_funcs.flux_funcs.fluxbin(binedge, binupp, sflux, sdata)
    binsizes[n] = len(bindata) #+ len(binsdata)
#    plt.vlines(binedge, 1e-2, 1e4, zorder=4)
    binmean[n] = np.nanmean(binflux)
#    if binsizes[n] < 9:
#        continue
    ### get chi values within bin ###
#    binflux, binfluxerr, bindata = vari_funcs.k_mag_flux.create_quad_error_array(sigtb, bindata)
#    binchisq = vari_funcs.vary_stat.my_chisquare_err(binflux, binfluxerr)
#    
    ### P value of 30 with dof=6 is 0.00003931 ###
    ### therefore false positives = binsize * P ###
P =  0.00003931#0.001
numfalpos = binsizes * P
plt.figure(figsize=[8,5])
plt.plot(binmean, numfalpos)
plt.xscale('log')
#plt.xlim(8e1,1e7)
plt.xlabel('Flux')
plt.ylabel('Expected number of sources with $\chi^{2}$ > 30')
plt.tight_layout()
totalfalpos = np.nansum(numfalpos)