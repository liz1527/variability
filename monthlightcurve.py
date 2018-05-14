#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 19 12:07:06 2018

@author: ppxee
"""

### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
#from astropy.table import Table #for handling tables
import numpy as np #for handling arrays
#import math
import vari_funcs #my module to help run code neatly
plt.close('all') #close any open plots


tbdata = fits.open('mag_flux_tables/month_mag_flux_table_best_err3.fits')[1].data
obnum = 291878


#mask = fitsdata['NUMBER_1'] == ob
obdata = tbdata[tbdata['NUMBER'] == obnum]

months = ['sep05','oct05','nov05','dec05', 'jan06', 'dec06', 'jan07',  
          'aug07', 'sep07', 'oct07', 'sep08', 'oct08', 'nov08', 'jul09',  
          'aug09', 'sep09', 'oct09', 'nov09', 'dec09', 'jan10', 'feb10', 
          'aug10', 'sep10', 'oct10', 'nov10', 'dec10', 'jan11', 'feb11', 
          'aug11', 'sep11', 'oct11', 'nov11', 'dec11', 'jan12', 'feb12', 
          'jul12', 'aug12', 'sep12', 'oct12', 'nov12']

for month in months:
    if month == 'sep05':
        flux = obdata['MAG_APER_'+month][:,4]
        fluxerr = obdata['MAGERR_APER_'+month][:,4]
    else:
        flux = np.append(flux, obdata['MAG_APER_'+month][:,4])
        fluxerr = np.append(fluxerr, obdata['MAGERR_APER_'+month][:,4])

flux[flux == 99] = np.nan
fluxerr[flux==99] = np.nan
fluxerrnew = fluxerr + obdata['errordiff']

avgflux = np.nanmean(flux, axis=0)
diff = avgflux - 1
fluxnorm = flux - diff
fluxerrnorm = fluxerr # fluxnorm * (fluxerr/flux)
fluxerrnormnew = fluxerrnew # fluxnorm * (fluxerrnew/flux)
avgfluxnorm = np.nanmean(fluxnorm, axis=0)

#baseerrsq = np.square(fluxerr)
meanerr = np.nanmean(fluxerr, axis=0)
N = np.size(flux, axis=0)
#numobs = np.size(flux, axis=1)
sig = ((flux- avgflux)**2 - meanerr**2)# 
sigsum = np.nansum(sig)
normsig = sigsum/N

meanerr2 = np.nanmean(fluxerrnew, axis=0)
N = np.size(flux, axis=0)
#numobs = np.size(flux, axis=1)
sig2 = ((flux- avgflux)**2 - meanerr2**2)# 
sigsum2 = np.nansum(sig2)
normsig2 = sigsum2/N

meanerr3 = np.nanmean(fluxerrnorm, axis=0)
N = np.size(fluxnorm, axis=0)
#numobs = np.size(flux, axis=1)
sig3 = ((fluxnorm- avgfluxnorm)**2 - meanerr3**2)# 
sigsum3 = np.nansum(sig3)
normsig3 = sigsum3/N

meanerr4 = np.nanmean(fluxerrnormnew, axis=0)
N = np.size(fluxnorm, axis=0)
#numobs = np.size(flux, axis=1)
sig4 = ((fluxnorm - avgfluxnorm)**2 - meanerr4**2)# 
sigsum4 = np.nansum(sig4)
normsig4 = sigsum4/N

#vari_funcs.lightcurve5months(obnum, tbdata)
#set up time variable for plot
nums = fits.open('monthly_numbers.fits')[1].data
t = np.linspace(1, len(nums), num=len(nums))
tdataind = np.isin(nums['Month'], months)
tdata = t[tdataind]

#Plot graph in new figure
plt.figure(figsize=[17,7])
plt.xticks(t, nums['Month'], rotation='vertical')
plt.errorbar(tdata, flux, yerr=fluxerr, fmt = 'ro', zorder=1)
plt.xlabel('Month')
plt.ylabel('K-band magnitude of object')
plt.title('Light curve for object number %i' % obnum)
plt.tight_layout()
axes = plt.gca()
ylims = axes.get_ylim()
ymid = (ylims[1]+ylims[0])/2
#plt.ylim(ymin=ymid-0.25, ymax=ymid+0.25)

#Plot graph in new figure
#plt.figure(figsize=[17,7])
plt.xticks(t, nums['Month'], rotation='vertical')
plt.errorbar(tdata, flux, yerr=fluxerrnew, fmt = 'bo', zorder=0)
plt.xlabel('Month')
plt.ylabel('K-band magnitude of object')
plt.title('Light curve for object number %i' % obnum)
plt.tight_layout()
axes = plt.gca()
ylims = axes.get_ylim()
ymid = (ylims[1]+ylims[0])/2
#plt.ylim(ymin=19.28-0.25, ymax=19.28+0.25)
#plt.savefig(str(obnum)+'month.png')


#Plot graph in new figure
plt.figure(figsize=[17,7])
plt.xticks(t, nums['Month'], rotation='vertical')
plt.errorbar(tdata, fluxnorm, yerr=fluxerrnorm, fmt = 'ro', zorder=1)
plt.xlabel('Month')
plt.ylabel('K-band magnitude of object')
plt.title('Light curve for object number %i' % obnum)
plt.tight_layout()
axes = plt.gca()
ylims = axes.get_ylim()
ymid = (ylims[1]+ylims[0])/2
#plt.ylim(ymin=ymid-0.25, ymax=ymid+0.25)

#Plot graph in new figure
#plt.figure(figsize=[17,7])
plt.xticks(t, nums['Month'], rotation='vertical')
plt.errorbar(tdata, fluxnorm, yerr=fluxerrnormnew, fmt = 'bo', zorder=0)
plt.xlabel('Month')
plt.ylabel('K-band magnitude of object')
plt.title('Light curve for object number %i' % obnum)
plt.tight_layout()
axes = plt.gca()
ylims = axes.get_ylim()
ymid = (ylims[1]+ylims[0])/2