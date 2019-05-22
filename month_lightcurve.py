#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 13:53:26 2019

code to create month lightcurves

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
#from astropy.cosmology import FlatLambdaCDM
#from astropy import units as u
plt.close('all') #close any open plots

def month_avg_lightcurve(avgflux, avgfluxerr):
    months = ['sep05','oct05','nov05','dec05', 'jan06', 'dec06', 'jan07',  
          'aug07', 'sep07', 'oct07', 'sep08', 'oct08', 'nov08', 'jul09',  
          'aug09', 'sep09', 'oct09', 'nov09', 'dec09', 'jan10', 'feb10', 
          'aug10', 'sep10', 'oct10', 'nov10', 'dec10', 'jan11', 'feb11', 
          'aug11', 'sep11', 'oct11', 'nov11', 'dec11', 'jan12', 'feb12', 
          'jul12', 'aug12', 'sep12', 'oct12', 'nov12']
   
       
    #set up time variable for plot
    nums = fits.open('monthly_numbers.fits')[1].data
    t = np.linspace(1, len(nums), num=len(nums))
    tdataind = np.isin(nums['Month'], months)
    tdata = t[tdataind]
    
    ticks = nums['Month']
    mask = np.zeros(91)
    inds = [0,4,14,16,23,26,36,38,46,53,59,65,71,77,82,86]
    mask[inds] = 1
    mask = mask.astype(bool)
    ticks[~mask] = ''
    
    #Plot graph in new figure
    plt.figure(figsize=[17,7])
    plt.xticks(t, ticks, rotation='vertical')
    plt.errorbar(tdata, avgflux, yerr=avgfluxerr, fmt = 'bo')
    plt.xlabel('Month')
    plt.ylabel('K-band flux of object')
    plt.title('Average Lightcurve')
    plt.tight_layout()
    return



tbdata = fits.open('mag_flux_tables/month_mag_flux_table.fits')[1].data
obnum = 102552


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
        
#    if month == 'sep05':
#        flux = obdata['FLUX_APER_'+month][:,0]
#        fluxerr = obdata['FLUXERR_APER_'+month][:,0]
#    else:
#        flux = np.append(flux, obdata['FLUX_APER_'+month][:,0])
#        fluxerr = np.append(fluxerr, obdata['FLUXERR_APER_'+month][:,0])
        
mask = flux == 99
#mask = flux <= 0
flux[mask] = np.nan
fluxerr[mask] = np.nan
month_avg_lightcurve(flux, fluxerr)

chisq = vari_funcs.my_chisquare_err([flux], [fluxerr])

#plt.ylim(ymin=40000)
#axes = plt.gca()
#ylims = axes.get_ylim()
#ymid = (ylims[1]+ylims[0])/2
#plt.ylim(ymin=ymid-0.25, ymax=ymid+0.25)
#plt.title('Light curve for object number %i' % obnum)
plt.ylabel('K-band magnitude')
plt.title('Lightcurve of Object '+str(obnum)+' '+r' $\chi^{2} = $'+str(round(chisq[0], 2)))
