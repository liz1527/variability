#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  8 16:44:00 2018

Code to look into properties of potential xtalks
@author: ppxee
"""

### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
from astropy.table import Table #for handling tables
import numpy as np #for handling arrays
from astropy.stats import median_absolute_deviation
import vari_funcs #my module to help run code neatly
plt.close('all') #close any open plots

font = {'family' : 'DejaVu Sans',
        'weight' : 'normal',
        'size'   : 16}

plt.rc('font', **font)

def mag2flux(mag):
    flux = 10**((mag-30)/(-2.5))
    return flux

def get_percentage(tbdata):
    
    ### Get DR11 K mag and convert ###
    kmag = tbdata['KMAG_20']
    kflux = mag2flux(kmag)
    
    ### Get predicted magnitude ###
    predmag = tbdata['magpred']
    predflux = mag2flux(predmag)
    
    ### Get % ###
    percent = (predflux/kflux) * 100
    
    return percent, kflux, predflux
    
### Open the fits files and get data ###
xtalk_1_noneg_varys = fits.open('variable_tables/K/variables_no06_chi30_DR11data_megaxtalk_1arcsec.fits')[1].data
xtalk_2_noneg_varys = fits.open('variable_tables/K/variables_no06_chi30_DR11data_megaxtalk_2arcsec.fits')[1].data
xtalk_1_neg_varys = fits.open('variable_tables/K/variables_no06_chi30_neg_DR11data_megaxtalk_1arcsec.fits')[1].data
xtalk_2_neg_varys = fits.open('variable_tables/K/variables_no06_chi30_neg_DR11data_megaxtalk_2arcsec.fits')[1].data

### Get percentage and fluxes ###

xtalk_1_noneg_percent, xtalk_1_noneg_kflux, xtalk_1_noneg_predflux = get_percentage(xtalk_1_noneg_varys)
xtalk_2_noneg_percent, xtalk_2_noneg_kflux, xtalk_2_noneg_predflux = get_percentage(xtalk_2_noneg_varys)
xtalk_1_neg_percent, xtalk_1_neg_kflux, xtalk_1_neg_predflux = get_percentage(xtalk_1_neg_varys)
xtalk_2_neg_percent, xtalk_2_neg_kflux, xtalk_2_neg_predflux = get_percentage(xtalk_2_neg_varys)

#### plot percentage distribution ###
#plt.figure()
##bins = np.linspace(0,100,50)
#bins = np.logspace(0,2,20)
##plt.hist(xtalk_1_noneg_percent, histtype='step')
##plt.hist(xtalk_2_noneg_percent, histtype='step')
##plt.hist(xtalk_1_neg_percent, histtype='step')
##plt.hist(xtalk_2_neg_percent, histtype='step')
#plt.hist(xtalk_1_noneg_percent, bins, histtype='step')
#plt.hist(xtalk_2_noneg_percent, bins, histtype='step')
#plt.hist(xtalk_1_neg_percent, bins, histtype='step')
#plt.hist(xtalk_2_neg_percent, bins, histtype='step')
#plt.xscale('log')

### plot flux vx predicted flux ###
x = np.logspace(0,5)
y = x
y_10 = 0.1*x
y_5 = 0.05*x

plt.figure()
plt.plot(xtalk_1_noneg_kflux, xtalk_1_noneg_predflux,'o')
#plt.scatter(xtalk_1_noneg_kflux, xtalk_1_noneg_predflux, c=xtalk_1_noneg_percent)
plt.plot(x,y,'k', label='100%')
plt.plot(x,y_10,'b', label='10%')
plt.plot(x,y_5,'r', label='5%')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Object flux')
plt.ylabel('Pred X-talk flux')
plt.title('No Negatives 1" Match')
plt.xlim(1,1e5)
plt.ylim(0.1,1e5)
plt.legend()
plt.tight_layout()

plt.figure()
plt.plot(xtalk_2_noneg_kflux, xtalk_2_noneg_predflux,'o')
#plt.scatter(xtalk_2_noneg_kflux, xtalk_2_noneg_predflux, c=xtalk_2_noneg_percent)
plt.plot(x,y,'k', label='100%')
plt.plot(x,y_10,'b', label='10%')
plt.plot(x,y_5,'r', label='5%')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Object flux')
plt.ylabel('Pred X-talk flux')
plt.title('No Negatives 2" Match')
plt.xlim(1,1e5)
plt.ylim(0.1,1e5)
plt.legend()
plt.tight_layout()

plt.figure()
plt.plot(xtalk_1_neg_kflux, xtalk_1_neg_predflux,'o')
#plt.scatter(xtalk_1_neg_kflux, xtalk_1_neg_predflux, c=xtalk_1_neg_percent)
plt.plot(x,y,'k', label='100%')
plt.plot(x,y_10,'b', label='10%')
plt.plot(x,y_5,'r', label='5%')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Object flux')
plt.ylabel('Pred X-talk flux')
plt.title('Incl. Negatives 1" Match')
plt.xlim(1,1e5)
plt.ylim(0.1,1e5)
plt.legend()
plt.tight_layout()

plt.figure()
plt.plot(xtalk_2_neg_kflux, xtalk_2_neg_predflux,'o')
#plt.scatter(xtalk_2_neg_kflux, xtalk_2_neg_predflux, c=xtalk_2_neg_percent)
plt.plot(x,y,'k', label='100%')
plt.plot(x,y_10,'b', label='10%')
plt.plot(x,y_5,'r', label='5%')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Object flux')
plt.ylabel('Pred X-talk flux')
plt.title('Incl. Negatives 2" Match')
plt.xlim(1,1e5)
plt.ylim(0.1,1e5)
plt.legend()
plt.tight_layout()




