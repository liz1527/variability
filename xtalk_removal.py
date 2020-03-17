#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  8 16:44:00 2018

Code to identify variables that are within 2" of an X-talk with a predicted
flux that is >10% of the actual object flux.

@author: ppxee
"""

### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
from astropy.table import Table, Column #for handling tables
import numpy as np #for handling arrays
from astropy.stats import median_absolute_deviation
import vari_funcs #my module to help run code neatly
from astropy import units as u
plt.close('all') #close any open plots

font = {'family' : 'DejaVu Sans',
        'weight' : 'normal',
        'size'   : 16}

plt.rc('font', **font)

def mag2flux(mag):
    flux = 10**((mag-30)/(-2.5))
    return flux

def get_percentage(xtalk, tbdata):
    
    ### Get DR11 K mag and convert ###
    kmag = tbdata['KMAG_20']
    kflux = mag2flux(kmag)
    
    ### Get predicted magnitude ###
    predmag = xtalk['magpred']
    predflux = mag2flux(predmag)
    
    ### Get % ###
    percent = (predflux/kflux) * 100
    
    return percent, kflux, predflux
    
### Open the fits files and get data ###
varys = Table.read('variable_tables/K/variables_no06_chi30_neg_DR11data.fits')
xtalk = Table.read('UDS_catalogues/UDS_DR11_pred_xtalk_mag.fits')

#%% match variables and xtalk catalogue ###
### Get coordinates ###
#x_varys = varys['X_IMAGE']
#y_varys = varys['Y_IMAGE']
x_xtalk = xtalk['X']
y_xtalk = xtalk['Y']

### Set match distance ###
dist = 14.9 # 2 arcsec dist in pixels 

### find dist from x-talks for all objects ###
tot = 0 #keep track of total number of matches
bad = 0 #keep track of number that do not meet criteria
contam = np.empty(len(varys), dtype=bool) # array to hold boolean of object contamination
contam_col = Column(contam, name='X-talk_contam')
varys.add_column(contam_col)

varys_noxtalk = []
for n, id in enumerate(varys['ID']):
    obdata = varys[varys['ID']==id]
    x_varys = obdata['X_IMAGE']
    y_varys = obdata['Y_IMAGE']
    
    ### define radial distance from variable to all x talks ###
    r = np.sqrt(np.square(x_varys-x_xtalk) + np.square(y_varys-y_xtalk))
    
    ### identify any xtalks within 2 arcsec ###
    mask = r <= dist
    near_xtalk = xtalk[mask]

    if len(near_xtalk) != 0:
        tot += 1
        ### Get percentage and fluxes ###
        percent, kflux, predflux = get_percentage(near_xtalk, obdata)
        
        mask = percent > 10 # mask those with >10% flux
        near_xtalk = near_xtalk[mask]
        if len(near_xtalk) != 0:
            bad += 1
            varys['X-talk_contam'][n] = True
            obdata['X-talk_contam'] = True
            if bad == 1:
                varys_xtalk = obdata#near_xtalk
            else:
                varys_xtalk.add_row(obdata[0])#near_xtalk[0])
#            print('ID = '+str(id)+' is near '+str(len(near_xtalk))+' X-talks with >10% flux')
        else:
            varys['X-talk_contam'][n] = False
            obdata['X-talk_contam'] = False
            if varys_noxtalk == []:
                varys_noxtalk = obdata
            else:
                varys_noxtalk.add_row(obdata[0])
    else:
        varys['X-talk_contam'][n] = False
        obdata['X-talk_contam'] = False
        if varys_noxtalk == []:
            varys_noxtalk = obdata
        else:
            varys_noxtalk.add_row(obdata[0])
            
varys.write('variable_tables/K/variables_no06_chi30_neg_DR11data_xtalkchecked.fits')
varys_xtalk.write('variable_tables/K/variables_no06_chi30_neg_DR11data_xtalkcontam.fits')
varys_noxtalk.write('variable_tables/K/variables_no06_chi30_neg_DR11data_noxtalkcontam.fits')
































