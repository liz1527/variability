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
#filename = 'J/K_extraction/J_variables_chi32_noneg_DR11data'
filename = 'bad_jan12_DR11data'
varys = Table.read('variable_tables/'+filename+'.fits')
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
        percent = percent[mask]
        if len(near_xtalk) != 0:
            bad += 1
            varys['X-talk_contam'][n] = True
            obdata['X-talk_contam'] = True
            if bad == 1:
                varys_xtalk = obdata#near_xtalk
                num = len(near_xtalk)
                if num == 1:
                    bad_percents = percent
                    xtalk_x = near_xtalk['X']
                    xtalk_y = near_xtalk['Y']
                else:
                    maxpercent = np.argmax(percent)
                    bad_percents = percent[maxpercent]
                    xtalk_x = near_xtalk['X'][maxpercent]
                    xtalk_y = near_xtalk['Y'][maxpercent]
            else:
                varys_xtalk.add_row(obdata[0])#near_xtalk[0])
                num = np.append(num, len(near_xtalk))
                if num[bad-1] == 1:
                    bad_percents = np.append(bad_percents, percent)
                    xtalk_x = np.append(xtalk_x, near_xtalk['X'])
                    xtalk_y = np.append(xtalk_y, near_xtalk['Y'])
                else:
                    maxpercent = np.argmax(percent)
                    bad_percents = np.append(bad_percents, percent[maxpercent])
                    xtalk_x = np.append(xtalk_x, near_xtalk['X'][maxpercent])
                    xtalk_y = np.append(xtalk_y, near_xtalk['Y'][maxpercent])
            print('ID = '+str(id)+' is near '+str(len(near_xtalk))+' X-talks with >10% flux')
        else:
            varys['X-talk_contam'][n] = False
            obdata['X-talk_contam'] = False
            if len(varys_noxtalk) == 0:
                varys_noxtalk = obdata
            else:
                varys_noxtalk.add_row(obdata[0])
    else:
        varys['X-talk_contam'][n] = False
        obdata['X-talk_contam'] = False
        if len(varys_noxtalk) == 0:
            varys_noxtalk = obdata
        else:
            varys_noxtalk.add_row(obdata[0])
            
### Add number, percent and coord columns to contaminated table ###
num_col = Column(num, name='X-talk_num')
varys_xtalk.add_column(num_col)
percent_col = Column(bad_percents, name='X-talk_percent')
varys_xtalk.add_column(percent_col)
x_col = Column(xtalk_x, name='X-talk_X')
varys_xtalk.add_column(x_col)
y_col = Column(xtalk_y, name='X-talk_Y')
varys_xtalk.add_column(y_col)
        

#### Save tables ###
#varys.write('variable_tables/'+filename+'_xtalkchecked.fits',
#            overwrite=True)
#varys_xtalk.write('variable_tables/'+filename+'_xtalkcontam.fits',
#            overwrite=True)
#varys_noxtalk.write('variable_tables/'+filename+'_noxtalkcontam.fits',
#            overwrite=True)
































