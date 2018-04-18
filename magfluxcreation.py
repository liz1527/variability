#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 15:43:38 2018

Code to try and streamline the creation of mag-flux tables for variability
analysis.

Update on 13/4/18 to allow creations of month stacks mag_flux

@author: ppxee
"""

from astropy.table import Table, join, vstack, Column
import numpy as np

#%% Create Semester mag-flux tables
#sems = ['05B', '06B', '07B', '08B', '09B', '10B', '11B', '12B']
#
#def sem_mag_flux(sem):
##    if sem == '12B':
##        semtb = Table.read('SE_outputs_yearstacks/newq_'+sem+'_output.fits')
##    else:
##        semtb = Table.read('SE_outputs_yearstacks/newS_newq_'+sem+'_output.fits')
#    semtb = Table.read('SE_outputs_yearstacks/'+sem+'_output.fits')
#    #extract column names
#    cols = np.asarray(semtb.colnames)
#    
#    #create a mask to identify the columns that are needed
#    colmask = np.zeros(64) +1
#    colmask[0:5] = 0 # keep ID, RA, DEC, and x and y coords
#    colmask[40] = 0 # keep flux radius column
#    colmask[57:61] = 0 # keep aperture columns
#    colmask[63] = 0 # keep FWHM world
#    colmask = colmask.astype(bool)
#    
#    #create an array of the column names that arent wanted
#    badnames = cols[colmask]
#    
#    # iterate through names to remove columns
#    for name in badnames:
#        del semtb[name]
#        
#    for name in semtb.colnames:
#        semtb.rename_column(name, name+'_'+sem)  
#        
#    return semtb
#
#%% Code for creating mag flux tables with month stack data
def month_mag_flux(month):
    semtb = Table.read('SE_outputs_monthstacks/'+month+'_output.fits')
    #extract column names
    cols = np.asarray(semtb.colnames)
    
    #create a mask to identify the columns that are needed
    colmask = np.zeros(64) +1
    colmask[0:5] = 0 # keep ID, RA, DEC, and x and y coords
    colmask[40] = 0 # keep flux radius column
    colmask[57:61] = 0 # keep aperture columns
    colmask[63] = 0 # keep FWHM world
    colmask = colmask.astype(bool)
    
    #create an array of the column names that arent wanted
    badnames = cols[colmask]
    
    # iterate through names to remove columns
    for name in badnames:
        del semtb[name]
        
    for name in semtb.colnames:
        if name == 'NUMBER':
            continue
        semtb.rename_column(name, name+'_'+month)  
        
    return semtb

#%% 
    
months = ['sep05','oct05','nov05','dec05', 'jan06', 'dec06', 'jan07', 'aug07', 
          'sep07', 'oct07', 'sep08', 'oct08', 'nov08', 'jul09', 'aug09', 
          'sep09', 'oct09', 'nov09', 'dec09', 'jan10', 'feb10', 'aug10',
          'sep10', 'oct10', 'nov10', 'dec10', 'jan11', 'feb11', 'aug11',
          'sep11', 'oct11', 'nov11', 'dec11', 'jan12', 'feb12', 'jul12',
          'aug12', 'sep12', 'oct12', 'nov12']

for i, month  in enumerate(months):
    #create month table with correctly renamed columns
    montb = month_mag_flux(month)
    print('Join ' + str(i))
    if i == 0:
        semcom = montb # for the first month just make the combined table the 
                        #month table
    else:
        semcom = join(semcom, montb, keys='NUMBER') #join tables together using 
                                                    # IDs in my catalogue
    del montb
    

#%% Match these with various catalogs to create final tables
from astropy.coordinates import match_coordinates_sky
from astropy.coordinates import SkyCoord
from astropy import units as u

# match with stars catalogue
print('Matching Stars')
stars = Table.read('UDS_catalogues/DR8-secure-stars.fits')
starscoord = SkyCoord(stars['RA']*u.degree, stars['DEC']*u.degree)
semcomcoord = SkyCoord(semcom['ALPHA_J2000_sep05'], semcom['DELTA_J2000_sep05'])
idx, d2d , _ = match_coordinates_sky(starscoord, semcomcoord)
mask = d2d<=1*u.arcsec #make sure match is within 1 arcsec (like in topcat)
idx = idx[mask]
ind = np.arange(len(semcom)) #create array of indicies
ind = np.delete(ind, idx) #remove those that relate to stars
starsmf = semcom[idx] #create table of just stars
semcomns = semcom[ind] #create table of no stars

#match with best catalogue
print('Matching Best')
best = Table.read('UDS_catalogues/uds_multicat_ap3.v5b_best.fits')
bestcoord = SkyCoord(best['ALPHA_J2000']*u.degree, best['DELTA_J2000']*u.degree)
semcomnscoord = SkyCoord(semcomns['ALPHA_J2000_sep05'], semcomns['DELTA_J2000_sep05'])
idx, d2d , _ = match_coordinates_sky(bestcoord, semcomnscoord)
mask = d2d<=1*u.arcsec #make sure match is within 1 arcsec (like in topcat)
idx = idx[mask]
bestmf = semcomns[idx] #create best table with no stars

bestmfcoord = SkyCoord(bestmf['ALPHA_J2000_sep05'], bestmf['DELTA_J2000_sep05'])

# match with xmm
print('Matching XMM')
xmm = Table.read('UDS_catalogues/XMM_not_Chandra.fits')
xmmcoord = SkyCoord(xmm['RAJ2000'], xmm['DEJ2000'])
idx, d2d , _ = match_coordinates_sky(xmmcoord, bestmfcoord)
mask = d2d<=5*u.arcsec #make sure match is within 5 arcsec (like in topcat)
idx = idx[mask]
xmmmf = bestmf[idx]

## match with chandra
print('Matching Chandra')
chan = Table.read('UDS_catalogues/chandra_catalogue.fits')
chan['RA'].unit = u.deg
chan['Dec'].unit = u.deg
chancoord = SkyCoord(chan['RA'], chan['Dec'])
idx, d2d , _ = match_coordinates_sky(chancoord, bestmfcoord)
mask = d2d<=1*u.arcsec #make sure match is within 1 arcsec (like in topcat)
idx = idx[mask]
chanmf = bestmf[idx]

# combine chandra and xmm
print('Joining xray table')
xraymf = vstack([chanmf, xmmmf])
#%%
# boolean whether a source is seen in x-rays
xray = np.isin(bestmf['NUMBER'], xraymf['NUMBER'])
xraycol = Column(xray, 'X-ray')
bestmf.add_column(xraycol)

#%% Save the tables
semcom.write('month_mag_flux_table.fits')
bestmf.write('month_mag_flux_table_best.fits')
starsmf.write('month_stars_mag_flux_table.fits')
xraymf.write('month_xray_mag_flux_table_best.fits')