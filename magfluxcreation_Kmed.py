#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 15:43:38 2018

Code to try and streamline the creation of mag-flux tables for variability
analysis.

Update on 13/4/18 to allow creations of month stacks mag_flux

@author: ppxee
"""
#%% Create Semester mag-flux tables
from astropy.table import Table, join, vstack, Column
import numpy as np

sems = ['0','1','2','3','4','5','6','7','8','9']

def sem_mag_flux(sem):
    semtb = Table.read('SE_outputs_yearstacks/08B_output_cam1_10bin'+sem+'_output.fits')
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
        semtb.rename_column(name, name+'_'+sem)  
        
    return semtb

# Create semester mag flux tables
sem0 = sem_mag_flux(sems[0])
sem1 = sem_mag_flux(sems[1])
sem2 = sem_mag_flux(sems[2])
sem3 = sem_mag_flux(sems[3])
sem4 = sem_mag_flux(sems[4])
sem5 = sem_mag_flux(sems[5])
sem6 = sem_mag_flux(sems[6])
sem7 = sem_mag_flux(sems[7])
sem8 = sem_mag_flux(sems[8])
sem9 = sem_mag_flux(sems[9])

#%% Join these tables using their IDs (will be fine for just this part as numbers alwasy match)
#print('Join 1')
#sem06B.rename_column('NUMBER_06B', 'NUMBER_05B')
#semcom = join(sem05B, sem06B, keys='NUMBER_05B')
semcom = sem0
print('Join 2')
sem1.rename_column('NUMBER_1', 'NUMBER_0')
semcom = join(semcom, sem1, keys='NUMBER_0')
print('Join 3')
sem2.rename_column('NUMBER_0', 'NUMBER_0')
semcom = join(semcom, sem2, keys='NUMBER_0')
print('Join 4')
sem3.rename_column('NUMBER_0', 'NUMBER_0')
semcom = join(semcom, sem3, keys='NUMBER_0')
print('Join 5')
sem4.rename_column('NUMBER_0', 'NUMBER_0')
semcom = join(semcom, sem4, keys='NUMBER_0')
print('Join 6')
sem5.rename_column('NUMBER_0', 'NUMBER_0')
semcom = join(semcom, sem5, keys='NUMBER_0')
print('Join 7')
sem6.rename_column('NUMBER_0', 'NUMBER_05B')
semcom = join(semcom, sem6, keys='NUMBER_05B')
print('Join 7')
sem7.rename_column('NUMBER_0', 'NUMBER_05B')
semcom = join(semcom, sem7, keys='NUMBER_05B')
print('Join 7')
sem8.rename_column('NUMBER_0', 'NUMBER_05B')
semcom = join(semcom, sem8, keys='NUMBER_05B')
print('Join 7')
sem9.rename_column('NUMBER_0', 'NUMBER_05B')
semcom = join(semcom, sem9, keys='NUMBER_05B')

#%% Match these with various catalogs to create final tables
from astropy.coordinates import match_coordinates_sky
from astropy.coordinates import SkyCoord
from astropy import units as u

# match with stars catalogue
print('Matching Stars')
stars = Table.read('UDS_catalogues/DR11-secure-stars.fits')
starscoord = SkyCoord(stars['RA']*u.degree, stars['DEC']*u.degree)
semcomcoord = SkyCoord(semcom['ALPHA_J2000_05B'], semcom['DELTA_J2000_05B'])
idx, d2d , _ = match_coordinates_sky(starscoord, semcomcoord)
mask = d2d<=1*u.arcsec #make sure match is within 1 arcsec (like in topcat)
idx = idx[mask]
ind = np.arange(len(semcom)) #create array of indicies
ind = np.delete(ind, idx) #remove those that relate to stars
starsmf = semcom[idx] #create table of just stars
semcomns = semcom[ind] #create table of no stars

#match with best catalogue
print('Matching Best')
best = Table.read('UDS_catalogues/DR11-2arcsec-Jan-1-2018_best.fits')
bestcoord = SkyCoord(best['RA']*u.degree, best['DEC']*u.degree)
semcomnscoord = SkyCoord(semcomns['ALPHA_J2000_05B'], semcomns['DELTA_J2000_05B'])
idx, d2d , _ = match_coordinates_sky(bestcoord, semcomnscoord)
mask = d2d<=1*u.arcsec #make sure match is within 1 arcsec (like in topcat)
idx = idx[mask]
bestmf = semcomns[idx] #create best table with no stars

bestmfcoord = SkyCoord(bestmf['ALPHA_J2000_05B'], bestmf['DELTA_J2000_05B'])

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
xray = np.isin(bestmf['NUMBER_05B'], xraymf['NUMBER_05B'])
xraycol = Column(xray, 'X-ray')
bestmf.add_column(xraycol)

#%% Save the tables
semcom.write('mag_flux_tables/mag_flux_table_extra_clean_no06.fits')
bestmf.write('mag_flux_tables/mag_flux_table_best_extra_clean_no06.fits')
starsmf.write('mag_flux_tables/stars_mag_flux_table_extra_clean_no06.fits')
xraymf.write('mag_flux_tables/xray_mag_flux_table_best_extra_clean_no06.fits')