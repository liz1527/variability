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
## Create semester mag flux tables
#sem05B = sem_mag_flux(sems[0])
#sem06B = sem_mag_flux(sems[1])
#sem07B = sem_mag_flux(sems[2])
#sem08B = sem_mag_flux(sems[3])
#sem09B = sem_mag_flux(sems[4])
#sem10B = sem_mag_flux(sems[5])
#sem11B = sem_mag_flux(sems[6])
#sem12B = sem_mag_flux(sems[7])
#
##%% Join these tables using their IDs (will be fine for just this part as numbers alwasy match)
#print('Join 1')
#sem06B.rename_column('NUMBER_06B', 'NUMBER_05B')
#semcom = join(sem05B, sem06B, keys='NUMBER_05B')
#print('Join 2')
#sem07B.rename_column('NUMBER_07B', 'NUMBER_05B')
#semcom = join(semcom, sem07B, keys='NUMBER_05B')
#print('Join 3')
#sem08B.rename_column('NUMBER_08B', 'NUMBER_05B')
#semcom = join(semcom, sem08B, keys='NUMBER_05B')
#print('Join 4')
#sem09B.rename_column('NUMBER_09B', 'NUMBER_05B')
#semcom = join(semcom, sem09B, keys='NUMBER_05B')
#print('Join 5')
#sem10B.rename_column('NUMBER_10B', 'NUMBER_05B')
#semcom = join(semcom, sem10B, keys='NUMBER_05B')
#print('Join 6')
#sem11B.rename_column('NUMBER_11B', 'NUMBER_05B')
#semcom = join(semcom, sem11B, keys='NUMBER_05B')
#print('Join 7')
#sem12B.rename_column('NUMBER_12B', 'NUMBER_05B')
#semcom = join(semcom, sem12B, keys='NUMBER_05B')

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

months = ['sep05','oct05','nov05','dec05', 'jan06', 'jan07']

sep05 = month_mag_flux(months[0])
oct05 = month_mag_flux(months[1])
nov05 = month_mag_flux(months[2])
dec05 = month_mag_flux(months[3])
jan06 = month_mag_flux(months[4])
jan07 = month_mag_flux(months[5])

#%% Join these tables using their IDs (will be fine for just this part as numbers alwasy match)
print('Join 1')
semcom = join(sep05, oct05, keys='NUMBER')
print('Join 2')
semcom = join(semcom, nov05, keys='NUMBER')
print('Join 3')
semcom = join(semcom, dec05, keys='NUMBER')
print('Join 4')
semcom = join(semcom, jan06, keys='NUMBER')
print('Join 5')
semcom = join(semcom, jan07, keys='NUMBER')
#print('Join 6')
#semcom = join(semcom, sem11B, keys='NUMBER')
#print('Join 7')
#semcom = join(semcom, sem12B, keys='NUMBER')

##%% Match these with various catalogs to create final tables
#from astropy.coordinates import match_coordinates_sky
#from astropy.coordinates import SkyCoord
#from astropy import units as u
#
## match with stars catalogue
#print('Matching Stars')
#stars = Table.read('DR8-secure-stars.fits')
#starscoord = SkyCoord(stars['RA']*u.degree, stars['DEC']*u.degree)
#semcomcoord = SkyCoord(semcom['ALPHA_J2000_05B'], semcom['DELTA_J2000_05B'])
#idx, d2d , _ = match_coordinates_sky(starscoord, semcomcoord)
#mask = d2d<=1*u.arcsec #make sure match is within 1 arcsec (like in topcat)
#idx = idx[mask]
#ind = np.arange(len(semcom)) #create array of indicies
#ind = np.delete(ind, idx) #remove those that relate to stars
#starsmf = semcom[idx] #create table of just stars
#semcomns = semcom[ind] #create table of no stars
#
##match with best catalogue
#print('Matching Best')
#best = Table.read('uds_multicat_ap3.v5b_best.fits')
#bestcoord = SkyCoord(best['ALPHA_J2000']*u.degree, best['DELTA_J2000']*u.degree)
#semcomnscoord = SkyCoord(semcomns['ALPHA_J2000_05B'], semcomns['DELTA_J2000_05B'])
#idx, d2d , _ = match_coordinates_sky(bestcoord, semcomnscoord)
#mask = d2d<=1*u.arcsec #make sure match is within 1 arcsec (like in topcat)
#idx = idx[mask]
#bestmf = semcomns[idx] #create best table with no stars
#
#bestmfcoord = SkyCoord(bestmf['ALPHA_J2000_05B'], bestmf['DELTA_J2000_05B'])
#
## match with xmm
#print('Matching XMM')
#xmm = Table.read('XMM_not_Chandra.fits')
#xmmcoord = SkyCoord(xmm['RAJ2000'], xmm['DEJ2000'])
#idx, d2d , _ = match_coordinates_sky(xmmcoord, bestmfcoord)
#mask = d2d<=5*u.arcsec #make sure match is within 5 arcsec (like in topcat)
#idx = idx[mask]
#xmmmf = bestmf[idx]
#
### match with chandra
#print('Matching Chandra')
#chan = Table.read('chandra_catalogue.fits')
#chan['RA'].unit = u.deg
#chan['Dec'].unit = u.deg
#chancoord = SkyCoord(chan['RA'], chan['Dec'])
#idx, d2d , _ = match_coordinates_sky(chancoord, bestmfcoord)
#mask = d2d<=1*u.arcsec #make sure match is within 1 arcsec (like in topcat)
#idx = idx[mask]
#chanmf = bestmf[idx]
#
## combine chandra and xmm
#print('Joining xray table')
#xraymf = vstack([chanmf, xmmmf])
##%%
## boolean whether a source is seen in x-rays
#xray = np.isin(bestmf['NUMBER_05B'], xraymf['NUMBER_05B'])
#xraycol = Column(xray, 'X-ray')
#bestmf.add_column(xraycol)
#
##%% Save the tables
#semcom.write('mag_flux_tables/month_mag_flux_table.fits')
##bestmf.write('mag_flux_tables/mag_flux_table_best.fits')
##starsmf.write('mag_flux_tables/stars_mag_flux_table.fits')
##xraymf.write('mag_flux_tables/xray_mag_flux_table_best.fits')