#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 29 10:10:13 2018

Code to check spatial seeing variations

@author: ppxee
"""

### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
#from astropy.table import Table #for handling tables
import numpy as np #for handling arrays
#import math
#from astropy.stats import median_absolute_deviation
import vari_funcs_no06 #my module to help run code neatly
plt.close('all') #close any open plots

sems = ['05B', '07B', '08B', '09B', '10B', '11B', '12B']
const = 3.727777965e-5 # constant that defines unit conversion between degrees and pixels

## for convolved
tbdata = fits.open('mag_flux_convS_best.fits')[1].data
sdata = fits.open('stars_mag_flux_convS.fits')[1].data
totfwhm1=np.array([])
totfwhm2=np.array([])
totfwhm3=np.array([])
totfwhm4=np.array([])
totf1=np.array([])
totf2=np.array([])
totf3=np.array([])
totf4=np.array([])

for sem in sems:
    ra = tbdata['ALPHA_J2000_'+sem]
    dec = tbdata['DELTA_J2000_'+sem]
    sra = sdata['ALPHA_J2000_'+sem]
    sdec = sdata['DELTA_J2000_'+sem]
    
    ### define bounds of quadrants ###
    maxra = np.max(ra)
    minra = np.min(ra)
    maxdec = np.max(dec)
    mindec = np.min(dec)
    
    midra = (maxra+minra)/2
    middec = (maxdec+mindec)/2
    
    ### create masks for quadrant ###
    mask1 = sra < midra
    mask2 = sdec < middec
    quad1data = sdata[mask1*mask2]
    
    mask1 = sra >= midra
    mask2 = sdec < middec
    quad2data = sdata[mask1*mask2]
    
    mask1 = sra < midra
    mask2 = sdec >= middec
    quad3data = sdata[mask1*mask2]
    
    mask1 = sra >= midra
    mask2 = sdec >= middec
    quad4data = sdata[mask1*mask2]
    
    ### find average fwhm in each quadrant - convert into arcsec###
    fr1 = np.median(quad1data['FWHM_WORLD_'+sem]) * 3600
    fr2 = np.median(quad2data['FWHM_WORLD_'+sem]) * 3600
    fr3 = np.median(quad3data['FWHM_WORLD_'+sem]) * 3600
    fr4 = np.median(quad4data['FWHM_WORLD_'+sem]) * 3600
    
    ### find average flux in each quadrant ###
    f1 = np.mean(quad1data['MAG_APER_5_'+sem])
    f2 = np.mean(quad2data['MAG_APER_5_'+sem])
    f3 = np.mean(quad3data['MAG_APER_5_'+sem])
    f4 = np.mean(quad4data['MAG_APER_5_'+sem])
#    
#    ### Plot FR ###
#    plt.figure()
#    plotras = np.array([np.mean([minra,midra]), np.mean([maxra,midra]), np.mean([minra,midra]),np.mean([maxra,midra])])
#    plotdecs = np.array([np.mean([mindec,middec]), np.mean([mindec,middec]), np.mean([maxdec,middec]), np.mean([maxdec,middec])])
#    plt.scatter(plotras, plotdecs, c=[fr1,fr2,fr3,fr4])
#    plt.ylim(mindec, maxdec)
#    plt.xlim(minra, maxra)
#    cbar = plt.colorbar()
#    plt.gca().invert_xaxis()
#    cbar.set_label('FWHM (arcsec)')
#    plt.ylabel('Dec')
#    plt.xlabel('RA')
#    plt.title(sem+' Convolved')
#    plt.savefig(sem+' Convolved')

    #add average to an array of averages
    totfwhm1 = np.append(totfwhm1, fr1)
    totfwhm2 = np.append(totfwhm2, fr2)
    totfwhm3 = np.append(totfwhm3, fr3)
    totfwhm4 = np.append(totfwhm4, fr4)
    #add average to an array of averages
    totf1 = np.append(totf1, f1)
    totf2 = np.append(totf2, f2)
    totf3 = np.append(totf3, f3)
    totf4 = np.append(totf4, f4)
 
#sdatafull = fits.open('stars_mag_flux_table.fits')[1].data
#sdata = fits.open('starsfwhm.fits')[1].data
#mask = sdatafull['MAG_APER_5_05B'] >= 12
#sdata = sdata[mask]
##totfwhm1=np.array([])
##totfwhm2=np.array([])
##totfwhm3=np.array([])
##totfwhm4=np.array([])
##
#for sem in sems:
#        
#    #for unconvolved
#    tbdata = fits.open('SE_outputs_yearstacks/'+sem+'_output.fits')[1].data
#    ra = tbdata['ALPHA_J2000']
#    dec = tbdata['DELTA_J2000']
#    sra = sdata['RA']
#    sdec = sdata['DEC']
#    
#    ### define bounds of quadrants ###
#    maxra = np.max(ra)
#    minra = np.min(ra)
#    maxdec = np.max(dec)
#    mindec = np.min(dec)
#    
#    midra = (maxra+minra)/2
#    middec = (maxdec+mindec)/2
#    
#    ### create masks for quadrant ###
#    mask1 = sra < midra
#    mask2 = sdec < middec
#    quad1data = sdata[mask1*mask2]
#    
#    mask1 = sra >= midra
#    mask2 = sdec < middec
#    quad2data = sdata[mask1*mask2]
#    
#    mask1 = sra < midra
#    mask2 = sdec >= middec
#    quad3data = sdata[mask1*mask2]
#    
#    mask1 = sra >= midra
#    mask2 = sdec >= middec
#    quad4data = sdata[mask1*mask2]
#    
#    fr1 = np.median(quad1data['FWHM_'+sem]) * 3600
#    fr2 = np.median(quad2data['FWHM_'+sem]) * 3600
#    fr3 = np.median(quad3data['FWHM_'+sem]) * 3600
#    fr4 = np.median(quad4data['FWHM_'+sem]) * 3600
##    ### find average flux in each quadrant ###
##    f1 = np.mean(quad1data['MAG_APER_5_'+sem])
##    f2 = np.mean(quad2data['MAG_APER_5_'+sem])
##    f3 = np.mean(quad3data['MAG_APER_5_'+sem])
##    f4 = np.mean(quad4data['MAG_APER_5_'+sem])
#    
#    ### Plot FR ###
##    plt.figure()
    plotras = np.array([np.mean([minra,midra]), np.mean([maxra,midra]), np.mean([minra,midra]),np.mean([maxra,midra])])
    plotdecs = np.array([np.mean([mindec,middec]), np.mean([mindec,middec]), np.mean([maxdec,middec]), np.mean([maxdec,middec])])
##    plt.scatter(plotras, plotdecs, c=[fr1,fr2,fr3,fr4])
##    plt.ylim(mindec, maxdec)
##    plt.xlim(minra, maxra)
##    cbar = plt.colorbar()
##    plt.gca().invert_xaxis()
##    cbar.set_label('FWHM (arcsec)')
##    plt.ylabel('Dec')
##    plt.xlabel('RA')
##    plt.title(sem+' Unconvolved')
##    plt.savefig(sem+' Unconvolved')
#    
#    #add average to an array of averages
#    totfwhm1 = np.append(totfwhm1, fr1)
#    totfwhm2 = np.append(totfwhm2, fr2)
#    totfwhm3 = np.append(totfwhm3, fr3)
#    totfwhm4 = np.append(totfwhm4, fr4)
###    add average to an array of averages
##    totf1 = np.append(totf1, f1)
##    totf2 = np.append(totf2, f2)
##    totf3 = np.append(totf3, f3)
##    totf4 = np.append(totf4, f4)
#    
#vari_funcs_no06.avg_lightcurve(totfwhm1)
#plt.title('Median FWHM in bottom right quadrant')
#plt.ylabel('FWHM (arcsec)')
#vari_funcs_no06.avg_lightcurve(totfwhm2)
#plt.title('Median FWHM in bottom left quadrant')
#plt.ylabel('FWHM (arcsec)')
#vari_funcs_no06.avg_lightcurve(totfwhm3)
#plt.title('Median FWHM in top right quadrant')
#plt.ylabel('FWHM (arcsec)')
#vari_funcs_no06.avg_lightcurve(totfwhm4)
#plt.title('Median FWHM in top left quadrant')
#plt.ylabel('FWHM (arcsec)')

#vari_funcs_no06.avg_lightcurve(totf1)
#plt.title('Median magnitude in bottom right quadrant')
#plt.ylabel('Magnitude')
#vari_funcs_no06.avg_lightcurve(totf2)
#plt.title('Median magnitude in bottom left quadrant')
#plt.ylabel('Magnitude')
#vari_funcs_no06.avg_lightcurve(totf3)
#plt.title('Median magnitude in top right quadrant')
#plt.ylabel('Magnitude')
#vari_funcs_no06.avg_lightcurve(totf4)
#plt.title('Median magnitude in top left quadrant')
#plt.ylabel('Magnitude')

avgfwhm = np.array([np.mean(totfwhm1), np.mean(totfwhm2), np.mean(totfwhm3), np.mean(totfwhm4)])
plt.figure()
plt.scatter(plotras, plotdecs, c=avgfwhm)
plt.ylim(mindec, maxdec)
plt.xlim(minra, maxra)
cbar = plt.colorbar()
plt.gca().invert_xaxis()
cbar.set_label('FWHM (arcsec)')
plt.ylabel('Dec')
plt.xlabel('RA')
plt.title('Convolved Average')
plt.savefig('plots/Spatial-Comparisons/Convolved')
