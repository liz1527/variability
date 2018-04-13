#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 31 09:55:25 2018

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
tbdata = fits.open('mag_flux_tables/mag_flux_convS_best.fits')[1].data
xdata = fits.open('mag_flux_tables/xray_mag_flux_convS_best.fits')[1].data
sdata = fits.open('mag_flux_tables/stars_mag_flux_convS.fits')[1].data


def quadrants(basedata, initdata):
    
    bra = basedata['ALPHA_J2000_05B']
    bdec = basedata['DELTA_J2000_05B']
    ira = initdata['ALPHA_J2000_05B']
    idec = initdata['DELTA_J2000_05B']

    ### define bounds of quadrants ###
    maxra = np.max(bra)
    minra = np.min(bra)
    maxdec = np.max(bdec)
    mindec = np.min(bdec)
    
    midra = (maxra+minra)/2
    middec = (maxdec+mindec)/2
    
    ### create masks for quadrant ###
    mask1 = ira < midra
    mask2 = idec < middec
    quad1data = initdata[mask1*mask2]
    
    mask1 = ira >= midra
    mask2 = idec < middec
    quad2data = initdata[mask1*mask2]
    
    mask1 = ira < midra
    mask2 = idec >= middec
    quad3data = initdata[mask1*mask2]
    
    mask1 = ira >= midra
    mask2 = idec >= middec
    quad4data = initdata[mask1*mask2]
    
    return quad1data, quad2data, quad3data, quad4data

### split data into 4 quadrants ###
quad1data, quad2data, quad3data, quad4data = quadrants(tbdata, tbdata)
xquad1data, xquad2data, xquad3data, xquad4data = quadrants(tbdata, xdata)
squad1data, squad2data, squad3data, squad4data = quadrants(tbdata, sdata)

def magvaryplot(quaddata, xquaddata, squaddata):
    # get mag stacks
    flux = vari_funcs_no06.mag5_stacks(quaddata)
    xflux = vari_funcs_no06.mag5_stacks(xquaddata)
    sflux = vari_funcs_no06.mag5_stacks(squaddata)
    
    # remove 99s
    flux, quaddata = vari_funcs_no06.no99(flux, quaddata)
    xflux, xquaddata = vari_funcs_no06.no99(xflux, xquaddata)
    sflux, squaddata = vari_funcs_no06.no99(sflux, squaddata)
    
    avg = np.median(flux, axis=0)
    vari_funcs_no06.avg_lightcurve(avg)
    
    # plot
    vari_funcs_no06.flux_variability_plot(flux, xflux, 'mad', starflux=sflux,
                                          stars=True)
    # Find outliers ###
    bins = np.array([13, 15])
    bins = np.append(bins, np.arange(16,24,0.2))
    bins = np.append(bins, [24, 25, 26])
    outliers, tbnew, allmodz = vari_funcs_no06.find_outliers(flux, quaddata, bins)

    return outliers, tbnew, allmodz
    
quad1outliers, quad1datanew, quad1allmodz = magvaryplot(quad1data, xquad1data, squad1data)
plt.ylim(2e-4, 3)
plt.title('Bottom Right Quadrant')
plt.savefig('plots/MAD-Magnitude/magmadBRconv')
# put correct title on light curve
plt.figure(1)
plt.title('Median Magnitude Curve fro Bottom Right Quadrant')
plt.savefig('plots/Lightcurves/fluxcurveBRconv')
quad2outliers, quad2datanew, quad2allmodz = magvaryplot(quad2data, xquad2data, squad2data)
plt.ylim(2e-4, 3)
plt.title('Bottom Left Quadrant')
plt.savefig('plots/MAD-Magnitude/magmadBLconv')
plt.figure(3)
plt.title('Median Magnitude Curve for Bottom Left Quadrant')
plt.savefig('plots/Lightcurves/fluxcurveBL')
quad3outliers, quad3datanew, quad3allmodz = magvaryplot(quad3data, xquad3data, squad3data)
plt.ylim(2e-4, 3)
plt.title('Top Right Quadrant')
plt.savefig('plots/MAD-Magnitude/magmadTRconv')
plt.figure(5)
plt.title('Median Magnitude Curve for Top Right Quadrant')
plt.savefig('plots/Lightcurves/fluxcurveTRconv')
quad4outliers, quad4datanew, quad4allmodz = magvaryplot(quad4data, xquad4data, squad4data)
plt.ylim(2e-4, 3)
plt.title('Top Left Quadrant')
plt.savefig('plots/MAD-Magnitude/magmadTLconv')
plt.figure(7)
plt.title('Median Magnitude Curve for Top Left Quadrant')
plt.savefig('plots/Lightcurves/fluxcurveTLconv')



#outliers, tbnew, allmodz = vari_funcs_no06.find_outliers()
    
    
    
    
    
    
    
    
    
    
    