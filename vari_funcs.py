#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 10:32:41 2017

Code containing all functions that will be useful in codes looking at 
variability and plotting light curves.

This should be imported at the start of the majority of my codes.

Updated on 27/2/18 to change the way the mag_flux tables are read and to bring 
in line with updates done on vari_funcs_no06 while 06 was not being used.

@author: ppxee
"""

### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
#from astropy.table import Table #for handling tables
import numpy as np #for handling arrays
import math
from astropy.stats import median_absolute_deviation
from scipy.stats import chisquare
font = {'family' : 'DejaVu Sans',
        'weight' : 'normal',
        'size'   : 14}

plt.rc('font', **font)

### Define Functions ###
def flux4_stacks(tbdata):
    ''' Function that takes a catalogue of flux data from the sextracor output
    and makes a np array containing only the 3 arcsec aperture data for each 
    epoch
    Input:
        tbdata = the original combined catalogue of flux data 
    Output:
        flux = an array with 8 columns containing flux values for each year '''
        
    flux = np.stack(([tbdata['FLUX_APER_05B'][:,3], tbdata['FLUX_APER_06B'][:,3],
                tbdata['FLUX_APER_07B'][:,3], tbdata['FLUX_APER_08B'][:,3],
                tbdata['FLUX_APER_09B'][:,3], tbdata['FLUX_APER_10B'][:,3], 
                tbdata['FLUX_APER_11B'][:,3], tbdata['FLUX_APER_12B'][:,3]]), axis=1)
    return flux

def fluxerr4_stacks(tbdata):
    ''' Function that takes a catalogue of flux data from the sextracor output
    and makes a np array containing only the 3 arcsec aperture error data for 
    each epoch
    Input:
        tbdata = the original combined catalogue of flux data 
    Output:
        fluxerr = an array with 8 columns containing flux error values for each
        year '''
        
    fluxerr = np.stack([tbdata['FLUXERR_APER_05B'][:,3],tbdata['FLUXERR_APER_06B'][:,3],
                tbdata['FLUXERR_APER_07B'][:,3],tbdata['FLUXERR_APER_08B'][:,3],
                tbdata['FLUXERR_APER_09B'][:,3],tbdata['FLUXERR_APER_10B'][:,3], 
                tbdata['FLUXERR_APER_11B'][:,3],tbdata['FLUXERR_APER_12B'][:,3]], axis=1)
    return fluxerr

def flux5_stacks(tbdata):
    ''' Function that takes a catalogue of flux data from the sextracor output
    and makes a np array containing only the 2 arcsec aperture data for each 
    epoch
    Input:
        tbdata = the original combined catalogue of flux data 
    Output:
        flux = an array with 8 columns containing flux values for each year '''
        
    flux = np.stack(([tbdata['FLUX_APER_05B'][:,4], tbdata['FLUX_APER_06B'][:,4],
                tbdata['FLUX_APER_07B'][:,4], tbdata['FLUX_APER_08B'][:,4],
                tbdata['FLUX_APER_09B'][:,4], tbdata['FLUX_APER_10B'][:,4], 
                tbdata['FLUX_APER_11B'][:,4], tbdata['FLUX_APER_12B'][:,4]]), axis=1)
    return flux

def fluxerr5_stacks(tbdata):
    ''' Function that takes a catalogue of flux data from the sextracor output
    and makes a np array containing only the 2 arcsec aperture error data for 
    each epoch
    Input:
        tbdata = the original combined catalogue of flux data 
    Output:
        fluxerr = an array with 8 columns containing flux error values for each
        year '''
        
    fluxerr = np.stack([tbdata['FLUXERR_APER_05B'][:,3],tbdata['FLUXERR_APER_06B'][:,3],
                tbdata['FLUXERR_APER_07B'][:,3],tbdata['FLUXERR_APER_08B'][:,3],
                tbdata['FLUXERR_APER_09B'][:,3],tbdata['FLUXERR_APER_10B'][:,3], 
                tbdata['FLUXERR_APER_11B'][:,3],tbdata['FLUXERR_APER_12B'][:,3]], axis=1)
    return fluxerr
   
def mag4_stacks(tbdata):
    ''' Function that takes a catalogue of magnitude data from the sextracor 
    output and makes a np array containing only the 3 arcsec aperture data for 
    each epoch
    Input:
        tbdata = the original combined catalogue of flux data 
    Output:
        flux = an array with 8 columns containing flux values for each year '''
        
    flux = np.stack(([tbdata['MAG_APER_05B'][:,3], tbdata['MAG_APER_06B'][:,3],
                tbdata['MAG_APER_07B'][:,3], tbdata['MAG_APER_08B'][:,3],
                tbdata['MAG_APER_09B'][:,3], tbdata['MAG_APER_10B'][:,3], 
                tbdata['MAG_APER_11B'][:,3], tbdata['MAG_APER_12B'][:,3]]), axis=1)
    return flux

def magerr4_stacks(tbdata):
    ''' Function that takes a catalogue of magnitude data from the sextracor 
    output and makes a np array containing only the 3 arcsec aperture error 
    data for each epoch
    Input:
        tbdata = the original combined catalogue of flux data 
    Output:
        fluxerr = an array with 8 columns containing flux error values for each
        year '''
        
    fluxerr = np.stack([tbdata['MAGERR_APER_05B'][:,5],tbdata['MAGERR_APER_06B'][:,5],
                tbdata['MAGERR_APER_07B'][:,5],tbdata['MAGERR_APER_08B'][:,5],
                tbdata['MAGERR_APER_09B'][:,5],tbdata['MAGERR_APER_10B'][:,5], 
                tbdata['MAGERR_APER_11B'][:,5],tbdata['MAGERR_APER_12B'][:,5]], axis=1)
    return fluxerr

def mag5_stacks(tbdata):
    ''' Function that takes a catalogue of magnitude data from the sextractor 
    output and makes a np array containing only the 2 arcsec aperture data for
    each epoch
    Input:
        tbdata = the original combined catalogue of flux data 
    Output:
        flux = an array with 8 columns containing flux values for each year '''
        
    flux = np.stack(([tbdata['MAG_APER_05B'][:,4], tbdata['MAG_APER_06B'][:,4],
                tbdata['MAG_APER_07B'][:,4], tbdata['MAG_APER_08B'][:,4],
                tbdata['MAG_APER_09B'][:,4], tbdata['MAG_APER_10B'][:,4], 
                tbdata['MAG_APER_11B'][:,4], tbdata['MAG_APER_12B'][:,4]]), axis=1)
    return flux

def magerr5_stacks(tbdata):
    ''' Function that takes a catalogue of magnitude data from the sextractor 
    output and makes a np array containing only the 2 arcsec aperture error 
    data for each epoch
    Input:
        tbdata = the original combined catalogue of flux data 
    Output:
        fluxerr = an array with 8 columns containing flux error values for each
        year '''
        
    fluxerr = np.stack([tbdata['MAGERR_APER_05B'][:,4],tbdata['MAGERR_APER_06B'][:,4],
                tbdata['MAGERR_APER_07B'][:,4],tbdata['MAGERR_APER_08B'][:,4],
                tbdata['MAGERR_APER_09B'][:,4],tbdata['MAGERR_APER_10B'][:,4], 
                tbdata['MAGERR_APER_11B'][:,4],tbdata['MAGERR_APER_12B'][:,4]], axis=1)
    return fluxerr 

def lightcurve4(ob, fitsdata)  :
    ''' Function that plots the light curve of an object in terms of its flux 
    in an aperture 4 pixels across (i.e. 2 arcsec in diameter) 
    Inputs:
        ob = the ID of the object that you want the lightcurve from
        fitsdata = the original catalogue of data that the curve will be 
                    plotted from 
    Output:
        None '''
        
    #Get data for the object called
    mask = fitsdata['NUMBER_05B'] == ob
    obdata = fitsdata[mask]
    if not obdata: #Reject if no object number matches the input value
        print('error- invalid object number')
        return
    
    #Create arrays of flux values and error values
    flux = np.array([obdata['FLUX_APER_05B'][:,4],obdata['FLUX_APER_06B'][:,4],
                    obdata['FLUX_APER_07B'][:,4],obdata['FLUX_APER_08B'][:,4],
                    obdata['FLUX_APER_09B'][:,4],obdata['FLUX_APER_10B'][:,4], 
                    obdata['FLUX_APER_11B'][:,4],obdata['FLUX_APER_12B'][:,4]])
    
    fluxerr = np.array([obdata['FLUXERR_APER_05B'][:,4],obdata['FLUXERR_APER_06B'][:,4],
                    obdata['FLUXERR_APER_07B'][:,4],obdata['FLUXERR_APER_08B'][:,4],
                    obdata['FLUXERR_APER_09B'][:,4],obdata['FLUXERR_APER_10B'][:,4], 
                    obdata['FLUXERR_APER_11B'][:,4],obdata['FLUXERR_APER_12B'][:,4]])
    
    # normalise and correct for seeing
    flux = psf_correct(flux, flux, 'mean')
    fluxerr = psf_correct(flux, fluxerr, 'mean')
    avgflux =np.mean(flux)
    normflux = flux/avgflux
    normerr = fluxerr/avgflux
    
    #set up time variable for plot
    t = np.linspace(1, 8, num=8)
    years = ['05B', '06B', '07B', '08B', '09B', '10B', '11B', '12B']
    
    #Plot graph in new figure
    plt.figure()
    plt.xticks(t, years)
    plt.errorbar(t, normflux, yerr=normerr, fmt = 'ro')
    plt.xlabel('Semester')
    plt.ylabel('Flux of object')
    plt.title('Light curve for object number %i' % ob)
    return

def lightcurve5(ob, fitsdata)  :
    ''' Function that plots the light curve of an object in terms of its flux 
    in an aperture 5 pixels across (i.e. 3 arcsec in diameter) 
    Inputs:
        ob = the ID of the object that you want the lightcurve from
        fitsdata = the original catalogue of data that the curve will be 
                    plotted from 
    Output:
        None '''
        
    #Get data for the object called
    mask = fitsdata['NUMBER_05B'] == ob
    obdata = fitsdata[mask]
    if not obdata: #Reject if no object number matches the input value
        print('error- invalid object number')
        return
    
    #Create arrays of flux values and error values
    flux = np.array([obdata['MAG_APER_05B'][:,4],obdata['MAG_APER_06B'][:,4],
                    obdata['MAG_APER_07B'][:,4],obdata['MAG_APER_08B'][:,4],
                    obdata['MAG_APER_09B'][:,4],obdata['MAG_APER_10B'][:,4], 
                    obdata['MAG_APER_11B'][:,4],obdata['MAG_APER_12B'][:,4]])
    
    fluxerr = np.array([obdata['MAGERR_APER_05B'][:,4],obdata['MAGERR_APER_06B'][:,4],
                    obdata['MAGERR_APER_07B'][:,4],obdata['MAGERR_APER_08B'][:,4],
                    obdata['MAGERR_APER_09B'][:,4],obdata['MAGERR_APER_10B'][:,4], 
                    obdata['MAGERR_APER_11B'][:,4],obdata['MAGERR_APER_12B'][:,4]])
    
    #create full arrays for zeropoint correction
    fullflux = mag5_stacks(fitsdata)
    fullflux, fitsdata = no99(fullflux, fitsdata)
    
    # normalise and correct for seeing
    fullflux = normalise_mag(fullflux)
    avgfluxperepoch = np.median(fullflux, axis=0)#for UDS
    avgflux = 1#np.median(fullflux)
    const = avgflux-avgfluxperepoch
    fluxcorr = flux + const[:,None]
#    fluxcorr = psf_correct_mag(fullflux, flux, 'median')
    fluxerrcorr = err_correct(flux, fluxerr, fluxcorr)
#    avgflux =np.mean(fluxcorr)
#    normflux = fluxcorr/avgflux
#    normerr = fluxerrcorr/avgflux
    
    #set up time variable for plot
    t = np.linspace(1, 8, num=8)
    years = ['05B', '06B', '07B', '08B', '09B', '10B', '11B', '12B']
    
    #Plot graph in new figure
    plt.figure()
    plt.xticks(t, years)
    plt.errorbar(t, fluxcorr, yerr=fluxerrcorr, fmt = 'ro')
#    plt.errorbar(t, flux, yerr=fluxerr, fmt = 'ro')
    plt.xlabel('Semester')
    plt.ylabel('K-band magnitude of object')
    plt.title('Light curve for object number %i' % ob)
    return

def chandra_only(tbdata):
    ''' Function that restricts the objects included in analysis to only 
    those within the Chandra footprint 
    Input:
        tbdata = original catalogue of data 
    Output:
        newtbdata = new catalogue of data which only includes objects within 
                    the chandra footprint '''
    ### Restrict objects to those in the Chandra field ###
    mask1 = tbdata['DELTA_J2000_05B'] < -4.93 #max Dec
    mask2 = tbdata['DELTA_J2000_05B'] > -5.403 #min Dec
    mask3 = tbdata['ALPHA_J2000_05B'] < 34.72 #max RA
    mask4 = tbdata['ALPHA_J2000_05B'] > 34.07 #min RA
    mask = mask1 * mask2 * mask3 * mask4
    newtbdata = tbdata[mask]
    return(newtbdata)
    
#Calculate varience for each row
def sigmasq(flux, baseerr):
    ''' Function that calculates the excess varience value for each row in an 
    array 
    Inputs:
        flux = array of fluxes from objects in a number of epochs 
        baseerr = array of errors that the mean error should be calculated from
    Output:
        sig = array of excess variance values for every object '''
    avgflux = np.nanmean(flux, axis=1)
#    meanerr = np.nanmean(baseerr, axis=1)
#    meanerrsq = np.square(meanerr)
    N = np.size(flux, axis=1)
    numobs = np.size(flux, axis=0)
    sig = [((flux[n, :]- avgflux[n])**2 - (baseerr[n,:])**2) for n in range(numobs)]# 
    sigsum = np.nansum(sig, axis=1)
    normsig = sigsum/(N)
    return normsig

def normsigmasq(flux, baseerr):
    ''' Function that calculates the excess varience value for each row in an 
    array 
    Inputs:
        flux = array of fluxes from objects in a number of epochs 
        baseerr = array of errors that the mean error should be calculated from
    Output:
        sig = array of excess variance values for every object '''
    avgflux = np.nanmean(flux, axis=1)
    N = np.size(flux, axis=1)
    numobs = np.size(flux, axis=0)
    sig = [((flux[n, :]- avgflux[n])**2 - (baseerr[n,:])**2) for n in range(numobs)]# 
    sigsum = np.nansum(sig, axis=1)
    normsig = sigsum/(N)
    return normsig

def fluxbin(min, max, flux, tbdata):
    ''' Separate a flux array into a bin depending on the average flux of the
    object 
    Inputs:
        min = the minimum flux value for the bin
        max = the maximum flux value for the bin
        flux = the array for fluxes for the object
    Output:
        fluxbin = the array of fluxes for objects whose average flux is within 
                    the limits of the bin '''
    avgflux = np.mean(flux, axis=1)
    bina = avgflux >= min
    binb = avgflux < max
    bin = bina*binb
    fluxbin = flux[bin]
    tbdata = tbdata[bin]
    return fluxbin, tbdata

def avg_lightcurve(avgfluxarr, errarr):
    ''' Plot light curves for a provided set of fluxes rather than extracting
    a certain objects fluxes from a larger array.
    Input:
        avgfluxarr = array of 8 semester flux values
    Output:
        ax = axes handle so the plot can be added to if necessary '''
        
    #Check the array has values
    if math.isnan(avgfluxarr[0]):
        print('Array is empty')
        return
    #set up time variable for plot
    t = np.linspace(1, 8, num=8)
    years = ['05B', '06B', '07B', '08B', '09B', '10B', '11B', '12B']
    
    #Plot graph in new figure
    plt.figure()
    ax = plt.axes()
    plt.xticks(t, years)
    ax.errorbar(t, avgfluxarr, errarr, fmt='ro')
    plt.xlabel('Semester')
    plt.ylabel('Flux of object')
    plt.title('Average light curve')
#    plt.ylim(ymin = 6.3, ymax=7.2)
    return ax

def normalise_flux(flux):
    ''' Normalise each objects flux to its average value
    Input:
        flux = array of object flux values 
    Output:
        array of object flux values normalised to the average flux of the 
        object '''
    avgflux = np.mean(flux, axis=1)
    return flux / avgflux[:,None]

def normalise_flux_and_errors(flux, fluxerr):
    ''' Normalise each objects flux to its average value
    Input:
        flux = array of object flux values 
    Output:
        array of object flux values normalised to the average flux of the 
        object '''
    avgflux = np.mean(flux, axis=1)
    flux = flux/avgflux[:,None]
    fluxerr = fluxerr/avgflux[:,None]
    return flux, fluxerr

def normalise_mag(mag):
    ''' Normalise each objects flux to its average value
    Input:
        flux = array of object flux values 
    Output:
        array of object flux values normalised to the average flux of the 
        object '''
    avgflux = np.mean(mag, axis=1)
    diff = avgflux - 1
    return mag - diff[:,None]

def no99(fluxn, tbdata):
    fluxn[fluxn == 99] = np.nan
    mask = ~np.isnan(fluxn).any(axis=1)
    fluxn = fluxn[mask]
    tbdata = tbdata[mask]
    return fluxn, tbdata

def noneg(fluxn, tbdata):
    fluxn[fluxn <= 0] = np.nan
    mask = ~np.isnan(fluxn).any(axis=1)
    fluxn = fluxn[mask]
    tbdata = tbdata[mask]
    return fluxn, tbdata

def onpickchanonly(event):
    ''' Function that plots the lightcurve of an object when it is clicked on 
    the vairiability v flux plot '''
        
    combined = fits.open('mag_flux_tables/mag_flux_table_best.fits')
    tbdata = combined[1].data
    tbdata = chandra_only(tbdata)
    
    ### remove values that are +/-99 ###
    fluxn = mag5_stacks(tbdata)
    fluxn[fluxn == 99] = np.nan
    mask = ~np.isnan(fluxn).any(axis=1)
    tbdata = tbdata[mask]
    
    ob = tbdata['NUMBER_05B'][event.ind] #Define the object number from the index of the selected point
    if len(ob) > 1: #Reject selection if more than one object has been selected
        print('Too many objects selected')
        return
    print('Object identified')
    lightcurve5(ob, tbdata) #Plot the lightcurve from the lightcurve function

def onpick(event):
    ''' Function that plots the lightcurve of an object when it is clicked on 
    the vairiability v flux plot '''
        
    combined = fits.open('mag_flux_tables/mag_flux_table_best.fits')
    tbdata = combined[1].data
    
    ### remove values that are +/-99 ###
    fluxn = mag5_stacks(tbdata)
    fluxn, tbdata = no99(fluxn, tbdata)
    
    ob = tbdata['NUMBER_05B'][event.ind] #Define the object number from the index of the selected point
    if len(ob) > 1: #Reject selection if more than one object has been selected
        print('Too many objects selected')
        return
    print('Object identified')
    lightcurve5(ob, tbdata) #Plot the lightcurve from the lightcurve function
    
def flux_variability_plot(flux, fluxchan, plottype, flux2 = [], fluxchan2 = [],
                          fluxerr = [], fluxerr2 = [], starflux=[], starfluxerr=[],
                          comparison = False, normalised = False, stars=False,
                          psfcorrect=False, chanerr = [], chanerr2 = []):
    ''' Function to plot the variability vs mean flux plot using either the
    MAD statistic or the Excess Variance statistic 
    If using MAD, specify plottype = 'mad' and supply just the flux and chandra 
    flux arrays
    If using Excess Variance, specify plottype = 'excess' and supply the error
    arrays as well as the flus ones
    If you with to have a comparison plot of either type make sure to specify 
    the greyed out one first and then the coloured one 
    Inputs:
        flux = array of flux values for UDS objects
        fluxchan = array of flux values for chandra objects
        plottype = either 'mad' or 'excess' which defines what statistic is
                    used to quantify variability
        flux2 = optional second array of UDS objects if comparing two sets 
                (e.g. a before and after) this array will be the top coloured
                layer.
        fluxchan2 = optional second array of chandra objects 
        fluxerr = optional array of flux errors for objects, used if 
                    'excess' is specified
        fluxerr2 = optional second array of flux errors for objects, used 
                    if 'excess' and 'comparision = True' are specified
        starflux = optional array of fluxes for stars
        starfluxerr = optional array of flux errors for objects, used if 
                        'excess' and 'stars=True' are specified
        comparison = True or False (default) depending on if a comparison plot
                        is required 
        normalised = True or False (default) depending if the fluxes should be
                        normalised to the objects average flux 
        stars = True or False (default) depending on if stars should be added
                to the plot
    Output:
        fig = figure handle to allow clicking for light curves to be enabled if
                required '''
    
    fig = plt.figure()
    avgfluxperob = np.mean(flux, axis=1) #for UDS
    avgfluxchanperob = np.mean(fluxchan, axis=1) #for non-stellar chandra
    if stars == True:
        savgfluxperob = np.mean(starflux, axis=1) #for stars

    ### Check if normalisation is true and normalise if necessary ###
    if normalised == True:
        fluxold = flux # need for error calc
        chanold = fluxchan
#        flux = normalise_mag(flux)
#        fluxchan = normalise_mag(fluxchan) 
#        if stars == True:
#            starflux = normalise_mag(starflux)
        flux = normalise_flux(flux)
        fluxchan = normalise_flux(fluxchan) 
        if stars == True:
            starold = starflux
            starflux = normalise_flux(starflux)
    if psfcorrect == True:
        fluxchan = psf_correct_mag(flux, fluxchan, 'median')
        if stars == True:
            starflux = psf_correct_mag(flux, starflux, 'median')
        flux = psf_correct_mag(flux, flux, 'median')
    ### Find out which plot type is specified and calculate appropriate statistic ###
    if plottype == 'mad':
        vary = median_absolute_deviation(flux, axis=1)
        varychan = median_absolute_deviation(fluxchan, axis=1)
        if stars == True:
            varystar = median_absolute_deviation(starflux, axis=1)
        plt.ylabel('MAD')
    elif plottype == 'excess':
        if normalised == True:
             #need to normalise the errors as well as the flux values
            fluxerr = err_correct_flux(fluxold, fluxerr)
            chanerr = err_correct_flux(chanold, chanerr)
            vary = normsigmasq(flux, fluxerr)
            varychan = normsigmasq(fluxchan, chanerr)
        else:
            vary = sigmasq(flux, fluxerr)
            varychan = sigmasq(fluxchan, chanerr)
        if stars == True:
            if normalised == True:
                starfluxerr = err_correct_flux(starold, starfluxerr)
                varystar = normsigmasq(starflux, starfluxerr)
            else:
                varystar = sigmasq(starflux, starfluxerr)
        plt.ylabel('Excess Variance')
    elif plottype == 'chisq':
        [vary, _] = chisquare(flux, axis=1)
        [varychan, _] = chisquare(fluxchan, axis=1)
        if stars == True:
            [varystar, _] = chisquare(starflux, axis=1)
    else:
        print('Invalid plottype') #returns if unrecognised value is entered
        return
    
    ### Plot the variability v mean as appropriate ###
    if comparison == True:     
        
        avgfluxperob2 = np.mean(flux2, axis=1) #for UDS
        avgfluxchanperob2 = np.mean(fluxchan2, axis=1) #for non-stellar chandra
        
        if normalised == True:
            flux2old = flux # need for error calc
            flux2 = normalise_mag(flux2)
            fluxchan2 = normalise_mag(fluxchan2)     
        if plottype == 'mad':
            varycorr = median_absolute_deviation(flux2, axis=1)
            varychancorr = median_absolute_deviation(fluxchan2, axis=1)
        elif plottype == 'excess':
#            if normalised == True:
#                # need to normalise the errors as well as the flux values
#                fluxerr2 = err_correct(flux2old, fluxerr2, flux2)
            varycorr = normsigmasq(flux2, fluxerr2)
            varychancorr = normsigmasq(fluxchan2, chanerr2)

        ### plot varibility v flux graph for original in gray ###
        plt.plot(avgfluxperob, vary, '+', color =  'tab:gray', label='UDS before correction', alpha = 0.5) 
        plt.plot(avgfluxchanperob, varychan, 'o', color =  'tab:gray',  mfc = 'none', markersize = 10,
                 label='X-ray Source', alpha = 0.5) #no picker as will be selected in the UDS point
        
        ### plot varibility v flux graph for new in colour ###
        line, = plt.plot(avgfluxperob2, varycorr, 'b+', label='UDS Source', picker=2, alpha = 0.5) #label tells 
        # the legend what to write, picker is the readius of pixels which a user must 
        # click within to call the object
        plt.plot(avgfluxchanperob2, varychancorr, 'ro', mfc = 'none', markersize = 10,
                 label='X-ray Source', alpha = 0.5) #no picker as will be selected in the UDS point
    else:
        if stars==True:
            plt.plot(savgfluxperob, varystar, 'm*', mfc = 'none', markersize = 10,
                     label='Secure Star') 
        line, = plt.plot(avgfluxperob, vary, 'b+', label='UDS Source', picker=2)
        plt.plot(avgfluxchanperob, varychan, 'ro', mfc = 'none', markersize = 10,
                 label='X-ray Source') #no picker as will be selected in the UDS point

        
    ### Apply required plot charateristics ###
    plt.xscale('log')
#    plt.yscale('symlog', linthreshy=0.001)
#    plt.yscale('log')
#    plt.ylim(2e-4, 3)
#    plt.xlim(13,26)
    plt.xlabel('Mean Flux')
    plt.legend()
#    plt.gca().invert_xaxis()
    return fig, vary


def psf_correct_flux(baseflux, initflux, avgtype):
    ''' Function that applies a fractional PSF correction to initflux, based on
    the constant required to correct an epochs average flux to the overall 
    average flux in the baseflux array.
    i.e. if doing a correction based on all objects for objects in the chandra
    data then baseflux = flux and initflux = fluxchan.
    If doing a correction based on the stellar fluxes for all objects in the 
    field then baseflux = sflux and initflux = flux. 
    Basetype must be defined as either 'all' or 'star' so that the mean value
    can be used for all objects and the median value can be used for stellar 
    objects 
    Inputs:
        baseflux = flux array to base the corrections on
        initflux = flux array to be corrected
        basetype = either 'mean' or 'median' which dictate if the mean or 
                    median value is used for the correction 
    Output:
        Flux array with values crudely corrected for differences in seeing
        (average flux should now be the same for each epoch). '''
    if avgtype == 'mean':
        avgfluxperepoch = np.mean(baseflux, axis=0)#for UDS
        avgflux = np.mean(baseflux)
        const = avgflux/avgfluxperepoch
    elif avgtype == 'median':
        avgfluxperepoch = np.median(baseflux, axis=0)#for UDS
        avgflux = np.median(baseflux)
        const = avgflux/avgfluxperepoch
    else:
        print('Invalid basetype')
        return
    return initflux * const[None,:]

def psf_correct_mag(basemag, initmag, avgtype):
    ''' Function that applies a fractional PSF correction to initflux, based on
    the constant required to correct an epochs average flux to the overall 
    average flux in the baseflux array.
    i.e. if doing a correction based on all objects for objects in the chandra
    data then baseflux = flux and initflux = fluxchan.
    If doing a correction based on the stellar fluxes for all objects in the 
    field then baseflux = sflux and initflux = flux. 
    Basetype must be defined as either 'all' or 'star' so that the mean value
    can be used for all objects and the median value can be used for stellar 
    objects 
    Inputs:
        baseflux = flux array to base the corrections on
        initflux = flux array to be corrected
        basetype = either 'mean' or 'median' which dictate if the mean or 
                    median value is used for the correction 
    Output:
        Flux array with values crudely corrected for differences in seeing
        (average flux should now be the same for each epoch). '''
    if avgtype == 'mean':
        avgfluxperepoch = np.mean(basemag, axis=0)#for UDS
        avgflux = np.mean(basemag)
        const = avgflux-avgfluxperepoch
    elif avgtype == 'median':
        avgfluxperepoch = np.median(basemag, axis=0)#for UDS
        avgflux = np.median(basemag)
        const = avgflux-avgfluxperepoch
    else:
        print('Invalid basetype')
        return
    return initmag + const[None,:]

def err_correct(flux, fluxerr, fluxnew):
    ''' Function that applies a correction to the array of error values that 
    matches the correction applied to the corresponding array of fluxes.
    Inputs:
        flux = initial flux array before any corrections were applied
        fluxerr = initial flux err array
        fluxcorr = array of fluxes that have been corrected
    Output:
        Flux error array with values crudely corrected '''
        
    return fluxnew * (fluxerr/flux)

def err_correct_flux(oldflux, fluxerr):
    ''' Function that applies a correction to the array of error values that 
    matches the correction applied to the corresponding array of fluxes.
    Inputs:
        flux = initial flux array before any corrections were applied
        fluxerr = initial flux err array
        fluxcorr = array of fluxes that have been corrected
    Output:
        Flux error array with values crudely corrected '''
    avgflux = np.mean(oldflux, axis=1)
    return fluxerr/avgflux[:,None]

def mod_z_score(arr):
    medx = np.median(arr)
    mad = median_absolute_deviation(arr)
    zvalues = np.array([(0.6745*(x-medx))/mad for x in arr])
    return zvalues

def find_outliers(flux, tbdata, bins, threshold=6):
    ### Bin data ###
    allmodz = []
    tbnew = np.recarray([0], dtype=tbdata.dtype, formats=tbdata.formats)
    for n, binedge in enumerate(bins):
        if n==np.size(bins)-1:
            break
        fluxbin1, tbbin1 = fluxbin(binedge, bins[n+1], flux, tbdata)
        #calulate mad values in the bins
        mad1 = median_absolute_deviation(fluxbin1, axis=1) #for UDS
        modz = mod_z_score(mad1)
        tbnew = np.hstack([tbnew, tbbin1])
        allmodz = np.append(allmodz, modz)
    outliers = allmodz>=threshold
    return outliers, tbnew, allmodz



