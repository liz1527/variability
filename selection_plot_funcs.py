#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 15:52:36 2019

Module to contain any functions relating to making the main interactive 
selection plot.

@author: ppxee
"""

### Import general libraries ###
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
from astropy.table import Table #for handling tables
import numpy as np #for handling arrays
from astropy.stats import median_absolute_deviation

### Import my modules ###
import field_funcs #for restricting field
import k_mag_flux #for creating lightcurves in Ks
import j_mag_flux #for creating lightcurves in J
import flux_funcs #for removing bad flux values
import lightcurve_funcs #for plotting lightcurves
import correction_funcs #for crude psf corrections
import vary_stats #for my variability functions

def onpickchanonly(event):
    ''' Function that plots the lightcurve of an object when it is clicked on 
    the vairiability v flux plot.
    This version is for plotting a 3 arcsec magnitude light curve when plot is 
    restricted to just the chandra region.
    Inputs:
        event = info from the click
    '''
        
    combined = fits.open('mag_flux_tables/mag_flux_table_best.fits')
    tbdata = combined[1].data
    tbdata = field_funcs.chandra_only(tbdata)
    
    ### remove values that are +/-99 ###
    fluxn = k_mag_flux.mag_stacks(tbdata, aper=5)
    fluxn[fluxn == 99] = np.nan
    mask = ~np.isnan(fluxn).any(axis=1)
    tbdata = tbdata[mask]
    
    ob = tbdata['NUMBER_05B'][event.ind] #Define the object number from the index of the selected point
    if len(ob) > 1: #Reject selection if more than one object has been selected
        print('Too many objects selected')
        return
    print('Object identified')
    lightcurve_funcs.lightcurve5(ob, tbdata) #Plot the lightcurve from the lightcurve function

def onpick(event):
    ''' Function that plots the lightcurve of an object when it is clicked on 
    the vairiability v flux plot
    This version is for plotting a 3 arcsec magnitude light curve.
    Inputs:
        event = info from the click
    '''
        
    combined = fits.open('mag_flux_tables/mag_flux_table_best.fits')
    tbdata = combined[1].data
    
    ### remove values that are +/-99 ###
    fluxn = k_mag_flux.mag_stacks(tbdata, aper=5)
    fluxn, tbdata = flux_funcs.no99(fluxn, tbdata)
    
    ob = tbdata['NUMBER_05B'][event.ind] #Define the object number from the index of the selected point
    if len(ob) > 1: #Reject selection if more than one object has been selected
        print('Too many objects selected')
        return
    print('Object identified')
    lightcurve_funcs.lightcurve5(ob, tbdata) #Plot the lightcurve from the lightcurve function
 
def onpickflux(event):
    ''' Function that plots the lightcurve of an object when it is clicked on 
    the vairiability v flux plot 
    This version creates 3 arcsec lightcurves using the self-calibrated 
    uncertainties and the final convolved K-band data
    Inputs:
        event = info from click
    '''
    print('Click registered')
    tbdata = fits.open('mag_flux_tables/mag_flux_table_best_extra_clean_no06.fits')[1].data
    sigtb = Table.read('sigma_tables/quad_epoch_sigma_table_extra_clean_no06.fits')
    tbdata = field_funcs.remove_edges(tbdata)

    ### remove values that are +/-99 ###
    fluxn = k_mag_flux.flux5_stacks(tbdata)
    fluxn, tbdata = flux_funcs.noneg(fluxn, tbdata)
    flux, fluxerr, tbdata = k_mag_flux.create_quad_error_array(sigtb, tbdata)
    ob = tbdata['NUMBER_05B'][event.ind] #Define the object number from the index of the selected point

    ### reset X-ray column as messed up by stacking ###
    tbdata['X-ray'][tbdata['X-ray']==70] = False 
    tbdata['X-ray'][tbdata['X-ray']==84] = True
    
    if len(ob) > 1: #Reject selection if more than one object has been selected
        print('Too many objects selected')
        return
    print('Object identified')
    
#    flux,fluxerr = normalise_flux_and_errors(flux, fluxerr)
    
    print('Finding object values')
    obflux = np.reshape(flux[tbdata['NUMBER_05B']==ob],np.shape(flux)[1])
    obfluxerr = np.reshape(fluxerr[tbdata['NUMBER_05B']==ob],np.shape(flux)[1])
    obdata = tbdata[tbdata['NUMBER_05B']==ob]
    print(obfluxerr)
    
    #set up time variable for plot
    print('setting up plot')
    t = np.linspace(1, 8, num=8)
    years = ['05B', '06B', '07B', '08B', '09B', '10B', '11B', '12B']
    x = [1,3,4,5,6,7,8]
#    chisq =my_chisquare_err(obflux, obfluxerr)
    
    plt.figure()
    if obdata['X-ray'] == True:
        print('Plotting x-ray')
        plt.errorbar(x, obflux, yerr=obfluxerr,fmt='o', color='r')
    else:
        print('Plotting non x-ray')
        plt.errorbar(x, obflux, yerr=obfluxerr,fmt='o', color='b')
#        print('Plotted non x-ray')
    plt.xlabel('Semester')
    plt.ylabel('Flux')
    plt.title('Lightcurve of Object '+str(obdata['NUMBER_05B']))#+' '+r' $\chi^{2} = $'+str(round(chisq, 2)))
    plt.xticks(t, years)
    plt.tight_layout()
    print('I got here')  
 
def onpickflux_2arcsec(event):
    ''' Function that plots the lightcurve of an object when it is clicked on 
    the vairiability v flux plot  
    This version creates 2 arcsec lightcurves using the self-calibrated 
    uncertainties and the final convolved K-band data
    Inputs:
        event = info from click
    '''
    print('Click registered')
    tbdata = fits.open('mag_flux_tables/K/mag_flux_table_best_extra_clean_no06.fits')[1].data
    sigtb = Table.read('sigma_tables/quad_epoch_sigma_table_extra_clean_no06_2arcsec.fits')
    tbdata = field_funcs.remove_edges(tbdata)

    ### remove values that are +/-99 ###
    fluxn = k_mag_flux.flux4_stacks(tbdata)
    fluxn, tbdata = flux_funcs.noneg(fluxn, tbdata)
#    fluxn, tbdata = remove_low_flux(fluxn, tbdata)
    flux, fluxerr, tbdata = k_mag_flux.create_quad_error_array(sigtb, tbdata, aper=4)
    ob = tbdata['NUMBER_05B'][event.ind] #Define the object number from the index of the selected point

    ### reset X-ray column as messed up by stacking ###
    tbdata['X-ray'][tbdata['X-ray']==70] = False 
    tbdata['X-ray'][tbdata['X-ray']==84] = True
    
    if len(ob) > 1: #Reject selection if more than one object has been selected
        print('Too many objects selected')
        return
    print('Object identified')
    
#    flux,fluxerr = normalise_flux_and_errors(flux, fluxerr)
    
    print('Finding object values')
    obflux = np.reshape(flux[tbdata['NUMBER_05B']==ob],np.shape(flux)[1])
    obfluxerr = np.reshape(fluxerr[tbdata['NUMBER_05B']==ob],np.shape(flux)[1])
    obdata = tbdata[tbdata['NUMBER_05B']==ob]
    print(obfluxerr)
    
    #set up time variable for plot
    print('setting up plot')
    t = np.linspace(1, 8, num=8)
    years = ['05B', '06B', '07B', '08B', '09B', '10B', '11B', '12B']
    x = [1,3,4,5,6,7,8]
#    chisq =my_chisquare_err(obflux, obfluxerr)
    
    plt.figure()
    if obdata['X-ray'] == True:
        print('Plotting x-ray')
        plt.errorbar(x, obflux, yerr=obfluxerr,fmt='o', color='r')
    else:
        print('Plotting non x-ray')
        plt.errorbar(x, obflux, yerr=obfluxerr,fmt='o', color='b')
#        print('Plotted non x-ray')
    plt.xlabel('Semester')
    plt.ylabel('Flux')
    plt.title('Lightcurve of Object '+str(obdata['NUMBER_05B']))#+' '+r' $\chi^{2} = $'+str(round(chisq, 2)))
    plt.xticks(t, years)
    plt.tight_layout()
    print('I got here')
    
    
def flux_variability_plot(flux, fluxchan, plottype, flux2 = [], fluxchan2 = [],
                          fluxerr = [], fluxerr2 = [], chanerr = [], chanerr2 = [], 
                          starflux=[], starfluxerr=[], comparison = False, 
                          normalised = False, stars=False, psfcorrect=False, 
                          scale=''):
    ''' Function to plot the variability vs mean flux plot. It can display MAD, 
    excess variance, chi squared, or variance as the variability statistic.
    If using MAD, specify plottype = 'mad' and supply just the flux and chandra 
    flux arrays
    If using Excess Variance, specify plottype = 'excess' and supply the error
    arrays as well as the flux ones
    If using chi squared, specify plottype = 'chisq' and supply the error
    arrays as well as the flux ones
    If using variance, specify plottype = 'var' and supply just the flux and 
    chandra flux arrays
    If you want to have a comparison plot of any type make sure to specify 
    the greyed out one first and then the coloured one 
    If you want to include stars make sure you specify this and supply the
    appropriate flux arrays.
    You can also specify what type of yscale you would like
    Inputs:
        flux = array of flux values for UDS objects, each row is one lightcurve
        fluxchan = array of flux values for chandra objects
        plottype = one of 'mad', 'excess', 'chisq', or 'var' which defines what 
                   statistic is used to quantify variability
        flux2 = optional second array of UDS objects if comparing two sets 
                (e.g. a before and after) this array will be the top coloured
                layer.
        fluxchan2 = optional second array of chandra objects 
        fluxerr = optional array of flux errors for objects, used if 
                    'excess' or 'chisq' is specified
        fluxerr2 = optional second array of flux errors for objects, used 
                    if 'excess' or 'chisq' and 'comparision = True' are specified
        chanerr = optional array of flux errors for X-ray objects, used if 
                    'excess' or 'chisq' is specified
        chanerr2 = optional second array of flux errors for X-ray objects, used 
                    if 'excess' or 'chisq' and 'comparision = True' are specified
        starflux = optional array of fluxes for stars
        starfluxerr = optional array of flux errors for objects, used if 
                        'excess' or 'chisq' and 'stars=True' are specified
        comparison = bool, specifies whether a comparison plot is required.
                     Default is False
        normalised = bool, specifies if the average flux of the light curves
                     should be normalised to 1. Default is False
        stars = bool, specifies if stars should be added to the plot. Default 
                is False
        psfcorrect = bool, specifies if a crude seeing correction should be
                     applied to the light curves. Default is False
        scale = string specifying how to scale the y axis:
                    '' = default, linear scale
                    'log' = log scale
                    'symlog' = symlog scale with limit threshold of 0.0001
    Output:
        fig = figure handle to allow clicking for light curves to be enabled if
                required 
        vary = array of the variability statistic for each object in the flux 
               array.
    '''
    
    fig = plt.figure(figsize=[8,8])
    avgfluxperob = np.nanmean(flux, axis=1) #for UDS
    avgfluxchanperob = np.nanmean(fluxchan, axis=1) #for non-stellar chandra
    if stars == True:
        savgfluxperob = np.nanmean(starflux, axis=1) #for stars

    ### Check if normalisation is true and normalise if necessary ###
    if normalised == True:
        if plottype == 'mad':
    #        flux = flux_funcs.normalise_mag(flux)
    #        fluxchan = flux_funcs.normalise_mag(fluxchan) 
    #        if stars == True:
    #            starflux = flux_funcs.normalise_mag(starflux)
            flux = flux_funcs.normalise_flux(flux)
            fluxchan = flux_funcs.normalise_flux(fluxchan) 
            if stars == True:
                starflux = flux_funcs.normalise_flux(starflux)
        else:
            flux, fluxerr = flux_funcs.normalise_flux_and_errors(flux, fluxerr)
            fluxchan, chanerr = flux_funcs.normalise_flux_and_errors(fluxchan, chanerr)
            if stars == True:
                starflux, starfluxerr = flux_funcs.normalise_flux_and_errors(starflux, starfluxerr)
    if psfcorrect == True:
        fluxchan = correction_funcs.psf_correct_mag(flux, fluxchan, 'median')
        if stars == True:
            starflux = correction_funcs.psf_correct_mag(flux, starflux, 'median')
        flux = correction_funcs.psf_correct_mag(flux, flux, 'median')
    ### Find out which plot type is specified and calculate appropriate statistic ###
    if plottype == 'mad':
        vary = median_absolute_deviation(flux, axis=1)
        varychan = median_absolute_deviation(fluxchan, axis=1)
        if stars == True:
            varystar = median_absolute_deviation(starflux, axis=1)
        plt.ylabel('MAD')
    elif plottype == 'excess':
        vary = vary_stats.normsigmasq(flux, fluxerr)
        varychan = vary_stats.normsigmasq(fluxchan, chanerr)
        if stars == True:
            varystar = vary_stats.normsigmasq(starflux, starfluxerr)
        plt.ylabel('Excess Variance')
    elif plottype == 'chisq':
        vary = vary_stats.my_chisquare_err(flux, fluxerr)
        varychan = vary_stats.my_chisquare_err(fluxchan, chanerr)
        if stars == True:
            varystar = vary_stats.my_chisquare_err(starflux, starfluxerr)
        plt.ylabel('Chi Squared')
    elif plottype == 'var':
        vary = np.var(flux, axis=1, ddof=1)
        varychan = np.var(fluxchan, axis=1, ddof=1)
        if stars == True:
            varystar = np.var(starflux, axis=1, ddof=1)
        plt.ylabel('Variance')
    else:
        print('Invalid plottype') #returns if unrecognised value is entered
        return
    
    ### Plot the variability v mean as appropriate ###
    if comparison == True:     
        
        avgfluxperob2 = np.nanmean(flux2, axis=1) #for UDS
        avgfluxchanperob2 = np.nanmean(fluxchan2, axis=1) #for non-stellar chandra

        if plottype == 'mad':
            if normalised == True:
    #            flux2 = flux_funcs.normalise_mag(flux2)
    #            fluxchan2 = flux_funcs.normalise_mag(fluxchan2)     
                flux2 = flux_funcs.normalise_flux(flux2)
                fluxchan2 = flux_funcs.normalise_flux(fluxchan2) 
            varycorr = median_absolute_deviation(flux2, axis=1)
            varychancorr = median_absolute_deviation(fluxchan2, axis=1)
        elif plottype == 'excess':
            if normalised == True:
                flux2, fluxerr2 = flux_funcs.normalise_flux_and_errors(flux2, fluxerr2)
                fluxchan2, chanerr2 = flux_funcs.normalise_flux_and_errors(fluxchan2, chanerr2)
            varycorr = vary_stats.normsigmasq(flux2, fluxerr2)
            varychancorr = vary_stats.normsigmasq(fluxchan2, chanerr2)

        ### plot varibility v flux graph for original in gray ###
        plt.plot(avgfluxperob, vary, '+', color =  'tab:gray', label='Galaxy', alpha = 0.5) 
        plt.plot(avgfluxchanperob, varychan, 'o', color =  'tab:gray',  mfc = 'none', markersize = 10,
                 label='X-ray detected', alpha = 0.5) #no picker as will be selected in the UDS point
        
        ### plot varibility v flux graph for new in colour ###
        line, = plt.plot(avgfluxperob2, varycorr, 'b+', label='Galaxy Convolved', picker=2, alpha = 0.5) #label tells 
        # the legend what to write, picker is the readius of pixels which a user must 
        # click within to call the object
        plt.plot(avgfluxchanperob2, varychancorr, 'ro', mfc = 'none', markersize = 10,
                 label='X-ray detected', alpha = 0.5) #no picker as will be selected in the UDS point
    else:
        if stars==True:
            plt.plot(savgfluxperob, varystar, 'm*', mfc = 'none', markersize = 10,
                     label='DR11 Star') 
        line, = plt.plot(avgfluxperob, vary, 'b+', label='Galaxy', picker=2)
        plt.plot(avgfluxchanperob, varychan, 'ro', mfc = 'none', markersize = 10,
                 label='X-ray detected') #no picker as will be selected in the UDS point

        
    ### Apply required plot charateristics ###
    plt.xscale('log')
    if scale == 'log':
        plt.yscale('log')
    elif scale == 'symlog':
        plt.yscale('symlog', linthreshy=0.0001)
        
    plt.ylim(3e-2,3e4)
    plt.xlim(8e1, 1e7)
#    plt.xlim(13,26)
    plt.xlabel('Mean Flux')
#    plt.xlabel('Mean Magnitude')
    plt.legend()
#    plt.gca().invert_xaxis()
    return fig, vary


def plot_median_line(fluxn, tbdata, statistic='MAD',createplot=True):
    ''' Function to find (and plot a line showing) the median value for a 
    variety of statistics across the flux range of the sample - useful when 
    calculating uncertainties.
    Inputs:
        fluxn = 2D array of flux values where each line is the light curve of
                an object
        tbdata = table of data that corresponds to the fluxes given (same length)
        statisitic = which statistic to find the median of. Options are MAD 
                     ('MAD', default), excess variance ('excess'), standard
                     deviation ('std'), variance ('var')
        createplot = bool, defines whether or not to actually plot the median
                     line onto the most recent plot. Default is True
    Outputs:
        bins = array denoting what the bin edges were
        allmedstat = array of the median values for each bin
    '''
    bins = np.array([13, 15])
    bins = np.append(bins, np.arange(16,24,0.2))
    bins = np.append(bins, [24])
    
    bins = 10**((30-bins)/2.5)
    bins = np.flip(bins, axis=0)
    #bins = bins[16:44] #because of flux limit
    
    ### Bin data ###
    allmedstat = np.array([])
    for n, binedge in enumerate(bins):
    #    print(binedge)
        if n==np.size(bins)-1:
            break
        mag, bindata = flux_funcs.fluxbin(binedge, bins[n+1], fluxn, tbdata) #bindata
        if statistic == 'excess':
            magerr = k_mag_flux.fluxerr5_stacks(bindata) #make error array
            nmag, nmagerr = flux_funcs.normalise_flux_and_errors(mag, magerr)
        else:
            nmag = flux_funcs.normalise_flux(mag)
        
        if statistic == 'std':
            binstat = np.std(nmag, axis=1)
        elif statistic == 'excess':
            binstat = vary_stats.normsigmasq(nmag, nmagerr)
        elif statistic == 'MAD':
            binstat = median_absolute_deviation(nmag, axis=1)
        elif statistic == 'var':
            binstat = np.var(nmag, axis=1, ddof=1)
        else:
            print('Unrecognised statistic entered')
            return
        statmed = np.nanmedian(binstat)
        allmedstat = np.append(allmedstat, statmed)
    
    if createplot==True:
        plt.plot(bins[0:42], allmedstat, 'k--')
    return bins, allmedstat

def plot_median_line_stars(fluxn, tbdata, sflux, sdata, statistic='MAD',
                           createplot=True):
    ''' Function to find (and plot a line showing) the median value for a 
    variety of statistics across the flux range of the sample - useful when 
    calculating uncertainties. This version includes the fluxes of stars in the
    calculation.
    Inputs:
        fluxn = 2D array of flux values where each line is the light curve of
                an object
        tbdata = table of data that corresponds to the fluxes given (same length)
        sflux = 2D array of flux values where each line is the light curve of
                a star
        sdata = table of data that corresponds to the star fluxes given (same 
                length)
        statisitic = which statistic to find the median of. Options are MAD 
                     ('MAD', default), excess variance ('excess'), standard
                     deviation ('std'), variance ('var')
        createplot = bool, defines whether or not to actually plot the median
                     line onto the most recent plot. Default is True
    Outputs:
        bins = array denoting what the bin edges were
        allmedstat = array of the median values for each bin
    '''
    bins = np.arange(13,26,0.2)
    bins = np.append(bins, [26,27,28])
#    bins = np.append([13,14], bins)
    
    bins = 10**((30-bins)/2.5)
    bins = np.flip(bins, axis=0)
    #bins = bins[16:44] #because of flux limit
    
    ### Bin data ###
    allmedstat = np.array([])
    for n, binedge in enumerate(bins):
    #    print(binedge)
        if n==np.size(bins)-1:
            break
        
        # bin both set of data
        gmag, bindata = flux_funcs.fluxbin(binedge, bins[n+1], fluxn, tbdata) #bindata
        smag, sbindata = flux_funcs.fluxbin(binedge, bins[n+1], sflux, sdata) #bindata
        
        # combine into one array
        mag = np.vstack((gmag, smag))
        
        if statistic == 'excess':
            gmagerr = k_mag_flux.fluxerr5_stacks(bindata) #make error array
            smagerr = k_mag_flux.fluxerr5_stacks(sdata)
            magerr = np.vstack((gmagerr, smagerr)) #combine
            
            nmag, nmagerr = flux_funcs.normalise_flux_and_errors(mag, magerr)
        else:
            nmag = flux_funcs.normalise_flux(mag)
        
        if statistic == 'std':
            binstat = np.std(nmag, axis=1)
        elif statistic == 'excess':
            binstat = vary_stats.normsigmasq(nmag, nmagerr)
        elif statistic == 'MAD':
            binstat = median_absolute_deviation(nmag, axis=1)
        elif statistic == 'var':
            binstat = np.var(nmag, axis=1, ddof=1)
        else:
            print('Unrecognised statistic entered')
            return
        statmed = np.nanmedian(binstat)
        allmedstat = np.append(allmedstat, statmed)
        
    if createplot==True:
        plt.plot(bins[0:np.size(bins)-1], allmedstat, 'k--')
    return bins, allmedstat

def plot_median_line_stars_J(fluxn, tbdata, sflux, sdata, statistic='MAD'):
    ''' Function to find (and plot a line showing) the median value for a 
    variety of statistics across the flux range of the sample - useful when 
    calculating uncertainties. This version includes the fluxes of stars in the
    calculation and is designed to be used on the J band data
    Inputs:
        fluxn = 2D array of flux values where each line is the light curve of
                an object
        tbdata = table of data that corresponds to the fluxes given (same length)
        sflux = 2D array of flux values where each line is the light curve of
                a star
        sdata = table of data that corresponds to the star fluxes given (same 
                length)
        statisitic = which statistic to find the median of. Options are MAD 
                     ('MAD', default), excess variance ('excess'), standard
                     deviation ('std'), variance ('var')
        createplot = bool, defines whether or not to actually plot the median
                     line onto the most recent plot. Default is True
    Outputs:
        bins = array denoting what the bin edges were
        allmedstat = array of the median values for each bin
    '''
    bins = np.arange(13,28,0.2)
    bins = np.append(bins, [28,29,30,31])
    
    bins = 10**((30-bins)/2.5)
    bins = np.flip(bins, axis=0)
    #bins = bins[16:44] #because of flux limit
    
    ### Bin data ###
    allmedstat = np.array([])
    for n, binedge in enumerate(bins):
    #    print(binedge)
        if n==np.size(bins)-1:
            break
        
        # bin both set of data
        gmag, bindata = flux_funcs.fluxbin(binedge, bins[n+1], fluxn, tbdata) #bindata
        smag, sbindata = flux_funcs.fluxbin(binedge, bins[n+1], sflux, sdata) #bindata
        
        # combine into one array
        mag = np.vstack((gmag, smag))
        
        if statistic == 'excess':
            gmagerr = j_mag_flux.fluxerr4_stacks(bindata) #make error array
            smagerr = j_mag_flux.fluxerr4_stacks(sdata)
            magerr = np.vstack((gmagerr, smagerr)) #combine
            
            nmag, nmagerr = flux_funcs.normalise_flux_and_errors(mag, magerr)
        else:
            nmag = flux_funcs.normalise_flux(mag)
        
        if statistic == 'std':
            binstat = np.std(nmag, axis=1)
        elif statistic == 'excess':
            binstat = vary_stats.normsigmasq(nmag, nmagerr)
        elif statistic == 'MAD':
            binstat = median_absolute_deviation(nmag, axis=1)
        elif statistic == 'var':
            binstat = np.var(nmag, axis=1, ddof=1)
        else:
            print('Unrecognised statistic entered')
            return
        statmed = np.nanmedian(binstat)
        allmedstat = np.append(allmedstat, statmed)
        
#    plt.plot(bins[0:np.size(bins)-1], allmedstat, 'k--')
    return bins, allmedstat