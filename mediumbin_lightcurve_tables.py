#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 17:23:43 2018

@author: ppxee
"""


### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
import matplotlib.animation as animation
#plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'
from astropy.io import fits #for handling fits
#from astropy.table import Table #for handling tables
from photutils import CircularAperture, aperture_photometry
import numpy as np #for handling arrays
#import math
#from astropy.stats import median_absolute_deviation
#import vari_funcs_no06 #my module to help run code neatly
plt.close('all') #close any open plots
hdr08B = fits.getheader('Images/extra_clean_no06_UDS_08B_K.fits') # random year (same in all)
const = -hdr08B['CD1_1'] # constant that defines unit conversion for FWHM

n=1
binarr = range(11)
#jmag = np.empty(len(binarr))
#flarejmag = np.empty(len(binarr))
#jflux = np.empty(len(binarr))
#flarejflux = np.empty(len(binarr))
sejmag = np.empty(len(binarr))
sejmagerr = np.empty(len(binarr))
for n in binarr:
    if n == 0:
        tbdata = fits.open('SE_outputs_yearstacks/07B_output_J.fits')[1].data
        obdata = tbdata[tbdata['NUMBER']==57824]
        sejmag[n] = obdata['MAG_APER'][0][4]
        sejmagerr[n] = obdata['MAGERR_APER'][0][4]
    elif n == binarr[-2]:
#        tbdata = fits.open('SE_outputs_yearstacks/09B_output_J.fits')[1].data
#        obdata = tbdata[tbdata['NUMBER']==62243]
#        sejmag[n] = obdata['MAG_APER'][0][4]
#        sejmagerr[n] = obdata['MAGERR_APER'][0][4]
        sejmag[n] = np.nan
        sejmagerr[n] = np.nan
    elif n == binarr[-1]:
#        tbdata = fits.open('SE_outputs_yearstacks/10B_output_J.fits')[1].data
##        obdata = tbdata[tbdata['NUMBER']==62243]
#        sejmag[n] = tbdata['MAG_APER'][62243][4]
#        sejmagerr[n] = tbdata['MAGERR_APER'][62243][4]
        sejmag[n] = np.nan
        sejmagerr[n] = np.nan
    else:
        tbdata = fits.open('SE_outputs_shortstacks/08B_J_output_cam1_8bin_'+str(n-1)+'.fits')[1].data
        obdata = tbdata[tbdata['NUMBER']==13363]
        sejmag[n] = obdata['MAG_APER'][0][4]
        sejmagerr[n] = obdata['MAGERR_APER'][0][4]


binarr = range(13)
#kmag = np.empty(len(binarr))
#flarekmag = np.empty(len(binarr))
#kflux = np.empty(len(binarr))
#flarekflux = np.empty(len(binarr))
sekmag = np.empty(len(binarr))
sekmagerr = np.empty(len(binarr))
for n in binarr:
    if n == 0:
        tbdata = fits.open('SE_outputs_yearstacks/07B_output.fits')[1].data
        obdata = tbdata[tbdata['NUMBER']==62242]
        sekmag[n] = obdata['MAG_APER'][0][4]
        sekmagerr[n] = obdata['MAGERR_APER'][0][4]
    elif n == binarr[-2]:
        tbdata = fits.open('SE_outputs_yearstacks/09B_output.fits')[1].data
        obdata = tbdata[tbdata['NUMBER']==62242]
        sekmag[n] = obdata['MAG_APER'][0][4]
        sekmagerr[n] = obdata['MAGERR_APER'][0][4]
    elif n == binarr[-1]:
        tbdata = fits.open('SE_outputs_yearstacks/10B_output.fits')[1].data
        obdata = tbdata[tbdata['NUMBER']==62242]
        sekmag[n] = obdata['MAG_APER'][0][4]
        sekmagerr[n] = obdata['MAGERR_APER'][0][4]
    else:
        tbdata = fits.open('SE_outputs_shortstacks/08B_output_cam1_10bin_'+str(n-1)+'.fits')[1].data
        obdata = tbdata[tbdata['NUMBER']==15100]
        sekmag[n] = obdata['MAG_APER'][0][4]
        sekmagerr[n] = obdata['MAGERR_APER'][0][4]
    

### Get H data ###
binarr = range(9)
sehmag = np.empty(len(binarr))
sehmagerr = np.empty(len(binarr))
for n in binarr:
    if n == 0:
#        tbdata = fits.open('SE_outputs_yearstacks/07B_output.fits')[1].data
#        obdata = tbdata[tbdata['NUMBER']==62242]
#        sekmag[n] = obdata['MAG_APER'][0][4]
#        sekmagerr[n] = obdata['MAGERR_APER'][0][4]
        sehmag[n] = np.nan
        sehmagerr[n] = np.nan
    elif n == binarr[-2]:
#        tbdata = fits.open('SE_outputs_yearstacks/09B_output.fits')[1].data
#        obdata = tbdata[tbdata['NUMBER']==62242]
#        sekmag[n] = obdata['MAG_APER'][0][4]
#        sekmagerr[n] = obdata['MAGERR_APER'][0][4]
        sehmag[n] = np.nan
        sehmagerr[n] = np.nan
    elif n == binarr[-1]:
#        tbdata = fits.open('SE_outputs_yearstacks/10B_output.fits')[1].data
#        obdata = tbdata[tbdata['NUMBER']==62242]
#        sekmag[n] = tbdata['MAG_APER'][62243][4]
#        sekmagerr[n] = tbdata['MAGERR_APER'][62243][4]
        sehmag[n] = np.nan
        sehmagerr[n] = np.nan
    else:
        tbdata = fits.open('SE_outputs_shortstacks/08B_H_output_cam1_6bin_'+str(n-1)+'.fits')[1].data
        obdata = tbdata[tbdata['NUMBER']==11900]
        sehmag[n] = obdata['MAG_APER'][0][4]
        sehmagerr[n] = obdata['MAGERR_APER'][0][4]
        
### figure out date stamps ###
### create list of date stamps between 20/9 and 25/11
initdates = np.empty(67).astype(str)
for x in range(67):
    if x <= 10:
        initdates[x] = str(20+x)+'/9'
    elif x > 10 and x <= 41:
        initdates[x] = str(x-10)+'/10'
    elif x > 41:
        initdates[x] = str(x-41)+'/11'
            
dates = initdates
#sem07 = np.repeat('07B',30)
#dates = np.append(sem07, initdates)
#sem09 = np.repeat('09B',30)
#dates = np.append(dates, sem09)
#sem10 = np.repeat('10B',30)
#dates = np.append(dates, sem10)
#Jdates = ['3/10','16/10','18/10','18/10','19/10','19/10',]
#Hdates = ['20/9','21/9','22/9','27/9','1/10','10/10','21/10','27/10','1/11',
#          '8/11','16/11','21/11']
Jxdata = np.array([18, 28, 29, 30, 33, 41, 42, 49])#+30
#Jxdata = np.append(0, Jxdata)
#Jxdata = np.append(Jxdata, 220)
Hxdata = np.array([0,3,16,35,45,61])#+30
#Hxdata = np.append(0, Hxdata)
#Hxdata = np.append(Hxdata, 220)
Kxdata = np.array([0, 1, 9, 11, 12, 13, 16, 21, 39, 56])#+30
#Kxdata = np.append(0, Kxdata)
#Kxdata = np.append(Kxdata, 220)
initticks = np.array([0, 10, 20, 30, 40, 50, 60])
ticks= initticks#+30
#ticks = np.append(0, ticks)
#ticks = np.append(ticks, [126])#,156])

#plt.figure(figsize=[17,7])
#plt.plot([0,260,290], [flarejflux[0],flarejflux[-2],flarejflux[-1]], 'bs')#, label='J semester flux')
#plt.plot([0,260,290], [flarehflux[0],flarehflux[-2],flarehflux[-1]], 'ks')#, label='H semester flux')
#plt.plot([0,260,290], [flarekflux[0],flarekflux[-2],flarekflux[-1]], 'rs')#, label='K semester flux')
#
#plt.plot(Jxdata, flarejflux[1:-2], 'bo', label='J flare flux')
#plt.plot(Hxdata, flarehflux[1:-2], 'ko', label='H flare flux')
#plt.plot(Kxdata, flarekflux[1:-2], 'ro', label='K flare flux')
#plt.legend()
#plt.xticks(ticks, dates[ticks], rotation='vertical')
#plt.ylabel('Flux of flare in 1.4" aperture')
#plt.xlabel('Date')
#plt.tight_layout()
#
#plt.figure(figsize=[17,7])
#plt.plot([0,260,290], [jflux[0],jflux[-2],jflux[-1]], 'bs')#, label='J semester flux')
#plt.plot([0,260,290], [hflux[0],hflux[-2],hflux[-1]], 'ks')#, label='H semester flux')
#plt.plot([0,260,290], [kflux[0],kflux[-2],kflux[-1]], 'rs')#, label='K semester flux')
#plt.plot(Jxdata, jflux[1:-2], 'bo', label='J object flux')
#plt.plot(Hxdata, hflux[1:-2], 'ko', label='H object flux')
#plt.plot(Kxdata, kflux[1:-2], 'ro', label='K object flux')
#plt.legend()
#plt.xticks(ticks, dates[ticks], rotation='vertical')
#plt.ylabel('Flux of object in 3" aperture')
#plt.xlabel('Date')
#plt.tight_layout()
##
##
##plt.figure(figsize=[17,7])
##plt.plot(Jxdata, flarejflux, 'o', label='J flare flux')
##plt.legend()
##plt.xticks(ticks, dates[ticks], rotation='vertical')
##
##plt.figure(figsize=[17,7])
##plt.plot(Hxdata, flarehflux, 'o', label='H flare flux')
##plt.legend()
##plt.xticks(ticks, dates[ticks], rotation='vertical')
##
##plt.figure(figsize=[17,7])
##plt.plot(Kxdata, flarekflux, 'o', label='K flare flux')
##plt.legend()
##plt.xticks(ticks, dates[ticks], rotation='vertical')
##
#plt.figure(figsize=[7,7])
#plt.plot([0,126,156], [flarejmag[0],flarejmag[-2],flarejmag[-1]], 'bs')#, label='J semester flux')
#plt.plot([0,126,156], [flarehmag[0],flarehmag[-2],flarehmag[-1]], 'ks')#, label='H semester flux')
#plt.plot([0,126,156], [flarekmag[0],flarekmag[-2],flarekmag[-1]], 'rs')#, label='K semester flux')
#plt.plot(Jxdata, flarejmag[1:-2], 'bo', label='J flare mag')
#plt.plot(Hxdata, flarehmag[1:-2], 'ko', label='H flare mag')
#plt.plot(Kxdata, flarekmag[1:-2], 'ro', label='K flare mag')
##plt.legend()
#plt.xticks(ticks, dates[ticks], rotation='vertical')
#plt.ylabel('Magnitude of flare in 1.4" aperture')
#plt.xlabel('Date')
#plt.gca().invert_yaxis()
#plt.tight_layout()
#
#plt.figure(figsize=[7,7])
##plt.errorbar([0,126,156], [sejmag[0],sejmag[-2],sejmag[-1]], 
##             [sejmagerr[0],sejmagerr[-2],sejmagerr[-1]], 'fmt=bs')#, label='J semester flux')
#plt.errorbar(0, sejmag[0], sejmagerr[0], fmt='bs')#, label='J semester flux')
##plt.plot([0,126,156], [hmag[0],hmag[-2],hmag[-1]], 'ks')#, label='H semester flux')
##plt.plot([0,126,156], [kmag[0],kmag[-2],kmag[-1]], 'rs')#, label='K semester flux')
#plt.errorbar([0,126,156], [sekmag[0],sekmag[-2],sekmag[-1]], 
#             [sekmagerr[0],sekmagerr[-2],sekmagerr[-1]], fmt='rs')#, label='J semester flux')
#plt.errorbar(Jxdata, sejmag[1:-2], sejmagerr[1:-2], fmt='bo', label='J object mag')
#plt.errorbar(Hxdata, sehmag[1:-2], sehmagerr[1:-2], fmt='ko', label='H object mag')
#plt.errorbar(Kxdata, sekmag[1:-2], sekmagerr[1:-2], fmt='ro', label='K object mag')
#plt.legend()
#plt.xticks(ticks, dates[ticks], rotation='vertical')
#plt.ylabel('Magnitude of object in 3" aperture')
#plt.xlabel('Date')
#plt.gca().invert_yaxis()
#plt.tight_layout()
#
#plt.figure(figsize=[7,7])
#plt.plot([0,126,156], [flarejmag[0],flarejmag[-2],flarejmag[-1]], 'bs')#, label='J semester flux')
#plt.plot(Jxdata, flarejmag[1:-2], 'bo', label='J flare mag')
#plt.legend()
#plt.xticks(ticks, dates[ticks], rotation='vertical')
#axes = plt.gca()
#ylims = axes.get_ylim()
#ymid = (ylims[1]+ylims[0])/2
#plt.ylim(ymin=ymid-1.5, ymax=ymid+1.5)
#plt.ylabel('J magnitude of flare in 1.4" aperture')
#plt.xlabel('Date')
#plt.gca().invert_yaxis()
#plt.tight_layout()
#
#plt.figure(figsize=[7,7])
#plt.plot([0,126,156], [flarehmag[0],flarehmag[-2],flarehmag[-1]], 'ks')#, label='H semester flux')
#plt.plot(Hxdata, flarehmag[1:-2], 'ko', label='H flare mag')
#plt.legend()
#plt.xticks(ticks, dates[ticks], rotation='vertical')
#axes = plt.gca()
#ylims = axes.get_ylim()
#ymid = (ylims[1]+ylims[0])/2
#plt.ylim(ymin=ymid-1.5, ymax=ymid+1.5)
#plt.ylabel('H magnitude of flare in 1.4" aperture')
#plt.xlabel('Date')
#plt.gca().invert_yaxis()
#plt.tight_layout()
#
#plt.figure(figsize=[7,7])
#plt.plot([0,126,156], [flarekmag[0],flarekmag[-2],flarekmag[-1]], 'rs')#, label='K semester flux')
#plt.plot(Kxdata, flarekmag[1:-2], 'ro', label='K flare mag')
#plt.legend()
#plt.xticks(ticks, dates[ticks], rotation='vertical')
#axes = plt.gca()
#ylims = axes.get_ylim()
#ymid = (ylims[1]+ylims[0])/2
#plt.ylim(ymin=ymid-1.5, ymax=ymid+1.5)
#plt.ylabel('K magnitude of flare in 1.4" aperture')
#plt.xlabel('Date')
#plt.gca().invert_yaxis()
#plt.tight_layout()

#%% calculate rest frame mags
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u

### Define cosmology ###
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

z = 1.51 #approximately
DL = cosmo.luminosity_distance(z)
DL = DL.to(u.pc)

M_g = sejmag - 5*(np.log10(DL.value)-1) - 2.5*np.log10(1+z)
M_r = sehmag - 5*(np.log10(DL.value)-1) - 2.5*np.log10(1+z)
M_z = sekmag - 5*(np.log10(DL.value)-1) - 2.5*np.log10(1+z)

M_g_err = M_g * (sejmagerr/sejmag)
M_r_err = M_r * (sehmagerr/sehmag)
M_z_err = M_z * (sekmagerr/sekmag)

### plot new curves ####
plt.figure(figsize=[7,7])
#plt.plot([0,126,156], [M_g[0],M_g[-2],M_g[-1]], 'bs')#, label='J semester flux')
#plt.plot([0,126,156], [M_r[0],M_r[-2],M_r[-1]], 'ks')#, label='H semester flux')
#plt.plot([0,126,156], [M_z[0],M_z[-2],M_z[-1]], 'rs')#, label='K semester flux')
#plt.plot([0,126], [M_g[0],M_g[-1]], 'bs', mfc='None')#, label='J semester flux')
#plt.plot([0,126], [M_r[0],M_r[-1]], 'ks', mfc='None')#, label='H semester flux')
#plt.plot([0,126], [M_z[0],M_z[-1]], 'rs', mfc='None')#, label='K semester flux')
plt.errorbar(Jxdata, M_g[1:-2], M_g_err[1:-2], fmt='bo', label='M_g object mag')
plt.errorbar(Hxdata, M_r[1:-2], M_r_err[1:-2], fmt='ko', label='M_r object mag')
plt.errorbar(Kxdata, M_z[1:-2], M_z_err[1:-2], fmt='ro', label='M_z object mag')
plt.legend()
plt.xticks(ticks, dates[ticks], rotation='vertical')
plt.ylabel('Restframe magnitude of flare in 3" aperture')
plt.xlabel('Date')
plt.gca().invert_yaxis()
plt.tight_layout()