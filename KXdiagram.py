#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 19 12:24:37 2018

Code to create a KX diagram

@author: lizzieelmer
"""
import matplotlib.colors
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
#from astropy.table import Table #for handling tables
import numpy as np #for handling arrays
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
#from astropy.stats import median_absolute_deviation
#import vari_funcs #my module to help run code neatly
plt.close('all') #close any open plots

### Define cosmology ###
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

### get the data ###
tbdata = fits.open('UDS_catalogues/DR11-pointlike.fits')[1].data
varydata = fits.open('variable_tables/no06_variables_chi30_DR11data_restframe.fits')[1].data
#varydata = fits.open('variable_tables/no06_variables_chi30_chandata_DR11data_restframe.fits')[1].data
varypoint = fits.open('variable_tables/no06_variables_chi30_smallpoint.fits')[1].data
alldata = fits.open('UDS_catalogues/DR11-extended-gals-22.fits')[1].data

def get_v_j_k(tbdata):
    kmag = tbdata['M_K_z_p'] - 1.9
    jmag = tbdata['M_J_z_p'] - 0.938
    vmag = tbdata['M_V_z_p']
    v_j = vmag - jmag
    j_k = jmag - kmag
    return v_j, j_k
    
# find the bright galaxies
#alldatabright = alldata[alldata['M_K_z_p']<19]

### Get V-J and J-K ###
v_j, j_k = get_v_j_k(tbdata)
varyv_j, varyj_k = get_v_j_k(varydata)
pointv_j, pointj_k = get_v_j_k(varypoint)
allv_j, allj_k = get_v_j_k(alldata)
varyxv_j = varyv_j[varydata['X-ray']]
varyxj_k = varyj_k[varydata['X-ray']]
varynv_j = varyv_j[~varydata['X-ray']]
varynj_k = varyj_k[~varydata['X-ray']]

### get line equation ###
x = np.linspace(-1,3,100)
y = (x-0.18)/0.36

### Get X-ray lum ###
#z = varydata['z_spec']#[mask]
#z[z==-1] = varydata['z_p'][z==-1]
#DL = cosmo.luminosity_distance(z)
#DL = DL.to(u.cm)
##F = varydata['RFLUX_20']
##L = F*4*np.pi*(DL.value**2)
##nL = L[~varydata['X-ray']]
##xL = L[varydata['X-ray']]
#xrayF = varydata['Full_flux']#[chanmask]
#xrayL = xrayF*4*np.pi*(DL.value**2)
#
#### get optical lum???? L=L0 * 10^(0.4M)###
#M = varydata['RMAG_20']
#L0 = 3.0128e28
#L = L0 * (10**(0.4*M))

### Plot ###
plt.figure(figsize=[8,8])
plt.plot(allj_k,allv_j,'+', color='tab:grey',label='UDS Galaxy',alpha=0.5)
plt.plot(j_k,v_j,'o', color='k',label='UDS Point Source',alpha=0.5, markersize=5)
plt.plot(varynj_k,varynv_j,'bo',label='Non X-ray Variable Source')
plt.plot(varyxj_k,varyxv_j,'ro',label='X-ray Variable Source')
#plt.scatter(varyxj_k,varyxv_j,c=xrayL/L,marker='o',label='X-ray Variable Source',
#            norm=matplotlib.colors.LogNorm(),zorder=5)
plt.plot(pointj_k,pointv_j,'ks',mfc='None',markersize=8, label='Centrally Point-Like Variable Source')
plt.plot(x,y,'k--')
plt.xlim(xmin=0,xmax=3)
plt.ylim(ymin=-1,ymax=6)
plt.gca().invert_yaxis()
plt.xlabel('J - K')
plt.ylabel('V - J')
plt.legend()
#plt.colorbar()
plt.tight_layout()