#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 15:48:46 2019

Code to look at the stellarity of the variable sample

@author: ppxee
"""

### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
from astropy.table import Table #for handling tables
import numpy as np #for handling arrays
import vari_funcs #my module to help run code neatly
plt.close('all')

### Import fits files ###
noxvarydata = Table.read('variable_tables/K/old_catalogue/no06_variables_chi30_2arcsec_noXray_DR11data_restframe_DR11_photometry.fits')
xvarydata = Table.read('variable_tables/K/old_catalogue/no06_variables_chi30_2arcsec_Xray_DR11data_restframe_DR11_photometry.fits')
fullxray = Table.read('UDS_catalogues/chandra_catalogue_DR11_photometry_novarys.fits')
dr11 = Table.read('UDS_catalogues/DR11_output_best.fits')
dr11stars = Table.read('UDS_catalogues/DR11_output_stars.fits')

### get data ##
def get_data(tbdata):
    small = tbdata['MAG_APER'][:,0] # 0.7 arcsec
    big = tbdata['MAG_APER'][:,3] # 2 arcsec
    
    ### set values where mag == 99 to nan ### 
    small[small==99] = np.nan
    big[big==99] = np.nan
    
    ### compute difference ###
    stellar = big - small
    
    return big, small, stellar

noxbig, noxsmall, noxstellar = get_data(noxvarydata)
xbig, xsmall, xstellar = get_data(xvarydata)
allxbig, allxsmall, allxstellar = get_data(fullxray)
big, small, stellar = get_data(dr11)
sbig, ssmall, sstellar = get_data(dr11stars)

#%%
bins=np.linspace(-2.1,0,50)
plt.figure(figsize=[9,6])
plt.hist(xstellar, bins , color='r', density=True, histtype='step', 
         label='Variable X-ray Source', linestyle='-', linewidth=1.5)
plt.hist(noxstellar, bins , color='b', density=True, histtype='step', 
         label='Variable Non-X-ray Source', linestyle='-.', linewidth=1.5)
plt.hist(allxstellar, bins, color='k', density=True, histtype='step',
         label='Non-Variable X-ray Source', linestyle='--', linewidth=1.5)
plt.hist(stellar, bins, color='tab:grey', density=True, histtype='stepfilled',
         label='Galaxy', alpha=0.5, zorder=0)
plt.hist(sstellar, bins, color='m', density=True, histtype='stepfilled',
         label='Star', alpha=0.5, zorder=0)
plt.xlabel(r'$K_{2^{\prime\prime}} - K_{0.7^{\prime\prime}}$')
plt.legend()
plt.tight_layout()

#%% plot stellarity vs mag figure
plt.figure(figsize=[8,6])
plt.plot(big, stellar, '.', color='tab:grey',
         label='Galaxy')
plt.plot(sbig, sstellar, '*', color='m',
         label='Star')
plt.plot(allxbig, allxstellar, '+', color='k',
         label='Non-Variable X-ray Source')
plt.plot(noxbig, noxstellar, 'o', color='b', 
         label='Variable Non-X-ray Source')
plt.plot(xbig, xstellar, 's', color='r', 
         label='Variable X-ray Source')
plt.ylabel(r'$K_{2^{\prime\prime}} - K_{0.7^{\prime\prime}}$')
plt.xlabel(r'$K_{2^{\prime\prime}}$')
plt.xlim(12, 24.5)
plt.ylim(-2.2, 0.15)
plt.legend()
plt.tight_layout()

#%% plot r-i vs g-r
def get_colour_data(tbdata):
    r = tbdata['RMAG_20']
    i = tbdata['iMAG_20']
    g = tbdata['BMAG_20'] # don't have g, used b last time but check
    
    r_i = r - i 
    g_r = g - r
    return r, i, g, r_i, g_r

r, i, g, r_i, g_r = get_colour_data(dr11)
sr, si, sg, sr_i, sg_r = get_colour_data(dr11stars)
allxr, allxi, allxg, allxr_i, allxg_r = get_colour_data(fullxray)
xr, xi, xg, xr_i, xg_r = get_colour_data(xvarydata)
noxr, noxi, noxg, noxr_i, noxg_r = get_colour_data(noxvarydata)

plt.figure(figsize=[8,6])
plt.plot(g_r, r_i, '.', color='tab:grey', markersize=1, alpha=0.2,
         label='Galaxy')
plt.plot(sg_r, sr_i, '*', color='m',
         label='Star')
plt.plot(allxg_r, allxr_i, '+', color='k',
         label='Non-Variable X-ray Source')
plt.plot(noxg_r, noxr_i, 'o', color='b', 
         label='Variable Non-X-ray Source')
plt.plot(xg_r, xr_i, 's', color='r', 
         label='Variable X-ray Source')

plt.ylabel('r - i')
plt.xlabel('b - r')
plt.xlim(-2, 4)
plt.ylim(-1.5, 4)
plt.legend(loc='upper left')
plt.tight_layout()


#%% plot z-K vs B-z
def get_colour_data_BzK(tbdata):
    z = tbdata['zMAG_20']
    K = tbdata['KMAG_20']
    B = tbdata['BMAG_20'] 
    
    ### restrict to z > 20 and B>20 ###
    mask1 = z > 20
    mask2 = B > 20
    mask3 = z != 99
    mask4 = B != 99
    mask5 = K != 99
    mask = mask1*mask2*mask3*mask4*mask5.astype(bool)
    
    z = z[mask]
    K = K[mask]
    B = B[mask]
    
    z_K = z - K 
    B_z = B - z
    return z, K, B, z_K, B_z, mask

z, K, B, z_K, B_z, _ = get_colour_data_BzK(dr11)
sz, sK, sB, sz_K, sB_z, _ = get_colour_data_BzK(dr11stars)
allxz, allxK, allxB, allxz_K, allxB_z, _ = get_colour_data_BzK(fullxray)
xz, xK, xB, xz_K, xB_z, _ = get_colour_data_BzK(xvarydata)
noxz, noxK, noxB, noxz_K, noxB_z, mask = get_colour_data_BzK(noxvarydata)

plt.figure(figsize=[8,6])
plt.plot(B_z, z_K, '.', color='tab:grey', markersize=1, alpha=0.2,
         label='Galaxy')
plt.plot(sB_z, sz_K, '*', color='m',
         label='Star')
plt.plot(allxB_z, allxz_K, '+', color='k',
         label='Non-Variable X-ray Source')
plt.plot(noxB_z, noxz_K, 'o', color='b', 
         label='Variable Non-X-ray Source')
plt.plot(xB_z, xz_K, 's', color='r', 
         label='Variable X-ray Source')

### plot defining line ###
x = np.linspace(-0.5,6,10)
y = 0.3*x - 0.5
plt.plot(x,y,'k--')

plt.ylabel('z - K')
plt.xlabel('B - z')
plt.xlim(-0.5, 6)
plt.ylim(-1.5, 6)
plt.legend(loc='upper right')
plt.tight_layout()

### extract source that is in stellar locus ###
#mask_source = noxz_K < ((0.3*noxB_z) - 0.5)
#data = noxvarydata[mask][mask_source]
#testr, testi, testg, testr_i, testg_r = get_colour_data(data)
#plt.figure(3)
#plt.plot(testg_r, testr_i, 'd', color='g', markersize = 10,
#         label='Variable Non-X-ray Source')
#plt.figure(2)
#testbig, testsmall, teststellar = get_data(data)
#plt.plot(testbig, teststellar, 'd', color='g',
#         label='Star')