#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 10:59:07 2018

Code to create UVJ diagram with correct filter shifts

@author: ppxee
"""

### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
from astropy.table import Table #for handling tables
import numpy as np #for handling arrays
#import math
#from astropy.stats import median_absolute_deviation
from scipy import stats
import vari_funcs #my module to help run code neatly
#from scipy.stats import chisquare
plt.close('all') #close any open plots
#plt.style.use('default')

### Get data ###
nontb = Table.read('variable_tables/no06_variables_chi30_2arcsec_noXray_DR11data_restframe.fits')
xraytb = Table.read('variable_tables/no06_variables_chi30_2arcsec_Xray_DR11data_restframe.fits')
sc = Table.read('variable_tables/no06_variables_chi30_2arcsec_DR11data_restframe_SC.fits')
fullxray = Table.read('mag_flux_tables/novarys_chanDR11data_restframe_mag_flux_table_best_extra_clean_no06.fits')
fullUDS = Table.read('UDS_catalogues/DR11-2arcsec-June24-2018+plusXY_UVJ-SC.fits')

### Limit to Chandra region for simplicity ###
#datatb = vari_funcs.chandra_only(datatb)
#fullxray = vari_funcs.chandra_only(fullxray)
#chandata = vari_funcs.chandra_only(chandata)
#fullUDS = vari_funcs.chandra_only(fullUDS)

def get_data(tbdata, meanu, meanv, meanj):
    
    ### Get UDS UVJ values ###
    umag = tbdata['M_u_z_p'] + meanu
    vmag = tbdata['M_V_z_p'] + meanv
    jmag = tbdata['M_J_z_p'] + meanj

    u_v = umag - vmag
    v_j = vmag - jmag
    
    return u_v, v_j


### Check the difference between best z in DR11 cat and z in SC cat ###
z = sc['z_spec']#[mask]
z[z==-1] = sc['z_p'][z==-1]
sc_z = sc['ZUSED']

### Ignore if |diff| > 0.1 ###
diffz = z - sc_z
testdata = sc[np.abs(diffz) < 0.1]
testdata = testdata[testdata['U'] < 0]

### Get UDS UVJ values ###
umag = testdata['M_u_z_p']
vmag = testdata['M_V_z_p']
jmag = testdata['M_J_z_p']

### Get Super Colour UVJ ##'
scumag = testdata['U']
scvmag = testdata['V']
scjmag = testdata['J']

### Find the difference between the two values ###
diffu = scumag - umag
diffv = scvmag - vmag
diffj = scjmag - jmag

meanu = np.mean(diffu) 
meanv = np.mean(diffv) 
meanj = np.mean(diffj) 

### remove those with very high phot z ###
#xraytb = xraytb[xraytb['z_p']<6]
#nontb = nontb[nontb['z_p']<6]

xraytbup = vari_funcs.flux_split(xraytb, 'upper')
nontbup = vari_funcs.flux_split(nontb, 'upper')
xraytblow = vari_funcs.flux_split(xraytb, 'lower')
nontblow = vari_funcs.flux_split(nontb, 'lower')

allxu_v, allxv_j = get_data(fullxray, meanu, meanv, meanj)
noxu_v, noxv_j = get_data(nontb, meanu, meanv, meanj)
xu_v, xv_j = get_data(xraytb, meanu, meanv, meanj)

### Use actual SC for UDS sources ###
fullumag = fullUDS['U'][~fullUDS['Stars-DR11']]
fullvmag = fullUDS['V'][~fullUDS['Stars-DR11']]
fulljmag = fullUDS['J'][~fullUDS['Stars-DR11']]
fullu_v = fullumag - fullvmag
fullv_j = fullvmag - fulljmag

noxu_vup, noxv_jup = get_data(nontbup, meanu, meanv, meanj)
xu_vup, xv_jup = get_data(xraytbup, meanu, meanv, meanj)

noxu_vlow, noxv_jlow = get_data(nontblow, meanu, meanv, meanj)
xu_vlow, xv_jlow = get_data(xraytblow, meanu, meanv, meanj)

#%%
plt.figure(figsize=[8,8])
plt.plot(fullv_j, fullu_v, '.',markersize=0.5, color='tab:grey', alpha=0.5, label='UDS Galaxies')
plt.plot(allxv_j, allxu_v, '+', color='k', label='Non Variable Chandra Sources')
plt.plot(noxv_j, noxu_v,'bo', label='Non X-ray Variable Sources')
plt.plot(xv_j, xu_v,'ro', label='X-ray Variable Sources')

#plt.plot(noxv_jup, noxu_vup,'go', label='High flux, not X-ray')
#plt.plot(xv_jup, xu_vup,'gd', label='High flux, X-ray')
#
#plt.plot(noxv_jlow, noxu_vlow,'mo', label='Low flux, not X-ray')
#plt.plot(xv_jlow, xu_vlow,'md', label='Low flux, X-ray')

# plot dividing lines => y = 0.875x + 0.6 #
x = np.linspace(0.8,1.6)
y = 0.875*x + 0.6
plt.plot(x,y,'k--')
plt.hlines(1.3,-1,0.8,linestyles='dashed')
plt.vlines(1.6,2.0,2.4,linestyles='dashed')

plt.xlim(xmin=-1, xmax=3)
plt.ylim(ymin=-0.5,ymax=2.4)
plt.xlabel('V - J')
plt.ylabel('U - V')
plt.legend(loc='upper left', frameon=False)
plt.tight_layout()
#plt.savefig('plots/Chi_variables/no06_extra_clean/uvjdiagram_Rmatched.png')
