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
plt.style.use('default')

### Get data ###
datatb = Table.read('variable_tables/no06_variables_chi30_DR11data_restframe.fits')
sc = Table.read('variable_tables/no06_variables_chi30_DR11data_restframe_SC.fits')
fullxray = Table.read('mag_flux_tables/novarys_chanDR11data_restframe_mag_flux_table_best_extra_clean_no06.fits')
xmmdata = fits.open('variable_tables/no06_variables_chi30_xmmdata_DR11data_restframe.fits')[1].data
chandata = Table.read('variable_tables/no06_variables_chi30_chandata_DR11data_restframe.fits')
fullUDS = Table.read('UDS_catalogues/DR11-2arcsec-June24-2018+plusXY_UVJ-SC.fits')

### Limit to Chandra region for simplicity ###
#datatb = vari_funcs.chandra_only(datatb)
#fullxray = vari_funcs.chandra_only(fullxray)
#chandata = vari_funcs.chandra_only(chandata)
#fullUDS = vari_funcs.chandra_only(fullUDS)

datatb = datatb[datatb['RMAG_20'] < 24]
datatb = datatb[datatb['RMAG_20'] > 20.3]


### Split the variable table into x-ray and non x-ray ###
xraytb = datatb[datatb['X-ray']]
nontb = datatb[~datatb['X-ray']]

### Plot R band distributions so can match for proposal ###
xray_R = xraytb['RMAG_20']
non_R = nontb['RMAG_20']
plt.hist([xray_R, non_R], color=['r','b'], histtype='step', 
         label = ['X-ray','Non X-ray'])
plt.xlabel('R Band Magnitude')
plt.ylabel('Counts')
plt.legend()
plt.tight_layout()

### Check the difference between best z in DR11 cat and z in SC cat ###
z = sc['z_spec']#[mask]
z[z==-1] = sc['z_p'][z==-1]
sc_z = sc['ZUSED']

### Ignore if |diff| > 0.1
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

#plt.figure()
#plt.hist([diffu,diffv,diffj],50)


### Get UDS UVJ values ###
umag = nontb['M_u_z_p'] + meanu
vmag = nontb['M_V_z_p'] + meanv
jmag = nontb['M_J_z_p'] + meanj
chanumag = fullxray['M_u_z_p'] + meanu
chanvmag = fullxray['M_V_z_p'] + meanv
chanjmag = fullxray['M_J_z_p'] + meanj
chanvaryumag = xraytb['M_u_z_p'] + meanu
chanvaryvmag = xraytb['M_V_z_p'] + meanv
chanvaryjmag = xraytb['M_J_z_p'] + meanj
fullumag = fullUDS['U'][~fullUDS['Stars-DR11']]
fullvmag = fullUDS['V'][~fullUDS['Stars-DR11']]
fulljmag = fullUDS['J'][~fullUDS['Stars-DR11']]


u_v = umag - vmag
v_j = vmag - jmag

chanu_v = chanumag - chanvmag
chanv_j = chanvmag - chanjmag

chanvaryu_v = chanvaryumag - chanvaryvmag
chanvaryv_j = chanvaryvmag - chanvaryjmag

fullu_v = fullumag - fullvmag
fullv_j = fullvmag - fulljmag

#%%
plt.figure(figsize=[8,8])
plt.plot(fullv_j, fullu_v, '.',markersize=0.5, color='tab:grey', alpha=0.35, label='UDS Galaxies')
plt.plot(chanv_j, chanu_v, '+', color='k', label='Non Variable Chandra Sources')
plt.plot(v_j, u_v,'bo', label='Non X-ray Variable Sources')
plt.plot(chanvaryv_j, chanvaryu_v,'ro', label='X-ray Variable Sources')

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
plt.legend()
plt.tight_layout()
plt.savefig('plots/Chi_variables/no06_extra_clean/uvjdiagram_Rmatched.png')
