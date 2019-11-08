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
import matplotlib.colors as colors
#from scipy.stats import chisquare
plt.close('all') #close any open plots
#plt.style.use('default')

### Get data ###
nontb = fits.open('variable_tables/no06_variables_chi30_2arcsec_nochanXray_DR11data_restframe.fits')[1].data
xraytb = fits.open('variable_tables/no06_variables_chi30_2arcsec_chandata_DR11data_restframe.fits')[1].data
#nontb = Table.read('variable_tables/no06_variables_chi30_2arcsec_noXray_DR11data_restframe.fits')
#xraytb = Table.read('variable_tables/no06_variables_chi30_2arcsec_Xray_DR11data_restframe.fits')
sc = Table.read('variable_tables/no06_variables_chi30_2arcsec_DR11data_restframe_SC.fits')
fullxray = Table.read('mag_flux_tables/K/novarys_chanDR11data_restframe_mag_flux_table_best_extra_clean_no06.fits')
fullUDS = Table.read('UDS_catalogues/DR11-2arcsec-June24-2018+plusXY_UVJ-SC.fits')

### Limit to Chandra region for simplicity ###
nontb = vari_funcs.chandra_only(nontb)
xraytb = vari_funcs.chandra_only(xraytb)
fullxray = vari_funcs.chandra_only(fullxray)
#chandata = vari_funcs.chandra_only(chandata)
#fullUDS = vari_funcs.chandra_only(fullUDS)

#
#### Plot R band distributions so can match for proposal ###
#xray_R = xraytb['RMAG_20']
#non_R = nontb['RMAG_20']
#plt.hist([xray_R, non_R], color=['r','b'], histtype='step', 
#         label = ['X-ray','Non X-ray'])
#plt.xlabel('R Band Magnitude')
#plt.ylabel('Counts')
#plt.legend()
#plt.tight_layout()

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

#plt.figure()
#plt.hist([diffu,diffv,diffj],50)


#xraytb = vari_funcs.flux_split(xraytb, 'upper')
#nontb = vari_funcs.flux_split(nontb, 'upper')



allxu_v, allxv_j = get_data(fullxray, meanu, meanv, meanj)
noxu_v, noxv_j = get_data(nontb, meanu, meanv, meanj)
xu_v, xv_j = get_data(xraytb, meanu, meanv, meanj)

### Use actual SC for UDS sources ###
fullumag = fullUDS['U'][~fullUDS['Stars-DR11']]
fullvmag = fullUDS['V'][~fullUDS['Stars-DR11']]
fulljmag = fullUDS['J'][~fullUDS['Stars-DR11']]
fullu_v = fullumag - fullvmag
fullv_j = fullvmag - fulljmag

### get flux for colours ###
noxflux = vari_funcs.flux4_stacks(nontb)
xflux = vari_funcs.flux4_stacks(xraytb)

noxmean = np.nanmean(noxflux, axis=1)
xmean = np.nanmean(xflux, axis=1)

#noxmean = vari_funcs.get_jansky_flux(nontb)
#xmean = vari_funcs.get_jansky_flux(xraytb)

### find max and min ###
cmax = np.nanmax([np.nanmax(noxmean), np.nanmax(xmean)])
cmin = np.nanmin([np.nanmin(noxmean), np.nanmin(xmean)])

#%%
plt.figure(figsize=[8,8])
plt.plot(fullv_j, fullu_v, '.',markersize=0.5, color='tab:grey', alpha=0.2, label='Galaxy', zorder=0)
plt.plot(allxv_j, allxu_v, 's', markersize=5, color='k', label='Non-Variable X-ray AGN', zorder=1)
plt.plot(noxv_j, noxu_v,'bo', label='Variable Non X-ray AGN', zorder=3)
plt.plot(xv_j, xu_v,'ro', label='Variable X-ray AGN', zorder=3)
plt.legend()

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
plt.tight_layout()

#%%
plt.figure(figsize=[9,8])
plt.plot(fullv_j, fullu_v, '.',markersize=0.5, color='tab:grey', alpha=0.35, label='UDS Galaxies', zorder=0)
#plt.plot(allxv_j, allxu_v, '+', color='k', label='Non Variable Chandra Sources', zorder=1)
#plt.plot(noxv_j, noxu_v,'bo', mfc='None', markersize=7,
#         label='Non X-ray Variable Sources', zorder=3)
#plt.plot(xv_j, xu_v,'ro', mfc='None', markersize=7,
#         label='X-ray Variable Sources', zorder=3)
plt.legend()

plt.scatter(noxv_j, noxu_v, c=noxmean, marker='x',
            norm=colors.LogNorm(vmin=cmin, vmax=cmax), zorder=2, label='Non X-ray Variable Sources')
plt.scatter(xv_j, xu_v, c=xmean, marker='o', 
            norm=colors.LogNorm(vmin=cmin, vmax=cmax), zorder=2, label='X-ray Variable Sources')
cbar = plt.colorbar()
cbar.set_label('Mean 2" Flux')
#cbar.set_label('2" Flux (jy)')
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
plt.tight_layout()
#plt.savefig('plots/Chi_variables/no06_extra_clean/uvjdiagram_Rmatched.png')


#%%
### get mstar for colours ###
noxm = nontb['Mstar_z_p']
xm = xraytb['Mstar_z_p']
allxm = fullxray['Mstar_z_p']

### ignore 0s ###
def ignore_zeros(m, v_j, u_v):
    m[m==0] = np.nan
    mask = ~np.isnan(m)
    v_j = v_j[mask]
    u_v = u_v[mask]
    m = m[mask]
    return m, v_j, u_v
noxm, noxv_j, noxu_v = ignore_zeros(noxm, noxv_j, noxu_v)
xm, xv_j, xu_v = ignore_zeros(xm, xv_j, xu_v)
#allxm, allxv_j, allxu_v = ignore_zeros(allxm, allxv_j, allxu_v)

### plot hist for info ###
plt.figure()
plt.hist([noxm, xm], bins=np.logspace(4.5,12), color=['b','r'], histtype='step', label=['Non X-ray','X-ray'])
#plt.hist([noxm, xm, allxm], bins=np.logspace(4.5,12), color=['b','r','k'], histtype='step')
plt.xscale('log')
plt.xlabel(r'$M_{star}$')
plt.ylabel('Number')
plt.legend(loc='upper left')
plt.tight_layout()


### find max and min ###
cmax = np.nanmax([np.nanmax(noxm), np.nanmax(xm)])
cmin = 1e8#np.nanmin([np.nanmin(noxm), np.nanmin(xm)])
plt.figure(figsize=[9,8])
plt.plot(fullv_j, fullu_v, '.',markersize=0.5, color='tab:grey', alpha=0.35, label='UDS Galaxies', zorder=0)
plt.plot(allxv_j, allxu_v, 's', color='k', markersize=5, label='Non Variable Chandra Sources', zorder=1)
#plt.plot(noxv_j, noxu_v,'bo', mfc='None', markersize=7,
#         label='Non X-ray Variable Sources', zorder=3)
#plt.plot(xv_j, xu_v,'ro', mfc='None', markersize=7,
#         label='X-ray Variable Sources', zorder=3)
plt.legend()

plt.scatter(noxv_j, noxu_v, c=noxm, marker='x',
            norm=colors.LogNorm(vmin=cmin, vmax=cmax), zorder=2, label='Non X-ray Variable Sources')
plt.scatter(xv_j, xu_v, c=xm, marker='o',
            norm=colors.LogNorm(vmin=cmin, vmax=cmax), zorder=2, label='X-ray Variable Sources')
cbar = plt.colorbar()
cbar.set_label(r'$M_{star}$')
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
plt.tight_layout()