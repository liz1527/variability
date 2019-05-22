#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 10:59:07 2018

Code to explore properties of variable sources

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

### Get data ###
datatb = Table.read('variable_tables/no06_variables_chi30_edged_DR11data_restframe.fits')
fullxray = Table.read('mag_flux_tables/chanDR11data_restframe_mag_flux_table_best_extra_clean_no06.fits')
chandata = Table.read('variable_tables/no06_variables_chi30_chandata_DR11data_restframe.fits')
fullUDS = Table.read('UDS_catalogues/DR11-2arcsec-June24-2018+plusXY.fits')

### Remove Edges ###
#datatb = vari_funcs.remove_edges(datatb)
#fullxray = vari_funcs.remove_edges(fullxray)
#chandata = vari_funcs.remove_edges(chandata)
#fullUDS = vari_funcs.remove_edges(fullUDS)

### Limit to Chandra region for simplicity ###
datatb = vari_funcs.chandra_only(datatb)
fullxray = vari_funcs.chandra_only(fullxray)
chandata = vari_funcs.chandra_only(chandata)
#fullUDS = vari_funcs.chandra_only(fullUDS)

### Split the variable table into x-ray and non x-ray ###
xraytb = datatb[datatb['X-ray']]
nontb = datatb[~datatb['X-ray']]

### Get stuff from table ###
mask1 = datatb['RMAG_20']!=99
mask2 = datatb['KMAG_20']!=99
allmask = mask1*mask2.astype(bool)
rmag = datatb['RMAG_20'][allmask]
kmag = datatb['KMAG_20'][allmask]

mask1 = nontb['RMAG_20']!=99
mask2 = nontb['KMAG_20']!=99
nonmask = mask1*mask2.astype(bool)
nonrmag = nontb['RMAG_20'][nonmask]
nonkmag = nontb['KMAG_20'][nonmask]

mask1 = xraytb['RMAG_20']!=99
mask2 = xraytb['KMAG_20']!=99
xraymask = mask1*mask2.astype(bool)
xrayrmag = xraytb['RMAG_20'][xraymask]
xraykmag = xraytb['KMAG_20'][xraymask]

mask1 = fullxray['RMAG_20']!=99
mask2 = fullxray['KMAG_20']!=99
fullmask = mask1*mask2.astype(bool)
fullrmag = fullxray['RMAG_20'][fullmask]
fullkmag = fullxray['KMAG_20'][fullmask]

mask1 = fullUDS['RMAG_20']!=99
mask2 = fullUDS['KMAG_20']!=99
udsmask = mask1*mask2.astype(bool)
udsrmag = fullUDS['RMAG_20'][udsmask]
udskmag = fullUDS['KMAG_20'][udsmask]
##%%
#plt.figure()
#plt.hist(udskmag,50)
##%%
#specz = datatb['z_spec'][allmask]
#xrayspecz = xraytb['z_spec'][xraymask]
#nonspecz = nontb['z_spec'][nonmask]
#z = datatb['z_p'][allmask]
#xrayz = xraytb['z_p'][xraymask]
#nonz = nontb['z_p'][nonmask]
#fullz = fullxray['z_p'][fullmask]
#
### R-K against K plot ###
plt.figure(figsize=[7,7])
nonrmink = nonrmag-nonkmag
xrayrmink = xrayrmag - xraykmag
fullrmink = fullrmag - fullkmag
udsrmink = udsrmag - udskmag

plt.plot(udskmag, udsrmink,'.', markersize=0.1, #alpha=0.35,
         color='tab:gray',label='Non Variable UDS Sources')
plt.plot(fullkmag, fullrmink,'k+',label='Non Variable X-ray Sources')
plt.plot(nonkmag, nonrmink, 'bo',label='Non X-ray Variable Sources')
plt.plot(xraykmag, xrayrmink, 'ro',label= 'X-ray Variable Sources')
plt.xlabel('K')
plt.ylabel('R-K')
plt.legend()
plt.tight_layout()
#plt.savefig('plots/Chi_variables/no06_extra_clean/RminusKvsK.png')
#
#plt.figure(figsize=[7,7])
##plt.plot(udskmag, udsrmink,'s', color='tab:gray',label='Non Variable UDS Sources')
#plt.plot(fullz, fullrmink,'o', color='tab:gray',label='Non Variable X-ray Sources')
#plt.plot(nonz, nonrmink, 'bo',label='Non X-ray Variable Sources')
#plt.plot(xrayz, xrayrmink, 'ro',label= 'X-ray Variable Sources')
#plt.xlabel('z')
#plt.ylabel('R-K')
#plt.legend()
#plt.tight_layout()
#plt.savefig('plots/Chi_variables/no06_extra_clean/RminusKvsz.png')
#
#### Get histogram of magnitudes ###
plt.figure()
_, rbins, _ = plt.hist(rmag, 25, range=(19,30))
plt.xlabel('R Band Magnitude')
plt.ylabel('Counts')
plt.tight_layout()
#plt.savefig('plots/Chi_variables/no06_extra_clean/Rbandmaghist.png')
#
#plt.figure()
#_, kbins, _ = plt.hist(kmag, 25, range=(17,26))
#plt.xlabel('K Band Magnitude')
#plt.ylabel('Counts')
#plt.tight_layout()
#plt.savefig('plots/Chi_variables/no06_extra_clean/Kbandmaghist.png')
#
### plot historgrams with xray and non xray separated ###
plt.figure()
plt.hist(nonrmag, rbins, color='b', label = 'Non X-ray',histtype='step')
plt.hist(xrayrmag, rbins, color='r', label='X-ray',histtype='step')
plt.xlabel('R Band Magnitude')
plt.ylabel('Counts')
plt.legend()
plt.tight_layout()
#plt.savefig('plots/Chi_variables/no06_extra_clean/Rbandmaghistcomp.pdf')
#
##
##plt.figure()
##plt.hist([nonrmag, xrayrmag], rbins, color=['b','r'], 
##         label = ['Non X-ray','X-ray'], histtype='barstacked')
##plt.xlabel('R Band Magnitude')
##plt.ylabel('Counts')
##plt.legend()
##plt.tight_layout()
###plt.savefig('50rbandmaghistcompstack.png')
#
#plt.figure()
#plt.hist(nonkmag, kbins, color='b', label = 'Non X-ray',histtype='step')
#plt.hist(xraykmag, kbins, color='r', label='X-ray',histtype='step')
#plt.xlabel('K Band Magnitude')
#plt.ylabel('Counts')
#plt.legend()
#plt.tight_layout()
#plt.savefig('plots/Chi_variables/no06_extra_clean/Kbandmaghistcomp.png')
#
##plt.figure()
##plt.hist([nonkmag, xraykmag], kbins, color=['b','r'], 
##         label = ['Non X-ray','X-ray'], histtype='barstacked')
##plt.xlabel('K Band Magnitude')
##plt.ylabel('Counts')
##plt.legend()
##plt.tight_layout()
###plt.savefig('50kbandmaghistcompstack.png')
#
#### Plot histogram of spec z ###
#plt.figure()
#_, zbins, _ = plt.hist(specz, bins=25, range=(0,np.nanmax(specz)))
#plt.xlabel('z_spec')
#plt.ylabel('Counts')
#plt.savefig('plots/Chi_variables/no06_extra_clean/speczhist.png')
#
#plt.figure()
#plt.hist(nonspecz, zbins, color='b', label = 'Non X-ray',histtype='step')
#plt.hist(xrayspecz, zbins, color='r', label='X-ray',histtype='step')
#plt.xlabel('z_spec')
#plt.ylabel('Counts')
#plt.savefig('plots/Chi_variables/no06_extra_clean/speczmaghistcomp.png')
#
#### Plot R/K band mag against X-ray flux ###
##xmmdata = Table.read('variable_tables/variables_chi40_xmmdata.fits')
#chanfull = fullxray['Full_flux']
#chanhard = fullxray['Hard_flux']
#chanrmag = fullxray['RMAG_20']
#chankmag = fullxray['KMAG_20']
#chanvaryfull = chandata['Full_flux']
#chanvaryhard = chandata['Hard_flux']
#chanvaryrmag = chandata['RMAG_20']
#chanvarykmag = chandata['KMAG_20']
#
#### remove any with rmag==99 ###
#mask = chanrmag==99
#chanfull = chanfull[~mask]
#chanhard = chanhard[~mask]
#chanrmag = chanrmag[~mask]
#chankmag = chankmag[~mask]
#mask = chanvaryrmag==99
#chanvaryfull = chanvaryfull[~mask]
#chanvaryhard = chanvaryhard[~mask]
#chanvaryrmag = chanvaryrmag[~mask]
#chanvarykmag = chanvarykmag[~mask]
#
#### have plotted both varying and non varying X-ray, now need to plot non x-ray
#### with x-ray upper limits from Chandra ###
#nonfull = np.full_like(nonrmag, 4.4e-16)
#nonhard = np.full_like(nonrmag, 6.5e-16)
#
#### Add fx/fv ~ 1 line to see if they lie on it ###
#fv = datatb['RFLUX_20']
#fx = 1e-11*fv
#cx = fv
#
#### Actually make the plot
#plt.figure()
#plt.plot(chanrmag, chanfull, 'o',color='tab:grey', label='X-ray Source')
#plt.plot(chanvaryrmag, chanvaryfull, 'ro', label='Variable X-ray Source')#, mfc='None',markersize=10)
#plt.errorbar(nonrmag, nonfull,yerr=1e-16, fmt='bo', #mfc='None', #markersize=10,
#             uplims=True, label='Non X-ray Source')
##plt.plot(rmag, fx)
#plt.yscale('log')
#plt.xlabel('R band magnitude')
#plt.ylabel('Full Chandra Flux (erg/s/cm2)')
##plt.ylabel('Full Chandra Flux (ergs/s/cm^2)')
#plt.legend()
#plt.tight_layout()
#plt.savefig('plots/Chi_variables/no06_extra_clean/RvsFullXray.png')
#
#
#### Actually make the plot
#plt.figure()
#plt.plot(chanrmag, chanhard, 'o',color='tab:grey', label='X-ray Source')
#plt.plot(chanvaryrmag, chanvaryhard, 'ro', label='Variable X-ray Source')#, mfc='None',markersize=10)
#plt.errorbar(nonrmag, nonhard,yerr=1e-16, fmt='bo', #mfc='None', #markersize=10,
#             uplims=True, label='Non X-ray Source')
##plt.plot(rmag, fx)
#plt.yscale('log')
#plt.xlabel('R band magnitude')
#plt.ylabel('Hard Chandra Flux (erg/s/cm2)')
##plt.ylabel('Full Chandra Flux (ergs/s/cm^2)')
#plt.legend()
#plt.tight_layout()
#plt.savefig('plots/Chi_variables/no06_extra_clean/RvsHardXray.png')
#
#### Actually make the k plot
#plt.figure()
#plt.plot(chankmag, chanfull, 'o',color='tab:grey', label='X-ray Source')
#plt.plot(chanvarykmag, chanvaryfull, 'ro', label='Variable X-ray Source')#, mfc='None',markersize=10)
#plt.errorbar(nonkmag, nonfull,yerr=1e-16, fmt='bo', #mfc='None', #markersize=10,
#             uplims=True, label='Non X-ray Source')
##plt.plot(rmag, fx)
#plt.yscale('log')
#plt.xlabel('K band magnitude')
#plt.ylabel('Full Chandra Flux (erg/s/cm2)')
##plt.ylabel('Full Chandra Flux (ergs/s/cm^2)')
#plt.legend()
#plt.tight_layout()
#plt.savefig('plots/Chi_variables/no06_extra_clean/KvsFullXray.png')
#
#### Actually make the k plot
#plt.figure()
#plt.plot(chankmag, chanhard, 'o',color='tab:grey', label='X-ray Source')
#plt.plot(chanvarykmag, chanvaryhard, 'ro', label='Variable X-ray Source')#, mfc='None',markersize=10)
#plt.errorbar(nonkmag, nonhard,yerr=1e-16, fmt='bo', #mfc='None', #markersize=10,
#             uplims=True, label='Non X-ray Source')
##plt.plot(rmag, fx)
#plt.yscale('log')
#plt.xlabel('K band magnitude')
#plt.ylabel('Hard Chandra Flux (erg/s/cm2)')
##plt.ylabel('Full Chandra Flux (ergs/s/cm^2)')
#plt.legend()
#plt.tight_layout()
#plt.savefig('plots/Chi_variables/no06_extra_clean/KvsHardXray.png')

#### Make UVJ diagram ###
#umag = datatb['M_u_z_p']
#vmag = datatb['M_V_z_p']
#jmag = datatb['M_J_z_p']
#chanumag = fullxray['M_u_z_p']
#chanvmag = fullxray['M_V_z_p']
#chanjmag = fullxray['M_J_z_p']
#chanvaryumag = chandata['M_u_z_p']
#chanvaryvmag = chandata['M_V_z_p']
#chanvaryjmag = chandata['M_J_z_p']
#fullumag = fullUDS['M_u_z_p'][~fullUDS['Stars-DR11']]
#fullvmag = fullUDS['M_V_z_p'][~fullUDS['Stars-DR11']]
#fulljmag = fullUDS['M_J_z_p'][~fullUDS['Stars-DR11']]
#
#u_v = umag - vmag
#v_j = vmag - jmag
#
#chanu_v = chanumag - chanvmag
#chanv_j = chanvmag - chanjmag
#
#chanvaryu_v = chanvaryumag - chanvaryvmag
#chanvaryv_j = chanvaryvmag - chanvaryjmag
#
#fullu_v = fullumag - fullvmag
#fullv_j = fullvmag - fulljmag
##%%
#plt.figure(figsize=[7,7])
#plt.plot(fullv_j, fullu_v, '.',markersize=0.5, color='tab:grey', alpha=0.35, label='UDS Galaxies')
#plt.plot(chanv_j, chanu_v, '+', color='k', label='Non Variable Chandra Sources')
#plt.plot(v_j, u_v,'bo', label='Non X-ray Variable Sources')
#plt.plot(chanvaryv_j, chanvaryu_v,'ro', label='X-ray Variable Sources')
#
### plot dividing lines => y = 0.875x + 0.6 #
##x = np.linspace(0.8,1.6)
##y = 0.875*x + 0.6
##plt.plot(x,y,'k--')
##plt.hlines(1.3,-0.7,0.8,linestyles='dashed')
##plt.vlines(1.6,2.0,2.1,linestyles='dashed')
#
#plt.xlim(xmin=-0.7, xmax=3)
#plt.ylim(ymin=-1,ymax=2.1)
#plt.xlabel('V - J')
#plt.ylabel('U - V')
#plt.legend()
#plt.tight_layout()
##plt.savefig('plots/Chi_variables/no06_extra_clean/uvjdiagram.png')
