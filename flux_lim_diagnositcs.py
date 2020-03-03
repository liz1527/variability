#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 15:06:17 2020

Code to diagnose where we should put a flux cut if including negatives

@author: ppxee
"""

### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
from astropy.table import Table #for handling tables
import numpy as np #for handling arrays
import vari_funcs #my module to help run code neatly
plt.close('all') #close any open plots

### Open the fits files and get data ###
Kvarys = Table.read('variable_tables/K/variables_no06_chi30_neg_DR11data.fits')
Jvarys = Table.read('variable_tables/J/K_extraction/J_variables_chi32_neg_DR11data_K.fits')
bothvarys = Table.read('variable_tables/variables_J_and_K_neg.fits')
Kvarysno = Table.read('variable_tables/K/variables_no06_chi30_DR11data.fits')
Jvarysno = Table.read('variable_tables/J/K_extraction/J_variables_chi32_noneg_DR11data.fits')
bothvarysno = Table.read('variable_tables/variables_J_and_K_noneg.fits')
deviant = Table.read('variable_tables/K/variables_no06_chi30_neg_deviant_DR11data.fits')
deviantno = Table.read('variable_tables/K/variables_no06_chi30_noneg_deviant.fits')

### Extract DR11 IDs from all tables ###
Kids = Kvarys['ID']
Jids = Jvarys['ID']
bothids = bothvarys['ID']
Kidsno = Kvarysno['ID']
Jidsno = Jvarysno['ID']
bothidsno = bothvarysno['ID']

### Restrict J and K to just those not in both ###
bothKdata = Kvarys[np.isin(Kids, bothids)]
bothKdatano = Kvarysno[np.isin(Kidsno, bothidsno)]

### Get fluxes ###
def get_flux(tbdata):
    flux = vari_funcs.k_mag_flux.flux4_stacks(tbdata)
    meanflux = np.nanmean(flux, axis=1)
    return meanflux

Kflux = get_flux(Kvarys)
Jflux = get_flux(Jvarys)
bothflux = get_flux(bothKdata)
Kfluxno = get_flux(Kvarysno)
Jfluxno = get_flux(Jvarysno)
bothfluxno = get_flux(bothKdatano)
devflux = get_flux(deviant)
devfluxno = get_flux(deviantno)

### Define bins ###
bins = np.array([13, 15])
bins = np.append(bins, np.arange(16,24,0.2))
bins = np.append(bins, [24])

bins = 10**((30-bins)/2.5)
bins = np.flip(bins, axis=0)

### get counts + plot hist ###
plt.figure()
Knum, _, _ = plt.hist(Kflux, bins, histtype='step', label='K-Band')
Jnum, _, _ = plt.hist(Jflux, bins, histtype='step', label='J-Band')
bothnum, _, _ = plt.hist(bothflux, bins, histtype='step', label='Both')
plt.xscale('log')
plt.xlabel('Mean K-band Flux')
plt.legend()
plt.tight_layout()
plt.figure()
Knumno, _, _ = plt.hist(Kfluxno, bins, histtype='step', label='K-Band')
Jnumno, _, _ = plt.hist(Jfluxno, bins, histtype='step', label='J-Band')
bothnumno, _, _ = plt.hist(bothfluxno, bins, histtype='step', label='Both')
plt.xscale('log')
plt.xlabel('Mean K-band Flux')
plt.legend()
plt.tight_layout()
plt.figure()
devnum, _, _ = plt.hist(devflux, bins, histtype='step')
devnumno, _, _ = plt.hist(devfluxno, bins, histtype='step')
plt.xscale('log')
plt.xlabel('Mean K-band Flux')
plt.legend()
plt.tight_layout()

### get fraction in each bin ###
totnum = (Knum + Jnum) - bothnum
frac = bothnum/totnum # fraction in both for neg
totnumno = (Knumno + Jnumno) - bothnumno
fracno = bothnumno/totnumno # fraction in both for no neg

devfrac = devnum/totnum # fraction deviant for neg
devfracno = devnumno/totnumno # fraction deviant for no neg

### plot frac as a function of magnitude ###
plt.figure()
x = bins[:len(bins)-1]
plt.plot(x, frac, label='No negatives')
plt.plot(x, fracno, label='Negatives')
plt.xscale('log')
plt.xlabel('Mean K-band Flux')
plt.ylabel('Fraction detected in J and K')
plt.xlim(8e1, 1e7)
plt.legend()
plt.tight_layout()

plt.figure()
plt.plot(x, devfracno, label='No negatives')
plt.plot(x, devfrac, label='Negatives')
plt.xscale('log')
plt.xlabel('Mean K-band Flux')
plt.ylabel('Fraction deviant')
plt.xlim(8e1, 1e7)
plt.legend()
plt.tight_layout()