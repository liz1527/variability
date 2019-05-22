#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  5 16:24:18 2018

Code to show error differences between correction methods

@author: ppxee
"""

import numpy as np
from astropy.table import Table
import vari_funcs
import matplotlib.pyplot as plt
plt.close('all')

tbdata = Table.read('mag_flux_tables/mag_flux_table_best_corr.fits')

#flux = vari_funcs.flux5_stacks(tbdata)
#fluxerr = vari_funcs.fluxerr5_stacks(tbdata)
#newfluxerr = vari_funcs.fluxerr5_stacks_corr(tbdata)

#ob = 67178

obnums = [67178, 227399, 203391, 187967,193325,243208]

for n, ob in enumerate(obnums):
    plt.figure(1)
    plt.subplot(2,3,n+1)
    vari_funcs.lightcurveflux5(int(ob), tbdata, corrected=True, new_fig=False)
#    plt.tight_layout()
    axes = plt.gca()
    ylims = axes.get_ylim()
#    ymid = (ylims[1]+ylims[0])/2
#    plt.ylim(ymin=ymid-5000, ymax=ymid+5000)    
    plt.figure(2)
    plt.subplot(2,3,n+1)
    plt.savefig(str(ob)+'correctedfluxcurve.png')
    vari_funcs.lightcurveflux5(int(ob), tbdata, new_fig=False)
#    plt.tight_layout()
#    axes = plt.gca()
#    ylims = axes.get_ylim()
#    ymid = (ylims[1]+ylims[0])/2
    plt.ylim(ymin=ylims[0], ymax=ylims[1])
    plt.savefig(str(ob)+'fluxcurve.png')
