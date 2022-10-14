#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 17:26:42 2022

Plot spectra from pypeit products

@author: 0cooper
"""

# the basics

import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from astropy.io import fits
from astropy.table import Table
from astropy.stats import sigma_clip
from scipy.signal import find_peaks
import read_mospy_files as rmp 
import glob
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve
from scipy.ndimage import gaussian_filter1d
plt.style.use('../cooper-presentation.mplstyle')

# select object
obj_name = sys.argv[1]
notes = Table.read('../cooper_cosmos_notes.csv',format='csv')
idx = np.where(notes['obj']==obj_name)[0]
zguess = float(notes['zguess'][idx])
lguess = float(notes['wave'][idx])
print('working on object',obj_name,'with zguess = ',zguess,'at lambda = ',lguess,'AA')

if notes['mask'][idx] == 'wmmc01':
    ddir = '../wmmc01/'
    fwhm = 0.456 # star_75
elif notes['mask'][idx] == 'wmmc02':
    ddir = '../wmmc02/pypeit_products/Science_coadd_1x1/'
    fwhm = 0.515 # star_0
elif notes['mask'][idx] == 'wmmc03':
    ddir = '../wmmc03/pypeit_products/Science_coadd_1x1/'
    fwhm = 0.648 #star_144
elif notes['mask'][idx] == 'wmmc05':
    ddir = '../wmmc05/pypeit_products/Science_coadd_1x1/'
    fwhm = 0.481 # star_575
elif notes['mask'][idx] == 'wmmc06':
    ddir = '../wmmc06/pypeit_products/Science_coadd_1x1/'
    fwhm = 0.548 # star_794
elif notes['mask'][idx] == 'wmmu01':
    ddir = '../wmmu01/pypeit_products/Science_coadd_1x1/'
    fwhm = 0.673 # star_62
elif notes['mask'][idx] == 'wmme01':
    ddir = '../wmme01/pypeit_products/Science_coadd_1x1/'
    fwhm = 0.697 #star_881
elif notes['mask'][idx] == 'wmme02':
    ddir = '../wmme02/pypeit_products/Science_coadd_1x1/'
    fwhm = 0.875 #star_513
elif notes['mask'][idx] == 'wmme03':
    ddir = '../wmme03/pypeit_products/Science_coadd_1x1/'
    fwhm = 0.558 #star_553
        
### add the other masks


# open 1D and 2D spectrum for an object
tab1 = Table.read(ddir+obj_name+'_1dspec.txt',format='ascii') # created 1d spec txt file from read_pypeit_files.py
file = glob.glob(ddir+'spec2d*.fits')[0] # 2d pypeit image
hdu2 = fits.open(file) # open image
pid = int(notes['pypeit_ref_idx'][idx])-1 # pypeit id to identify the slit
header = hdu2[0].header # header info
signal = hdu2[1].data # signal image data
wave = hdu2[8].data # wave image data
ref = hdu2[10].data[pid] # reference to slit in 2d image
print('check this is the right source:',ref)
x1 = ref[4][0] # left side of slit
x2 = ref[5][0] # right side of slit
xcen = ref[2] # object
sig2d = signal[0:,int(x1):int(x2)].T # signal for the slit (transposed)

pixscale = hdu2[0].header['PSCALE'] # arcsec/pix
nodamparc = 2.5 # arcsec
nodamp = nodamparc/pixscale # nod amplitude in pixels

# gaussian smooth 2D spec
kernel = Gaussian2DKernel(x_stddev=1.)
conv_im = convolve(sig2d, kernel)

### START FIG!
fig, (ax1, ax2) = plt.subplots(2, sharex=True)
fig.subplots_adjust(wspace=0,hspace=0)

# 1D spec
ax1.plot(tab1['pix'],gaussian_filter1d(tab1['opt_counts'],sigma=1),c='#C1CCEE',label=obj_name)
#ax1.plot(tab1['pix'],gaussian_filter1d(tab1['box_counts'],sigma=1),c='m',label=obj_name)
ax1.plot(tab1['pix'],tab1['opt_counts'],c='#C1CCEE',ls='--',alpha=0.6) 
#ax1.plot(tab1['pix'],tab1['box_counts'],c='m',ls='--',alpha=0.6) 
ax1.fill_between(tab1['pix'], y1=tab1['opt_sigma'], y2=-1*tab1['opt_sigma'], color='gray')
#ax1.fill_between(tab1['pix'], y1=tab1['box_sigma'], y2=-1*tab1['opt_sigma'], color='gray', alpha=0.2)
ax1.set_ylabel(r'Flux')

# 2D spec
lower=-1; upper=2
sample = sigma_clip(sig2d) 
vmin = sample.mean() + lower * sample.std()
vmax = sample.mean() + upper * sample.std()
ax2.imshow(conv_im, origin='lower', cmap='bone', aspect='auto', vmin=vmin, vmax=vmax)
ax2.set_xlabel('Column [Pix]')
ax2.set_ylabel('Row [Pix]')
ax2.grid(False)

# block out sky lines 2D
ipk = find_peaks(tab1['opt_sigma'], height=np.mean(tab1['opt_sigma']))[0]
for i in ipk:
    ax2.axvline(tab1['pix'][i],c='#4D4D6B',lw=15,alpha=0.85)

    
# secax = secondary axis for wavelength space
def forward(x):
    return np.interp(x, tab1['pix'], tab1['lambda'])

def inverse(x):
    return np.interp(x, tab1['lambda'], tab1['pix'])

secax = ax1.secondary_xaxis('top', functions=(forward, inverse))
secax.xaxis.set_minor_locator(AutoMinorLocator())
secax.set_xlabel(r'Wavelength [$\AA$]')

# mark measured line at pixel guess based on lambda guess
pguess = tab1['pix'][rmp.closest(tab1['lambda'],lguess)[0]] # pixel closest to wavelength of line to center on
s = tab1['opt_counts'][rmp.closest(tab1['lambda'],lguess)[0]] # signal closest to wavelength of line
n = tab1['opt_sigma'][rmp.closest(tab1['lambda'],lguess)[0]] # error closest to wavelength of line
sn = np.round(s/n,2)
print("candidate line at ",lguess,'AA; ',pguess,'pix')
print("detected at ",sn,' sigma')
ax1.axvline(pguess,lw=1,c='k',ls='-.',alpha=0.4)
ax2.axvline(pguess,lw=1,c='k',ls='-.',alpha=0.4)
ax2.axhline(xcen,lw=1,c='k',ls='-.',alpha=0.4)
ax2.axhline(xcen+nodamp,lw=1,c='k',ls=':',alpha=0.4)
ax2.axhline(xcen-nodamp,lw=1,c='k',ls=':',alpha=0.4)


ax1.set_xlim(pguess-70,pguess+70)
ax2.set_xlim(pguess-70,pguess+70)

ax1.set_ylim(-2*np.nanmedian(tab1['opt_sigma']),4*np.nanmedian(tab1['opt_sigma']))

plt.savefig('/Users/oc4858/Library/CloudStorage/Box-Box/werls/plots/'+obj_name+'_pypeit_1D_2D.png',dpi=500)
plt.show()
