#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 12:39:35 2022

@author: 0cooper
"""

# the basics
import sys
import numpy as np 
from matplotlib import pyplot as plt 
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import scipy.io
from astropy import units as u 
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.cosmology import Planck15 as cosmo
import warnings
warnings.filterwarnings("ignore")
plt.style.use('../cooper-presentation.mplstyle')


# select object
obj_name = sys.argv[1]
notes = Table.read('../cooper_cosmos_notes.csv',format='csv')
idx = np.where(notes['obj']==obj_name)[0]
zguess = float(notes['zguess'][idx])
wra = float(notes['objra'][idx]); wdec = float(notes['objdec'][idx])
print('working on object',obj_name,'with zguess = ',zguess)


if notes['mask'][idx][0][:-2] == 'wmmc':
    file = '/Users/oc4858/werls/COSMOS_F21.idl'
    idl = scipy.io.readsav(file).cat
    # load best fit sed
    wvl = idl.photz[0].templates[0].lam_rest[0]*u.AA
    fnu_mods = idl.photz[0].templates[0].model[0]*u.erg/u.s/u.Hz/u.cm**2
elif notes['mask'][idx][0][:-2] == 'wmmu':
    #file = '/Users/oc4858/Library/CloudStorage/Box-Box/treasurechest/uds/UDS_candidates_all.idl'
    file = '/Users/oc4858/werls/UDS_candidates_all.idl'
    idl = scipy.io.readsav(file).fobj
    # load best fit sed
    wvl = idl.photz[0].lam_rest_model[0]*u.AA
    fnu_mods = idl.photz[0].model[0]*u.erg/u.s/u.Hz/u.cm**2
elif notes['mask'][idx][0][:-2] == 'wmme':
    file = '/Users/oc4858/Library/CloudStorage/Box-Box/treasurechest/egs/EGS_F22_trimmed.idl'
    pzfile = '/Users/oc4858/Library/CloudStorage/Box-Box/treasurechest/egs/EGS_F22_pz.idl'
    ezfile = '/Users/oc4858/Library/CloudStorage/Box-Box/treasurechest/egs/EGS_F22_EAZYSED.idl'
    pz = scipy.io.readsav(pzfile)
    ez = scipy.io.readsav(ezfile)
    idl = scipy.io.readsav(file).cat
    # load best fit sed
    wvl = ez.lam_rest*u.AA
    fnu_mods = ez.model*u.erg/u.s/u.Hz/u.cm**2

# load catalog ra/dec
catra = idl.ra[0]; catdec = idl.dec[0]
# crossmatch
c = SkyCoord(ra=catra*u.degree, dec=catdec*u.degree)
catalog = SkyCoord(ra=wra*u.degree, dec=wdec*u.degree)
ind, d2d, d3d = catalog.match_to_catalog_sky(c)
if d2d > 1*u.arcsec:
    print("!! xmatch to source is > 1 arcsec !!")

print('distance to xmatch is ',d2d.to(u.arcsec))

# load fluxes
fnu_obs = [idl.flux[0][0][i][ind] for i in range(len(idl.flux[0][0]))]*u.nJy
fnu_obserr = [idl.dflux[0][0][i][ind] for i in range(len(idl.dflux[0][0]))]*u.nJy
filt_name = idl.flux[0].dtype.names
filt_lambda = [.592188,.804553,1.055025,1.248607,1.392321,1.537034,3.537841,4.478049]*u.um
filt_fwhm = [.232293,.185835,.291703,.300520,.394088,.287418,.743171,1.009682]*u.um

# find mags
Jmag = (idl.flux[0].FJ[0][ind]*u.nJy).to(u.ABmag) # apparent J mag
MUV = (Jmag.value)-(cosmo.distmod(zguess).value)+(2.5*np.log10(1+zguess)) # abs UV mag

# load best fit sed
fnu_mod = fnu_mods[ind]

if notes['mask'][idx][0][:-2] == 'wmme':
    # load zpdf
    za = ez.za[ind]
    zarr = pz.zgrid
    zpdf = pz.pz[ind]
else:
    # load zpdf
    za = idl.photz[0].za[0][ind]
    zarr = idl.photz[0].zgrid[0]
    zpdf = idl.photz[0].pz[0][ind]


### Begin figure!
fig, ax = plt.subplots(1, figsize=[7,5])

ax.set_ylabel('Flux [$10^{-29}$ erg/cm$^2$/s/Hz]')
ax.errorbar(filt_lambda.to(u.AA),fnu_obs.to(u.erg/u.s/u.Hz/u.cm**2)*1e29,xerr=filt_fwhm.to(u.AA),yerr=fnu_obserr.to(u.erg/u.s/u.Hz/u.cm**2)*1e29,fmt='.w', ecolor = 'w', capsize=3, elinewidth=1,zorder=2)
ymin = (np.min(fnu_obs.to(u.erg/u.s/u.Hz/u.cm**2))*1e29*0.5).value
ymax = (np.max(fnu_obs.to(u.erg/u.s/u.Hz/u.cm**2))*1e29*10).value


# Plot the best-fit model 
ax.plot(wvl*(1+za),fnu_mod*1e29,color='c',alpha=1,label='model',zorder=1) 

# Show where nebular emission lines would potentially boost the flux
ax.vlines(3727*(1+za),ymin,ymax,label='[OII]',zorder=0,color='0.3',ls=':')
ax.vlines(5007*(1+za),ymin,ymax,label='[OIII]b',zorder=0,color='0.3',ls=':')
ax.vlines(4861*(1+za),ymin,ymax,label='Hb',zorder=0,color='0.3',ls=':') # H_beta
ax.vlines(6563*(1+za),ymin,ymax,label='Ha',zorder=0,color='0.3',ls=':') # H_alpha 
ax.vlines(1216*(1+zguess),ymin,ymax,label='Lyman-alpha break',zorder=0,color='chocolate',ls='-.')

# set plot limits
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(1000,300000) 
ymin1 = [ymin if (ymin>1e-5)==True else 1e-5][0]
ax.set_ylim(ymin1,ymax)
ax.set_xlabel('wavelength [Å]')

# annotate
ax.annotate(text='photo-z='+str(np.round(za,2)),xy=(0.15,0.83),xycoords='figure fraction',fontsize='small',color='w')
ax.annotate(text='spec-z='+str(zguess),xy=(0.15,0.8),xycoords='figure fraction',fontsize='small',color='w')
ax.annotate(text='MUV='+str(np.round(MUV,2)),xy=(0.15,0.77),xycoords='figure fraction',fontsize='small',color='w')
ax.annotate(text='J mag='+str(np.round(Jmag.value,2)),xy=(0.15,0.74),xycoords='figure fraction',fontsize='small',color='w')


# plot zpdf
# Create inset of width 30% and height 40% of the parent axes' bounding box
axin = inset_axes(ax, width="33%", height="33%", loc='lower right')
axin.plot(zarr,zpdf,c='teal',lw=3)
axin.axvline(zguess,c='chocolate',lw=3)
xmin = [za-1.5 if (za > 1.5)==True else 0][0]
axin.set_xlim(0,zguess+1.5)
axin.xaxis.tick_top()
axin.tick_params(size=10)
axin.yaxis.set_ticks([])
axin.xaxis.set_label_position('top') 
axin.grid(False)


plt.savefig('/Users/oc4858/Library/CloudStorage/Box-Box/werls/plots/'+obj_name+'_candels_sed.png',dpi=500)


