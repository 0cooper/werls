#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 11:39:36 2022

@author: 0cooper
"""

# the basics
import sys
import numpy as np 
from matplotlib import pyplot as plt 
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import h5py # used in the Data Visualization section 
from astropy.coordinates import SkyCoord
import astropy
from astropy.io import fits,ascii,votable
from astropy import units as u 
from astropy import constants as const
from astropy.table import Table
from astropy.cosmology import Planck15 as cosmo

plt.style.use('../cooper-presentation.mplstyle')


# Specify the version of the catalog and the folder with the input/output files
catversion = 'Farmer'  # this string can be either 'Classic' or 'Farmer'
dir_in = '/Users/oc4858/Library/CloudStorage/Box-Box/treasurechest/cosmos/COSMOS2020_R1/'  
#dir_out = './'  # the directory where the output of this notebook will be stored
fitversion = 'lp'  
# Which type of photometric estimates to use? (suffix of the column name)
# This choice must be consistent with `catversion`,
# choices for Classic are: '_FLUX_APER2', '_FLUX_APER3', '_MAG_APER2,', '_MAG_APER3'
# choices for Farmer are '_FLUX' or '_MAG' 
flx = '_FLUX'  
flxerr = '_FLUXERR'  # catalog column for flux/mag error, just add 'ERR'
outflx = 'cgs' # 'cgs' or 'uJy'
# Filter names, mean wavelength, and other info (see Table 1 in W+21)
filt_name = ['GALEX_FUV', 'GALEX_NUV','CFHT_u','CFHT_ustar','HSC_g', 'HSC_r', 'HSC_i', 'HSC_z', 'HSC_y', 'UVISTA_Y', 'UVISTA_J', 'UVISTA_H', 'UVISTA_Ks', 'SC_IB427', 'SC_IB464', 'SC_IA484', 'SC_IB505', 'SC_IA527', 'SC_IB574', 'SC_IA624', 'SC_IA679', 'SC_IB709', 'SC_IA738', 'SC_IA767', 'SC_IB827', 'SC_NB711', 'SC_NB816', 'UVISTA_NB118', 'SC_B', 'SC_gp', 'SC_V', 'SC_rp', 'SC_ip','SC_zp', 'SC_zpp', 'IRAC_CH1', 'IRAC_CH2', 'IRAC_CH3','IRAC_CH4']  
filt_lambda = [0.1526,0.2307,0.3709,0.3858,0.4847,0.6219,0.7699,0.8894,0.9761,1.0216,1.2525,1.6466,2.1557,0.4266,0.4635,0.4851,0.5064,0.5261,0.5766,0.6232,0.6780,0.7073,0.7361,0.7694,0.8243,0.7121,0.8150,1.1909,0.4488,0.4804,0.5487,0.6305,0.7693,0.8978,0.9063,3.5686,4.5067,5.7788,7.9958]
filt_fwhm = [0.0224,0.07909,0.05181,0.05976,0.1383,0.1547,0.1471,0.0766,0.0786,0.0923,0.1718,0.2905,0.3074,0.02073,0.02182,0.02292,0.0231,0.02429,0.02729,0.03004,0.03363,0.03163,0.03235,0.03648,0.0343,0.0072,0.01198,0.01122,0.0892,0.1265,0.0954,0.1376,0.1497,0.0847,0.1335,0.7443,1.0119,1.4082,2.8796] 
# corresponding MW attenuation from Schelgel 
AlambdaDivEBV = [8.31,8.742,4.807,4.674,3.69,2.715,2.0,1.515,1.298,1.213,0.874,0.565,0.365,4.261,3.844,3.622,3.425,3.265,2.938,2.694,2.431,2.29,2.151,1.997,1.748,2.268,1.787,0.946,4.041,3.738,3.128,2.673,2.003,1.436,1.466,0.163,0.112,0.075,0.045]
# photometric offsets (not available for all filters, see Table 3 in W+21)
zpoff1 = [0.000,-0.352,-0.077,-0.023,0.073,0.101,0.038,0.036,0.086,0.054,0.017,-0.045,0.000,-0.104,-0.044,-0.021,-0.018,-0.045,-0.084,0.005,0.166,-0.023,-0.034,-0.032,-0.069,-0.010,-0.064,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,-0.212,-0.219,0.000,0.000]  # Farmer+LePhare
zpoff2 = [0.000,-0.029,-0.006,0.053,0.128,0.127,0.094,0.084,0.100,0.049,0.025,-0.044,0.000,-0.013,-0.008,0.022,0.025,0.033,-0.032,0.031,0.208,-0.009,0.003,-0.015,-0.001,0.023,-0.021,-0.017,-0.075,0.000,0.123,0.035,0.051,0.000,0.095,-0.087,-0.111,0.000,0.000]  # Classic+LePhare
zpoff3 = [0.000,0.000,-0.196,-0.054,0.006,0.090,0.043,0.071,0.118,0.078,0.047,-0.034,0.000,-0.199,-0.129,-0.084,-0.073,-0.087,-0.124,0.004,0.154,-0.022,-0.030,-0.013,-0.057,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,-0.102,-0.044,0.000,0.000] # Farmer+EAZY
zpoff4 = [0.000,0.000,0.000,-0.021,0.055,0.124,0.121,0.121,0.145,0.085,0.057,-0.036,0.000,-0.133,-0.098,-0.046,-0.037,-0.038,-0.062,0.038,0.214,0.024,0.022,0.01,0.022,0.000,0.000,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.021,0.025,0.000,0.000] # Classic+EAZY
# create the dictionary
filt_dict = {filt_name[i]:(filt_lambda[i]*1e4,filt_fwhm[i]*1e4,AlambdaDivEBV[i],[zpoff1[i],zpoff2[i],zpoff3[i],zpoff4[i]]) for i in range(len(filt_name))}
# read in sed template file
hdf = h5py.File(dir_in+"../cosmos2020-readcat-main/COSMOS2020_LePhare_v2_20210507_LIB_COMB.hdf5","r")

def dust_ext(w,law=0,ebv=0.):
    
    law1 = np.loadtxt("../SB_calzetti.dat").T
    law2 = np.loadtxt("../extlaw_0.9.dat").T
    ext_w = [law1[0],law2[0]]
    ext_k = [law1[1],law2[1]]
    if ebv>0.:
        k = np.interp(w,ext_w[law],ext_k[law])
        return np.power(10.,-0.4*ebv*k)
    else:
        return 1.
    
# select object
obj_name = sys.argv[1]

if len(sys.argv) == 2:
    # no need to crossmatch, we have the cosmos id
    # read in the catalog for only WERLS sources
    cat0 = Table.read('WERLSv1_COSMOS2020_FARMER_R1_v2.0.fits',format='fits',hdu=1)
    idx = np.where(cat0['obj']==obj_name)[0]
    targ = cat0[idx]
    zguess = float(targ['zguess'])
    print('working on object',obj_name,'with zguess = ',zguess)
else:
    # do the crossmatch
    # read in the catalog for only WERLS sources
    cat0 = Table.read(dir_in+'COSMOS2020_FARMER_R1_v2.0.fits',format='fits',hdu=1)
    # load catalog ra/dec
    catra = cat0['ALPHA_J2000']; catdec = cat0['DELTA_J2000']
    # load target ra/dec
    notes = Table.read('cooper_full_notes.csv',format='csv')
    ind = np.where(notes['obj']==obj_name)[0]
    zguess = float(notes[ind]['zguess'][0])
    print('working on object',obj_name)
    # crossmatch
    c = SkyCoord(ra=catra, dec=catdec)
    catalog = SkyCoord(ra=notes[ind]['objra'][0]*u.degree, dec=notes[ind]['objdec'][0]*u.degree)
    idx, d2d, d3d = catalog.match_to_catalog_sky(c)
    if d2d > 1*u.arcsec:
        print("!! xmatch to source is > 1 arcsec !!")
    print('distance to xmatch is ',d2d.to(u.arcsec))
    targ = cat0[idx]
    print('COSMOS ID is: ',targ['ID'])


# print best fits params of interest
print(targ['lp_zBEST'])
print('')
print(targ['lp_MFUV'])

# find UV mag
Jmag = targ['UVISTA_J_MAG'] # apparent J mag
MUV = Jmag-cosmo.distmod(zguess)+(2.5*np.log10(1+zguess)) # abs UV mag

print("check phot MUV: ",targ['UVISTA_J_MAG']-cosmo.distmod(targ['lp_zBEST'])+(2.5*np.log10(1+targ['lp_zBEST'])))
print("")
print("spec MUV: ",MUV)

# optional: change the filters you want to use
filt_use = ['CFHT_ustar', 'CFHT_u', 'HSC_g', 'HSC_r', 'HSC_i', 'HSC_z', 'HSC_y', 'UVISTA_Y', 'UVISTA_J', 'UVISTA_H', 'UVISTA_Ks', 'IRAC_CH1', 'IRAC_CH2', 'IRAC_CH3', 'IRAC_CH4']

# flux conversion from uJy to erg/cm2/s/Hz
if outflx=='cgs':
    for b in filt_use:
        targ[b+'_FLUX'] *= 1e-29
        targ[b+'_FLUX'].unit = u.erg/u.cm/u.cm/u.s/u.Hz
        targ[b+'_FLUXERR'] *= 1e-29
        targ[b+'_FLUXERR'].unit = u.erg/u.cm/u.cm/u.s/u.Hz
        
# grab just the photometry
photcat = targ[([i+'_FLUX' for i in filt_use]+[i+'_FLUXERR' for i in filt_use])]

# open zpdf catalog
pdfcat = fits.open(dir_in+'../COSMOS2020_{}_R1_v2.0_LEPHARE_PZ.fits'.format(catversion.upper()))
zpdf = pdfcat[0].data[targ['ID']][0]

wl_obs = np.array([filt_dict[i][0] for i in filt_use]).flatten()  # wavelength center of the filter used 
wl_obserr = (np.array([filt_dict[i][1] for i in filt_use])/2.).flatten()
fnu_obs = np.array([targ[i+'_FLUX'] for i in filt_use]).flatten() # Reads the measured magnitude at that wavelength
fnu_obserr = np.array([targ[i+'_FLUXERR'] for i in filt_use]).flatten() #Magnitude associated +/-error
sel = (fnu_obs>0.).flatten()


### Begin figure!
fig, ax = plt.subplots(1, figsize=[7,5])


if targ['{}_FLUX'.format(filt_use[0])].unit.to_string()=='uJy':
    ax.set_ylabel('Flux [$\mu$Jy]')
    ax.errorbar(wl_obs[sel],fnu_obs[sel],xerr=wl_obserr[sel],yerr=fnu_obserr[sel],fmt='.w', ecolor = 'w', capsize=3, elinewidth=1,zorder=2)
    ymin = min(fnu_obs[sel])*0.5
    ymax = max(fnu_obs[sel]+fnu_obserr[sel])*6
else: # assuming it's cgs
    ax.set_ylabel('Flux [$10^{-29}$ erg/cm$^2$/s/Hz]')
    ax.errorbar(wl_obs[sel],fnu_obs[sel]*1e29,xerr=wl_obserr[sel],yerr=fnu_obserr[sel]*1e29,fmt='.w', ecolor = 'w', capsize=3, elinewidth=1,zorder=2)
    ymin = min(fnu_obs[sel])*1e29*0.5
    ymax = max(fnu_obs[sel]+fnu_obserr[sel])*1e29*6
    
# Using the redshift of best-fit template
zp = targ['lp_zBEST']
m = int(targ['lp_model'])
wvl = hdf['/model{}/spectra'.format(m)].attrs['lambda[AA]'] *u.AA 
t = np.abs(hdf['/model{}'.format(m)].attrs['age']-targ['lp_age']).argmin()
flam_mod = hdf['/model{}/spectra'.format(m)][t,:] *u.erg/u.cm/u.cm/u.s/u.AA 
fnu_mod = flam_mod*(wvl**2)/const.c 
# Calculates the flux in units of [uJy] also applying dust ext
fnu_mod = fnu_mod.to(u.erg/u.cm/u.cm/u.s/u.Hz) * dust_ext(wvl.value,law=int(targ['lp_Attenuation']),ebv=targ['lp_dust'])
# Rescale the template
mscal = hdf['/model{}'.format(m)].attrs['mass'][t]/10**targ['lp_mass_best']  # luminosity/mass resc
dm = cosmo.luminosity_distance(zp)/(10*u.pc)  # distance modulus
offset = dm.decompose()**2*mscal/(1+zp) # all together * (1+z) factor


# Plot the best-fit model 
ax.plot(wvl*(1+zp),fnu_mod.to(u.uJy).value/offset,color='c',alpha=1,label='model',zorder=1) 

# Show where nebular emission lines would potentially boost the flux
ax.vlines(3727*(1+zp),ymin,ymax,label='[OII]',zorder=0,color='0.3',ls=':')
ax.vlines(5007*(1+zp),ymin,ymax,label='[OIII]b',zorder=0,color='0.3',ls=':')
ax.vlines(4861*(1+zp),ymin,ymax,label='Hb',zorder=0,color='0.3',ls=':') # H_beta
ax.vlines(6563*(1+zp),ymin,ymax,label='Ha',zorder=0,color='0.3',ls=':') # H_alpha 
ax.vlines(1216*(1+zguess),ymin,ymax,label='Lyman-alpha break',zorder=0,color='chocolate',ls='-.')

# set plot limits
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(1000,100000) 
ax.set_ylim(ymin,ymax)
ax.set_xlabel('wavelength [Å]')

# annotate
ax.annotate(text='photo-z='+str(zp[0]),xy=(0.15,0.83),xycoords='figure fraction',fontsize='small',color='w')
ax.annotate(text='spec-z='+str(zguess),xy=(0.15,0.8),xycoords='figure fraction',fontsize='small',color='w')
ax.annotate(text='MUV='+str(np.round(MUV[0],2)),xy=(0.15,0.77),xycoords='figure fraction',fontsize='small',color='w')
ax.annotate(text='J mag='+str(np.round(Jmag[0],2)),xy=(0.15,0.74),xycoords='figure fraction',fontsize='small',color='w')


# plot zpdf
# Create inset of width 30% and height 40% of the parent axes' bounding box
axin = inset_axes(ax, width="33%", height="33%", loc='lower right')
zarr = np.linspace(0,10,num=len(zpdf[1:]))
axin.plot(zarr,zpdf[1:],c='teal',lw=3,label=r'$z_{phot}=$'+str(targ['lp_zBEST'][0]))
axin.axvline(zguess,c='chocolate',lw=3,label=r'$z_{Ly-\alpha}=$'+str(zguess))
#axin.set_xlabel('Redshift',fontsize='small')
#axin.set_ylabel('LePhare PDF')
xmin = [zp-1.5 if (zp > 1.5)==True else 0][0]
axin.set_xlim(0,zp+1.5)
axin.xaxis.tick_top()
axin.tick_params(size=10)
axin.yaxis.set_ticks([])
axin.xaxis.set_label_position('top') 
#axin.set_ylim(-0.1,6.7)
axin.grid(False)

print("The COSMOS fitted model is model number",m)
print('The offset applied is',offset,'and a redshift of',zp)
#plt.show()  
plt.savefig('/Users/oc4858/Library/CloudStorage/Box-Box/werls/plots/'+obj_name+'_c2020_sed.png',dpi=500)









