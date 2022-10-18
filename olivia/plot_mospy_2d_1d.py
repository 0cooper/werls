"""
Plot spectra from mospy products

added:
match pix/wavelength
mark lines
filled in sky lines 1d
blocked out sky lines 2d
gaussian smooth 2d
manual extract 1d

need to do:
mark slit pos in y axis


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
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve
from scipy.ndimage import gaussian_filter1d
import extract_1d_spec as exsp 
plt.style.use('../cooper-presentation.mplstyle')

# select object
obj_name = sys.argv[1]
notes = Table.read('../cooper_cosmos_notes.csv',format='csv')
idx = np.where(notes['obj']==obj_name)[0]
zguess = float(notes['zguess'][idx])
lguess = float(notes['wave'][idx])
ycen = int(notes['ypix_slit'][idx])
print('working on object',obj_name,'with zguess = ',zguess,'at lambda = ',lguess,'AA and mospy ypix = ',ycen)

if notes['mask'][idx] == 'wmmc01':
    ddir = '../wmmc01/wmmc01_Y_'
    fwhm = 0.456 # star_75
elif notes['mask'][idx] == 'wmmc02':
    ddir = '../wmmc02/mospy_products/COMBINED/wmmc02_COMBINED_'
    fwhm = 0.515 # star_0
elif notes['mask'][idx] == 'wmmc03':
    ddir = '../wmmc03/mospy_products/COMBINED/wmmc03b_COMBINED_'
    fwhm = 0.648 #star_144
elif notes['mask'][idx] == 'wmmc05':
    ddir = '../wmmc05/mospy_products/COMBINED/wmmc05b_COMBINED_'
    fwhm = 0.481 # star_575
elif notes['mask'][idx] == 'wmmc06':
    ddir = '../wmmc06/mospy_products/wmmc06b_Y_'
    fwhm = 0.548 # star_794
elif notes['mask'][idx] == 'wmmu01':
    ddir = '../wmmu01/mospy_products/wmmu01_COMBINED_'
    fwhm = 0.673 # star_62
elif notes['mask'][idx] == 'wmme01':
    ddir = '../wmme01/mospy_products/wmme01_COMBINED_'
    fwhm = 0.697 #star_881
elif notes['mask'][idx] == 'wmme02':
    ddir = '../wmme02/mospy_products/wmme02b_Y_'
    fwhm = 0.875 #star_513
elif notes['mask'][idx] == 'wmme03':
    ddir = '../wmme03/mospy_products/wmme03_Y_'
    fwhm = 0.558 #star_553
    

# open 1D and 2D spectrum for an object
hdu2 = fits.open(ddir+obj_name+'_eps.fits') # mospy 2D spectrum file
image = hdu2[0].data # 2D spectrum array

if len(sys.argv) == 2:
    # don't include the rms addition to spectrum
    wav, pix, spec, _, err, _ = exsp.extract1d(obj_name, ddir, ycen, FWHM=fwhm) # optimally extract 1D spectrum
    rms = np.zeros_like(wav)
else:
    # include rms addition to spectrum to fix oversubtraction
    wav, pix, spec, _, err, rms = exsp.extract1d(obj_name, ddir, ycen, FWHM=fwhm) # optimally extract 1D spectrum

pixscale = hdu2[0].header['PSCALE'] # arcsec/pix
slitwidth = hdu2[0].header['YOFFSET'] # arcsec
nodamp = slitwidth/pixscale # nod amplitude in pixels

# gaussian smooth 2D spec
kernel = Gaussian2DKernel(x_stddev=1)
conv_im = convolve(hdu2[0].data, kernel)

### START FIG!
fig, (ax1, ax2) = plt.subplots(2, sharex=True)
fig.subplots_adjust(wspace=0,hspace=0)

# 1D spec
ax1.plot(pix,gaussian_filter1d(spec+rms,sigma=1),c='chocolate',label=obj_name)
ax1.plot(pix,spec+rms,c='chocolate',ls='--',alpha=0.6) 
ax1.fill_between(pix, y1=err, y2=-1*err, color='gray', label='sky')
ax1.set_ylabel(r'Flux')

# 2D spec
lower=-1; upper=2
sample = sigma_clip(hdu2[0].data) 
vmin = sample.mean() + lower * sample.std()
vmax = sample.mean() + upper * sample.std()
ax2.imshow(conv_im, origin='lower', cmap='afmhot', aspect='auto', vmin=vmin, vmax=vmax)
ax2.set_xlabel('Column [Pix]')
ax2.set_ylabel('Row [Pix]')
ax2.grid(False)

# block out sky lines 2D
ipk = find_peaks(err, height=np.mean(err))[0]
for i in ipk:
    ax2.axvline(pix[i],c='#9E420E',lw=15,alpha=0.85)

    
# secax = secondary axis for wavelength space
def forward(x):
    return np.interp(x, pix, wav)

def inverse(x):
    return np.interp(x, wav, pix)

secax = ax1.secondary_xaxis('top', functions=(forward, inverse))
secax.xaxis.set_minor_locator(AutoMinorLocator())
secax.set_xlabel(r'Wavelength [$\AA$]')

# mark measured line at pixel guess based on lambda guess
pguess = pix[rmp.closest(wav,lguess)[0]] # pixel closest to wavelength of line to center on
print("candidate line at ",lguess,'AA; ',pguess,'pix')
s = (spec+rms)[rmp.closest(wav,lguess)[0]] # signal closest to wavelength of line
n = err[rmp.closest(wav,lguess)[0]] # error closest to wavelength of line
sn = np.round(s/n,2)
print("detected at ",sn,' sigma')
ax1.axvline(pguess,lw=1,c='k',ls='-.',alpha=0.4)
ax2.axvline(pguess,lw=1,c='k',ls='-.',alpha=0.4)

# mark obj position on slit and nods
ax2.axhline(ycen,lw=1,c='k',ls='-.',alpha=0.4)


# limits
ymin = [ycen-25 if (ycen > 25)==True else 0][0]
ymax = [ycen+25 if (len(conv_im) > ycen+25)==True else len(conv_im)-1][0]
ax1.set_xlim(pguess-70,pguess+70)
ax1.set_ylim(-2*np.nanmedian(err),4*np.nanmedian(err))
ax2.set_xlim(pguess-70,pguess+70)
ax2.set_ylim(ymin,ymax) 

#plt.show()
# save fig
plt.savefig('/Users/oc4858/Library/CloudStorage/Box-Box/werls/plots/'+obj_name+'_mospy_1D_2D.png',dpi=500)
