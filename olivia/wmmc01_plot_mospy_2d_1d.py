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

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import astropy
from astropy import units as u
from astropy.io import fits
from astropy.table import Table
from astropy.visualization import ZScaleInterval
from astropy.stats import sigma_clip
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
import read_mospy_files as rmp ## script I made to read in mospy files
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve
plt.style.use('../../cooper-presentation.mplstyle')

# select object
notes = Table.read('../wmmc01/cooper_notes_wmmc01.csv',format='csv')
idx = 19
obj_name = notes['obj'][idx]; zguess = float(notes['zmeasured'][idx]); lguess = float(notes['wave'][idx])
print('working on object',obj_name,'with zguess = ',zguess,'at lambda = ',lguess,'AA')


# open 1D and 2D spectrum for an object
hdu2 = fits.open('../wmmc01/wmmc01_Y_'+obj_name+'_eps.fits') # 2D spectrum
hdu1 = fits.open('../wmmc01/wmmc01_Y_'+obj_name+'_1D_00.fits') # 1D spectrum
pix_obj, wav_obj, spec1d_obj = rmp.make_1d_array(ext=0,hdu=hdu1) # object spectrum
pix_sky, wav_sky, spec1d_sky = rmp.make_1d_array(ext=1,hdu=hdu1) # sky spectrum
image = hdu2[0].data
#obj_pos = int(notes['spat_fracpos'][idx]*image.shape[0])-18
pixscale = hdu2[0].header['PSCALE'] # arcsec/pix
slitwidth = 2.5 # arcsec
nodamp = slitwidth/pixscale

# manually extract 1D spec
ypix = 19 ### change me!
obj_pos = 19 ### change me!
row1, row2 = ypix-4, ypix+4 # define the target aperture range
spec = np.sum(image[row1:row2, :], axis=0)

# gaussian smooth 2D spec
kernel = Gaussian2DKernel(x_stddev=1)
conv_im = convolve(hdu2[0].data, kernel)

### START FIG!
fig, (ax1, ax2) = plt.subplots(2, sharex=True)
fig.subplots_adjust(wspace=0,hspace=0)

# 1D spec
ax1.plot(pix_obj,spec,c='chocolate',label=obj_name)
ax1.fill_between(pix_sky, y1=spec1d_sky, y2=0, color='gray', label='sky')
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
ipk = find_peaks(spec1d_sky, height=np.mean(spec1d_sky))[0]
for i in ipk:
    ax2.axvline(pix_sky[i],c='#9E420E',lw=15,alpha=0.85)

    
# secax = secondary axis for wavelength space
def forward(x):
    return np.interp(x, pix_obj, wav_obj)

def inverse(x):
    return np.interp(x, wav_obj, pix_obj)

secax = ax1.secondary_xaxis('top', functions=(forward, inverse))
secax.xaxis.set_minor_locator(AutoMinorLocator())
secax.set_xlabel(r'Wavelength [$\AA$]')

# mark measured line at pixel guess based on lambda guess
pguess = pix_obj[rmp.closest(wav_obj,lguess)[0]] # pixel closest to wavelength of line to center on
print("candidate line at ",lguess,'AA; ',pguess,'pix')

# mark obj position on slit and nods
ax2.axhline(obj_pos,xmin=0,xmax=0.1,lw=2,c='chocolate')
ax2.axhline(obj_pos,xmin=0.9,xmax=1,lw=2,c='chocolate')
ax2.axhline(obj_pos+nodamp,xmin=0,xmax=0.1,lw=2,c='k')
ax2.axhline(obj_pos+nodamp,xmin=0.9,xmax=1,lw=2,c='k')
ax2.axhline(obj_pos-nodamp,xmin=0,xmax=0.1,lw=2,c='k')
ax2.axhline(obj_pos-nodamp,xmin=0.9,xmax=1,lw=2,c='k')


# limits
ax1.set_xlim(pguess-70,pguess+70)
ax1.set_ylim(-0.04,0.08)
ax2.set_xlim(pguess-70,pguess+70)
ax2.set_ylim(0,50) 
ax1.annotate('z = '+str(zguess),xy=(pguess-45,0.065),fontsize=10,c='chocolate')

fig.savefig('plots/wmmc01/'+obj_name+'_manual_boxcar_1D_2D_lines.png',dpi=500)
