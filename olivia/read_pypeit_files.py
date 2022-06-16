"""
Grab spectra from pypeit products

@author: 0cooper

"""

# the basics

import numpy as np
import astropy
from astropy.io import fits
from astropy.table import Table
plt.style.use('../cooper-paper.mplstyle')


def pypeit_1d_to_txt(mask,file):
    """
    ex:
    mask = 'wmmc01'
    file = 'spec1d_m220212-m220212-wmmc01.fits'
    """

    pyp1 = fits.open('../'+mask+'/'+file)
    nspec = pyp1[0].header['NSPEC']

    for n in range(1,nspec+1):
        obj = pyp1[n].header['HIERARCH MASKDEF_OBJNAME']

        if obj == 'SERENDIP':
            slit = pyp1[n].header['SLITID']
            obj = 'SERENDIP'+'_slit_'+str(slit)
        else:
            pass

        print("pulling 1D spec for obj",obj)
        wav = []
        cts = []
        sig = []
        sky = []
        bwav = []
        bcts = []
        bsig = []
        bsky = []
        for j in range(len(pyp1[n].data)):
            wav.append(pyp1[n].data[j][2])
            cts.append(pyp1[n].data[j][3])
            sig.append(pyp1[n].data[j][5])
            sky.append(pyp1[n].data[j][7])
            bwav.append(pyp1[n].data[j][11])
            bcts.append(pyp1[n].data[j][12])
            bsig.append(pyp1[n].data[j][14])
            bsky.append(pyp1[n].data[j][16])
        tab = Table([wav,cts,sig,sky,bwav,bcts,bsig,bsky],\
            names=('lambda','opt_counts','opt_sigma','sky_counts','box_lambda','box_counts','box_sigma','box_sky_counts'))
        tab = tab[tab['lambda']>1]
        fname = obj+'_1dspec.txt'
        tab.write('../'+mask+'/pypeit_spec1d/'+fname,format='ascii',overwrite=True)
        
    return