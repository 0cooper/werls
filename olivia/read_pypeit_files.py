"""
Grab spectra from pypeit products

@author: 0cooper

"""

# the basics

import numpy as np
import astropy
from astropy.io import fits
from astropy.table import Table
import glob

def pypeit_1d_to_txt(mask):
    """
    ex:
    mask = 'wmmc01'
    file = 'spec1d_m220212-m220212-wmmc01.fits'
    """

    file = glob.glob('../'+mask+'/pypeit_products/Science_coadd_1x1/spec1d*.fits')[0]
    #file = glob.glob('../'+mask+'/spec1d*.fits')[0]
    pyp1 = fits.open(file)
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
        pix = []
        cts = []
        sig = []
        sky = []
        bwav = []
        bcts = []
        bsig = []
        bsky = []
        c=0
        for j in range(len(pyp1[n].data)):
            wav.append(pyp1[n].data[j][2])
            pix.append(c)
            cts.append(pyp1[n].data[j][3])
            sig.append(pyp1[n].data[j][5])
            sky.append(pyp1[n].data[j][7])
            bwav.append(pyp1[n].data[j][11])
            bcts.append(pyp1[n].data[j][12])
            bsig.append(pyp1[n].data[j][14])
            bsky.append(pyp1[n].data[j][16])
            c +=1
        tab = Table([wav,pix,cts,sig,sky,bwav,bcts,bsig,bsky],\
            names=('lambda','pix','opt_counts','opt_sigma','sky_counts','box_lambda','box_counts','box_sigma','box_sky_counts'))
        tab = tab[tab['lambda']>1]
        fname = obj+'_1dspec.txt'
        tab.write('../'+mask+'/pypeit_products/Science_coadd_1x1/'+fname,format='ascii',overwrite=True)
        #tab.write('../'+mask+'/'+fname,format='ascii',overwrite=True)
        
    return