"""
Grab spectra from mospy products

@author: 0cooper

"""

# the basics

import numpy as np
import astropy
from astropy.io import fits
from astropy.table import Table

# translate the 1D header data into a wavelength array

def make_1d_array(ext,hdu):
    """
    object ext = 0
    sky ext = 1
    """
    im1d = hdu[ext] # point to the extension we want
    wavelength_start = im1d.header['CRVAL1'] # starting wavelength at first pixel
    wavelength_logdisp = im1d.header['CD1_1'] # delta wavelength per pixel
    num_wavelength = im1d.header['NAXIS1'] # length of data array

    pix_start = im1d.header['CRPIX1']
    pix_logdisp = im1d.header['CD2_2']
    num_pix = im1d.header['NAXIS1']

    pix = pix_start + np.arange(0, pix_logdisp*num_pix, pix_logdisp)
    wavelength = wavelength_start + np.arange(0, wavelength_logdisp*num_wavelength, wavelength_logdisp) # wavelength
    spec1d = im1d.data # spectrum
    
    return pix,wavelength,spec1d

# function to find closest value in array to specific number

def closest(lst, val):
    lst = np.asarray(lst) 
    idx = (np.abs(lst - val)).argmin() 
    return idx,lst[idx]