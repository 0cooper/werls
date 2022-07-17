"""
Optimal and boxcar 1D spectral extraction for MOSFIRE 2D spectra.

Optimal extraction adapted from Taylor Hutchison's code hosted here: https://github.com/aibhleog/simply-spectra/blob/master/quick_vis.py

@author: 0cooper

"""
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import os

# for the optimized extraction shape

def gaussian(xaxis, mean, A, sig, offset): 
	'''
	Simple Gaussian function, to be used in the quick optimized extraction
	'''
	return A * np.exp(-np.power(xaxis-mean, 2.) / (2*np.power(sig, 2.))) + offset

def extract1d(obj,path,ycen,aper=7,width=4):
	'''
	Takes a mospy reduced MOS 2D image for a single object and manually extracts the 1D spectrum at 	a given y pixel using both boxcar and optimal extraction methods.
	
	INPUTS ---- obj:		str, name of the reduced 2D MOSFIRE file to be read in
		        path:		str, points to data directory
		        ycen:		int, row to extract spectrum on (pix)
				aper:		int, number of rows to extract the 1D spectrum over
				width:		int, pixel width for boxcar box
				
	RETURNS --- wave:		wavelength array
		        spec:		optimal extracted spectral array
                err:		error for optimal extracted spectral array
				specbox:	boxcar extracted spectral array
                errbox:		error for boxcar extracted spectral array
	'''

	# making sure the aperture is an odd number
	assert aper%2 == 1, "Aperture size needs to be an odd number of pixels. "\
		f"Currently, the aperture size is: \n \t\t  aper = {aper} pixels (default is 7 pixels)."

	# reading in data
	print(f'\nReading in data for {obj}')
	header = fits.getheader(path + obj + '_eps.fits')
	signal = fits.getdata(path + obj + '_eps.fits')
	error = fits.getdata(path + obj + '_sig.fits')
	print(f'Dimensions: \t signal spectrum {signal.shape}\n' +
			f'\t\t error spectrum {error.shape}')
	
	wavelength_start = header['CRVAL1'] # starting wavelength at first pixel
	wavelength_logdisp = header['CD1_1'] # delta wavelength per pixel
	num_wavelength = header['NAXIS1'] # length of data array
	wave = wavelength_start + np.arange(0, wavelength_logdisp*num_wavelength, wavelength_logdisp) # wavelength

	# defining optimized extraction gaussian
	pixscale = header['PSCALE'] # arcsec/pix
	fwhm = 0.8 / pixscale # arcsec / [arcsec/pixel]
	gauss = gaussian(np.arange(aper),mean=3.,A=1.,sig=fwhm/2.35,offset=0.)
	gauss /= sum(gauss) # to make it sum to 1 to use as weights
	gauss_2D = np.zeros((len(gauss),len(wave))) # making 2D array of weights
	for i in range(aper):
		gauss_2D[i] = gauss[i]

	# optimally-extracting 1D spectra
	half = int(aper/2) # to make the cut out of the 2D image
	spec = np.nansum(signal[ycen-half:ycen+half+1].copy()*gauss_2D,axis=0)
	err = np.nansum(error[ycen-half:ycen+half+1].copy()*gauss_2D,axis=0)
    
	# boxcar extract 1D spec
	row1, row2 = ycen-width, ycen+width # define the target aperture range
	specbox = np.sum(signal[row1:row2, :], axis=0)
	errbox = np.sum(error[row1:row2, :], axis=0)
    
	return wave, spec, err, specbox, errbox