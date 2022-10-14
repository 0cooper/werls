#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 16:03:47 2022

@author: oc4858
"""

# the basics
import sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy import units as u 
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.visualization import ZScaleInterval
import aplpy
plt.style.use('../cooper-presentation.mplstyle')


# select object
obj_name = sys.argv[1]
notes = Table.read('cooper_full_notes.csv',format='csv')
idx = np.where(notes['obj']==obj_name)[0]
targ = notes[idx]

if notes['mask'][idx] == 'wmmc01':
    file = '../wmmc01/chimean_wmmc01.fits' 
    #file = '/Users/oc4858/Library/CloudStorage/Box-Box/werls/wmmc01/0001_150.14783000_2.19219000_hlsp_candels_hst_wfc3_cos-tot_f160w_v1.0_drz.fits'
    #source = 'candels_f160w'
    source = 'chimean'
    reg = '../wmmc01/wmmc01_SlitRegions.reg'
elif notes['mask'][idx] == 'wmmc02':
    file = '/Users/oc4858/Library/CloudStorage/Box-Box/werls/wmmc02/0001_150.10250000_2.23876000_hlsp_candels_hst_wfc3_cos-tot_f160w_v1.0_drz.fits'
    source = 'candels_f160w'
    reg = '../wmmc02/mask_design/wmmc02_SlitRegions.reg'
elif notes['mask'][idx] == 'wmmc03':
    #file = '../wmmc03/mask_design/chimean_wmmc03.fits'
    file = '/Users/oc4858/Library/CloudStorage/Box-Box/werls/wmmc03/0001_150.09307000_2.36599000_hlsp_candels_hst_wfc3_cos-tot_f160w_v1.0_drz.fits'
    source = 'candels_f160w'
    reg = '../wmmc03/mask_design/wmmc03_SlitRegions.reg'
elif notes['mask'][idx] == 'wmmc05':
    file = '../wmmc05/mask_design/chimean_wmmc05.fits'
    #file = '/Users/oc4858/Library/CloudStorage/Box-Box/werls/wmmc05/0001_150.17318000_2.40106000_hlsp_candels_hst_wfc3_cos-tot_f160w_v1.0_drz.fits'
    source = 'candels_f160w'
    reg = '../wmmc05/mask_design/wmmc05_SlitRegions.reg'
elif notes['mask'][idx] == 'wmmc06':
    file = '../wmmc06/mask_design/chimean_wmmc06.fits'
    #file = '/Users/oc4858/Library/CloudStorage/Box-Box/werls/wmmc06/0001_150.06722000_2.48847000_hlsp_candels_hst_wfc3_cos-tot_f160w_v1.0_drz.fits'
    #source = 'candels_f160w'
    source = 'chimean'
    reg = '../wmmc06/mask_design/wmmc06_SlitRegions.reg'
elif notes['mask'][idx] == 'wmmu01':
    file = '/Users/oc4858/Library/CloudStorage/Box-Box/treasurechest/uds/hlsp_candels_hst_wfc3_uds-tot_f160w_v1.0_drz.fits'  
    source = 'candels_f160w'
    reg = '../wmmu01/mask_design/wmmu01_SlitRegions.reg'
elif notes['mask'][idx] == 'wmme01':
    #file = '/Users/oc4858/Library/CloudStorage/Box-Box/treasurechest/egs/egs_all_wfc3_ir_f160w_030mas_v1.9_drz.fits'
    file = '../egs_all_wfc3_ir_f160w_030mas_v1.9_drz.fits'
    source = 'candels_f160w'
    reg = '../wmme01/mask_design/wmme01_SlitRegions.reg'
elif notes['mask'][idx] == 'wmme02':
    file = '../egs_all_wfc3_ir_f160w_030mas_v1.9_drz.fits'
    source = 'candels_f160w'
    reg = '../wmme02/mask_design/wmme02_SlitRegions.reg'
elif notes['mask'][idx] == 'wmme03':
    file = '../egs_all_wfc3_ir_f160w_030mas_v1.9_drz.fits'
    source = 'candels_f160w'
    reg = '../wmme03/mask_design/wmme03_SlitRegions.reg'


# target
position = SkyCoord(targ['objra'][0],targ['objdec'][0],unit=(u.deg,u.deg))

# the aplpy way :) 
gc = aplpy.FITSFigure(file,north=False)
gc.show_grayscale(pmax=85)
#gc.show_grayscale(pmin=0.05,pmax=95,stretch='power',exponent=3) #pmin=0.25,pmax=90,
#gc.show_markers(position.ra, position.dec, edgecolor='green', facecolor='g',
     #marker='*', s=100 )
gc.recenter(position.ra, position.dec, radius=5/3600)
gc.show_regions(reg)
gc.add_grid()
gc.remove_grid()
gc.set_theme('publication')

plt.savefig('/Users/oc4858/Library/CloudStorage/Box-Box/werls/plots/'+obj_name+'_'+source+'_cutout.png',dpi=500)


### the old 2D cutout way

#wcs = WCS(image[0].header)

# make cutout centered on object
#cutout = Cutout2D(image[0].data, position, size=15*u.arcsec, wcs=wcs)
#ax = plt.subplot(projection=wcs)
#ax.grid(False)
#interval=ZScaleInterval()
#vmin,vmax = interval.get_limits(cutout.data)
#ax.imshow(cutout.data, origin='lower', vmin=vmin,vmax=vmax)
#ax.scatter(150.19379, 2.17021, transform=ax.get_transform('fk5'), s=51, 
 #          edgecolor='m', facecolor='m',lw=5,zorder=30)
#plt.xlabel('RA')
#plt.ylabel('Dec')
#plt.show()





