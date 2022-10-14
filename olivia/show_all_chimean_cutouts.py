#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 11:40:41 2022

@author: oc4858
"""

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
keep = Table.read('/Users/oc4858/werls/cosmos_keep.csv',format='csv')
objs = keep['obj']
notes = Table.read('cooper_full_notes.csv',format='csv')

for obj_name in objs:
    
    idx = np.where(notes['obj']==obj_name)[0]
    targ = notes[idx]
    
    if notes['mask'][idx] == 'wmmc01':
        file = '../wmmc01/chimean_wmmc01.fits' 
        source = 'chimean'
        reg = '../wmmc01/wmmc01_SlitRegions.reg'
    elif notes['mask'][idx] == 'wmmc02':
        file = '../wmmc02/mask_design/chimean_wmmc02.fits'
        source = 'chimean'
        reg = '../wmmc02/mask_design/wmmc02_SlitRegions.reg'
    elif notes['mask'][idx] == 'wmmc03':
        file = '../wmmc03/mask_design/chimean_wmmc03.fits'
        source = 'chimean'
        reg = '../wmmc03/mask_design/wmmc03_SlitRegions.reg'
    elif notes['mask'][idx] == 'wmmc05':
        file = '../wmmc05/mask_design/chimean_wmmc05.fits'
        source = 'chimean'
        reg = '../wmmc05/mask_design/wmmc05_SlitRegions.reg'
    elif notes['mask'][idx] == 'wmmc06':
        file = '../wmmc06/mask_design/chimean_wmmc06.fits'
        source = 'chimean'
        reg = '../wmmc06/mask_design/wmmc06_SlitRegions.reg'
    
    
    # target
    position = SkyCoord(targ['objra'][0],targ['objdec'][0],unit=(u.deg,u.deg))
    
    # the aplpy way :) 
    gc = aplpy.FITSFigure(file,north=False)
    gc.show_grayscale(pmax=85)
    gc.recenter(position.ra, position.dec, radius=5/3600)
    gc.show_regions(reg)
    gc.add_grid()
    gc.remove_grid()
    gc.set_theme('publication')
    
    plt.savefig('/Users/oc4858/Library/CloudStorage/Box-Box/werls/plots/'+obj_name+'_'+source+'_cutout.png',dpi=500)





