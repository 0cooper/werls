"""
Make regions file from table

@author: 0cooper

"""


# the basics
import astropy
from astropy.table import Table
from astropy.io import fits



def table_to_reg(filepath,outname,ra,dec,aprad=20,color='green',fmt='ascii'):
    
    """
    Takes csv or fits file table of objects and converts it into a regions file
    
    INPUT:       
            filepath = filepath for sample to make regions for
            
            outname = file name for output regions file
            
            ra = column name in csv file for RA
                
            dec = column name in csv file for Dec
           
            aprad = aperture radius for circles around targets in pixels
                Default: 20 pixels
            
            color = color of regions (green, red, blue, magenta, cyan, white)
                Default: green
                
            fmt = format of input file
                Default: 'ascii'
          
    OUTPUT:
            returns ds9 regions file
    
    
    """
    # make ds9 regions file for visual checking of object list
    ot = Table.read(filepath,format=fmt)

    regfile=str(outname)+'.reg'
    r = open(regfile,'w')
    x = ot[ra]
    y = ot[dec]
    for j in range(0,len(x)):
        r.write('circle '+str(x[j])+' '+str(y[j])+' '+str(aprad)+' '+'# text={'+str(j)+'} color='+str(color)+'\n')  # draw circles for apertures
    r.close()
    print('     Wrote ds9 regions file ', regfile)
    print('     To use in ds9: Regions > Load Regions > <filename> ')
    print('')
    