import os
from pathlib import Path
from astropy.io import fits

local = (str(os.path.abspath(os.getcwd())))
path = Path(local)
print('local path:',path)
fits_files = list(path.glob('*.fits'))
print('found',len(fits_files),'.fits-files')

r_count = 0
v_count = 0

move = True # if True, files are moved according to filter

for x in fits_files:
    y = str(x)
    z = y.replace(str(path),'')   
    fitsname = z[1:]
    data = fits.open(fitsname)
    hdr = data[0].header
    filtercolor = hdr['FILTER']
    data.close()
        
    print(fitsname,filtercolor)
        
    if filtercolor == 'R':
        if move:
            os.rename(str(path)+os.sep+str(fitsname),str(path)+os.sep+'R'+os.sep+str(fitsname))
        r_count = r_count + 1
    elif filtercolor == 'V':
        if move:
            os.rename(str(path)+os.sep+str(fitsname),str(path)+os.sep+'V'+os.sep+str(fitsname))
        v_count = v_count + 1
    else:
        print('error: filter',filtercolor,'in .fits-file',fitsname)

if move:        
    print('result: moved',r_count,'to',os.sep+'R and',v_count,'to',os.sep+'V')