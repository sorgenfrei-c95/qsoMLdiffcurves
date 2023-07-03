import numpy as np
from astropy.io import fits
import os
from pathlib import Path

targetsize = 3900 # target image shape: targetsize*targetsize pixel

local = (str(os.path.abspath(os.getcwd())))
path = Path(local)
raw = Path(str(local+'/raw'))
fits_files = list(raw.glob('*.fits'))
print('number of total fits-files:',len(fits_files))

count = 1
shape_count = 0
fail_count = 0
for x in fits_files:
    y = str(x)
    z = y.replace(str(path),'')   
    rawfitsname = z[1:]
    fitsname = rawfitsname[4:]
    #print(rawfitsname,fitsname)
    data = fits.open(rawfitsname)
    shape = data[0].data.shape
    if shape != (4096,4096):
        shape_count = shape_count + 1
        print('WARNING: unusual shape:',shape,'in',fitsname)
        
    if shape[0] <= targetsize or shape[1] <= targetsize:
        fail_count = fail_count + 1
        print('WARNING: size smaller than targetsize!',shape,'in',fitsname)
        print('--> skip image!')
        count = count + 1 # skip image and go to next image
        continue
    
    xtrim = int((shape[0]-targetsize)/2.0)
    ytrim = int((shape[1]-targetsize)/2.0)
    #print(xtrim,ytrim)
  
    # rotate elp-data:
    if fitsname[0:3] == 'elp':
        array = np.flip(data[0].data,(0,1))
        print('elp-rot.: '+fitsname)
        
    # flip early cpt-data:
    elif fitsname[0:8] == 'cpt1m013' and int(fitsname[14:22]) < 20190501: #20190701: happend two months earlier!?
        array = np.flip(data[0].data,0)
        print('cpt-flip: '+fitsname)
    
    # correct for missing ds9-y-line 1986 i.e. python-x-line 1985 in old e90-images for better aligning later:
    elif fitsname[-8:-5] == 'e90' and int(fitsname[14:22]) < 20160501:
        #print(fitsname)
        array = data[0].data.copy()                                                 # so the array has right shape and the first 1985 lines will be kept unchanged
        array[1985,:] = np.mean([data[0].data[1984,:],data[0].data[1985,:]],axis=0) # new line 1986 made with mean of neighboring lines
        array[1986:,:] = data[0].data[1985:-1,:]                                    # rest unchanged, but without the last line to keep the shape unchanged
        print('e90-line: '+fitsname)                                   
    
    # in any case save data in array:
    else:
        array = data[0].data
    
    # trimming (only if wanted, can be switched of with the following switch):
    trimming = True
    if trimming == True:  
        # eigentliches trimming (wobei sehr selten mögliche ungerade sizes beachtet werden müssen):
        if shape[0] % 2 == 0 and shape[1] % 2 == 0:
            trimmed_array = array[xtrim:-xtrim,ytrim:-ytrim]
        elif shape[0] % 2 == 0 and shape[1] % 2 != 0:
            trimmed_array = array[xtrim:-xtrim,ytrim+1:-ytrim]
        elif shape[0] % 2 != 0 and shape[1] % 2 == 0:
            trimmed_array = array[xtrim+1:-xtrim,ytrim:-ytrim]
        else:
            trimmed_array = array[xtrim+1:-xtrim,ytrim+1:-ytrim]
            
        # check size:
        #print(trimmed_array.shape)
        if trimmed_array.shape != (targetsize,targetsize):
            # cannot happen!!!
            print('WARNING: IMPOSSIBLE POST-TRIMMING-SHAPE-ERROR:',trimmed_array.shape,'in',fitsname)
            break
    else:
        print('!!!!! WARNING: TRIMMING WAS SWITCHED OFF !!!!!')
        trimmed_array = array
    
    data[0].data = trimmed_array
    
    data.writeto(str(path)+os.sep+'trimmed'+os.sep+'t_'+str(fitsname),overwrite=True)
    print('trimmed output',count,'of',len(fits_files),': ...'+os.sep+'trimmed'+os.sep+'t_'+str(fitsname))
    count = count + 1
    data.close()
    
print('unusual shapes in',shape_count,'cases.')
print('too small image in',fail_count,'cases.')
print('trimming.py is finished!')
