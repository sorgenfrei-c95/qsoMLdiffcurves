import numpy as np
import os

# read information from quasarposition.txt (copied and renamed output file from original galfit)

file = np.genfromtxt('../quasarpositions.txt',skip_header=2,usecols=(0),dtype='str')
time,xpos,xpos_err,ypos,ypos_err = np.loadtxt('../quasarpositions.txt',skiprows=2,usecols=(1,2,3,4,5),unpack=True)
print(file,time)
print(xpos,ypos)

N = len(file)

# rename all psf.fits (copied via softlinks into savelog from original galfit) to be identifyable by there date!

renaming = True # WARNING: IF TRUE, ALL PSF_FILES WILL BE RENAMED!

for i in range(N):
    date = file[i]#[14:24] (with this it would only be the date, which makes problems with duplicate nights from different telescops!) 
    print(date,time[i])
    
    copy_and_rename_psf = 'cp psf_'+str(i+1)+'.fits psf_'+str(date)+'.fits'
    print(copy_and_rename_psf)
    if renaming == True:
        os.system(copy_and_rename_psf) # psf-renaming

print('#')
print('#')
print('#')
print('#')
if len(time) != N or len(xpos) != N or len(ypos) != N:
    print('fatal error!')
else:        
    print(N,'psf-files renamed')