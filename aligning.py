import numpy as np
from astropy.io import fits
import os
from pathlib import Path
import sys

# before running this code: mkdir of interp_trimmed and aligned!
# and select good reference file and name it astrometry_ref.fits in the filter-directory!

print('')
print('#######')
print('aligning.py:')
print('#######')
print('')

aligning_bool = True # if True, images will be aligned with registration code
trimming_bool = True # if True, images will be trimmed once more after aligning
# If both False, esentially only log_interp2 will be read and a list with problematic files generated!

# possibilty to correct the daily refstar positions of isis with respect to good.data with gaia propermotions
gaiastarpositioncorrection_bool = True # if True, gaiastarpositioncorrection_bool.txt will say True,
                                       # which then allows the correction program (see next line) to do its job. 
                                       # Check gaiaposcorrectionforaligning.py for more details.
                                       
txt = open("gaiastarpositioncorrection_bool.txt", "w")
if gaiastarpositioncorrection_bool == True:
    txt.write('True')
    print('star positions will be corrected with gaia data')
else:
    txt.write('False')
    print('star positions will NOT be corrected with gaia data')
txt.close()

print('')
print('#######')
print('Step 1: preparation of registration program to aligne all images')
print('#######')
print('')

local = (str(os.path.abspath(os.getcwd())))
local_path = Path(local)
trimmed_path = Path(local+os.sep+'trimmed')
registered_path = Path(local+os.sep+'interp_trimmed')
print('local path:       ',local_path)
print('trimmed data path:',trimmed_path)
print('registered path:  ',registered_path)

fits_files = list(trimmed_path.glob('*.fits'))
print('number of total fits-files:',len(fits_files))

if aligning_bool == True:

    # verify trimmed files shape or interrupt script:
    trimmed_files_shape = 3900
    for x in fits_files:
        y = str(x)
        z = y.replace(str(local_path),'')   
        fitsname = z[1:]
        data = fits.open(fitsname)
        shape = data[0].data.shape
        data.close()
        if shape != (trimmed_files_shape,trimmed_files_shape):
            print('warning: wrong shape:',shape,'in:',fitsname)
            sys.exit('FileShapeError: interrupted aligning.py!') # stop program directly if trimmed shapes are wrong!
    
    # make process_config file needed for registration:    
    process_config = open("process_config", "w")
    IM_DIR='../../'+str(str(local_path).split('/lco/')[1])
    #print(IM_DIR)
    process_config.write(
    "IM_DIR          "+IM_DIR+"                           Directory where the images are \n"+
    "MRJ_DIR         ../../isis                           Installation directory \n"+
    "REFERENCE       astrometry_ref.fits                  Reference image for astrometry \n"+
    "REF_SUB         astrometry_ref.fits                  Reference image for subtraction \n"+
    "INFILE          ../../isis/register/dates            Dates of the frames \n"+
    "VARIABLES       phot.data                            coordinates of objects for which we want to make light curves \n"+
    "DEGREE          2                                    The degree of the polynomial astrometric transform between frames \n"+
    "CONFIG_DIR      ../../isis/register                  Where to find the configuration files \n"+
    "SIG_THRESH      2.0 \n"+
    "COSMIC_THRESH   60000.0 \n"+
    "REF_STACK       astrometry_ref.fits \n"+
    "N_REJECT        2 \n"+
    "MESH_SMOOTH     1 \n")
    process_config.close()
    os.system('mv process_config ../../isis/register')
    
    # make dates file needed for registration:
    os.system('/bin/ls trimmed/*.fits > ../../isis/register/dates')
    
    # check whether a (nice, prefarable low seeing and good) reference image
    # was selected and copied from the trimmed folder to the filter folder 
    # (i.e. the local directory) and renamed to astrometry_ref.fits:
    if os.path.isfile('astrometry_ref.fits') != True:
        sys.exit('RefFileError: astrometry_ref.fits is missing!') # stop program directly if there is no ref-image!
    else:
        ref_data = fits.open('astrometry_ref.fits')
        ref_shape = ref_data[0].data.shape
        ref_data.close()
        if ref_shape != (trimmed_files_shape,trimmed_files_shape):
            print('warning: wrong shape:',ref_shape,'in astrometry_ref.fits')
            sys.exit('FileShapeError: interrupted aligning.py!') # stop program directly if the shape of the ref-image is wrong!
 
print('')
print('#######')
print('Step 2: actual alignement of the images')
print('#######')
print('')

if aligning_bool == True:
    
    # change directory to execute interp.csh i.e. start isis-registration:
    print(os.getcwd())
    os.chdir('../../isis/register/')
    print(os.getcwd())
    
    # execute interp.csh and therefore start aligning with isis-registration-program:
    print('executing registration')
    os.system('./interp.csh > alignpy_consollog.txt')
    os.system('mv alignpy_consollog.txt '+str(local_path)) # move hidden consol output as txt-file to local

    # back to the original local filter directory:
    print(os.getcwd())
    os.chdir(local_path)
    print(os.getcwd())

####################################################################################################################
################## NOW: SECOND TRIMMING OF THE ALIGNED IMAGES VON INTERP_TRIMMED NACH ALIGNED ######################
####################################################################################################################

# ursprünglich "trimming2_ofinterp" --> muss nicht mehr zusätzlich ausgeführt werden!

print('')
print('#######')
print('Step 3: trimming of aligned images')
print('#######')
print('')

if trimming_bool == True:
    
    # wirft die ersten und letzten "trimming" Pixel in Zeilen und Spalten weg
    trimming = 100               
    print('trimming =',trimming)
    # wenn t_*fits 3900*3900 pixel shape haben (wie bei he2149 und q2237) muss jetzt 3700*3700 rauskommen!
    
    data_path_ALIGNED = Path(local+os.sep+'interp_trimmed')
    output_path_ALIGNED = Path(local+os.sep+'aligned')
    print('local path:         ',local_path)
    print('registerd data path:',data_path_ALIGNED)
    print('final output path:  ',output_path_ALIGNED)
    
    fits_files_ALIGNED = list(data_path_ALIGNED.glob('*.fits'))
    print('number of total fits-files:',len(fits_files_ALIGNED))
    
    count = 1
    for x in fits_files_ALIGNED:
        y = str(x)
        z = y.replace(str(local_path),'')   
        fitsname = z[1:]
        data = fits.open(fitsname)
        shape = data[0].data.shape
        
        if shape != (trimmed_files_shape,trimmed_files_shape):
            print('warning: unusual shape:',data[0].data.shape,'in:',fitsname)
            break
        
        array = data[0].data
        trimmed_array = array[trimming:-trimming,trimming:-trimming]
        data[0].data = trimmed_array
        new_shape = data[0].data.shape
        
        if new_shape != (trimmed_files_shape-2*trimming,trimmed_files_shape-2*trimming):
            print('warning: unusual new shape:',data[0].data.shape,'in:',fitsname)
            break # cannot happen!
        
        data.writeto(str(output_path_ALIGNED)+os.sep+'a_'+str(fitsname[15:]),overwrite=True)
        print('trimmed output',count,'of',len(fits_files_ALIGNED),': ...'+os.sep+'aligned'+os.sep+'a_'+str(fitsname[15:]))
        count = count + 1
        data.close()
    
    # Since the aligned files from interp_trimmed are now once more trimmed and saved in aligned,
    # the files in interp_trimmed are removed now (only if and because trimming_bool = True here) to save space:
    print('removing interp_trimmed-files after all images are trimmed and saved in the aligned-folder!')
    os.system('rm interp_trimmed/t_*.fits')
    
######################## ADDITIONAL STEP ##########################
    
print('')
print('#######')
print('Step 4: identify bad images with bad sigmas from log_interp list')
print('#######')
print('')

error_file_count = 0
error_file_list = []
log = open("log_interp2","r")
for line in log:
    Line = line.split()
    if Line[0] != "sigmax:":
        print('WARNING:',line)
        error_file_count = error_file_count + 1
        try:
            filenamestart = line.find("trimmed/")
            filenameend = line.find(".fits")+5
            filename = line[filenamestart:filenameend]
            #print(filename)
            error_file_list.append(filename)
        except:
            print('WARNING: no filename was extracted... Look at log_interp2')        
    else:
        #print(Line)
        log_vals = np.array([float(Line[1]),float(Line[3]),int(Line[5]),int(Line[7])])
        #print(log_vals)
        bad_sigma_1 = log_vals[0] > 0.4 or log_vals[1] > 0.4
        bad_sigma_2 = log_vals[0] > 0.3 and log_vals[1] > 0.3
        bad_sigma_3 = log_vals[0] > 0.2 or log_vals[1] > 0.2
        bad_ndata = log_vals[2] < 100 
        bad_nrest = log_vals[3] < 100
        if bad_sigma_1 or bad_sigma_2 or (bad_sigma_3 and (bad_ndata or bad_nrest)) or (bad_ndata and bad_nrest):
            print('WARNING:',log_vals,'in',Line[8])
            error_file_count = error_file_count + 1
            error_file_list.append(Line[8])
log.close()
    
print('')
print('WARNING from log_interp2:',error_file_count,'potentially bad files')
print(error_file_list)
print('')
print('string of',error_file_count,'bad images in aligned:')
bad_string = '-'.join(error_file_list).replace('-trimmed/',' a_').replace('trimmed/','a_')
print(bad_string)
print('')
print('aligning.py is finished!')
