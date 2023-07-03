import numpy as np
import os
import time
from pathlib import Path
from astropy.io import fits

############################################################################################################################
#### IMPORTANT!!! ##########################################################################################################
# First use for a new quasar "True" and follow the steps, then only use "False", which is the actual difference imageing!
make_stamps = False
# if True, stamps_*.reg will be made using the refimage. 
#               --> Then load all stamps_*.reg onto the refernce image in ds9.
#               --> Exclude by visual inspection bad stamps and thus selecte the stamps, needed and used in the main step.
#               --> Save all the selected stamps in ds9 as stamps_lco.reg in the band-folder (not in the diff-folder).
#               --> Run the python-script stampregionfiletotxtconvertion.py in the band folder, which makes stamps_lco.txt.
# if False, the actual diffimaging with hotpants of all images (using stamps_lco.txt) will be done, i.e. the main step.
############################################################################################################################
############################################################################################################################

parallel = True # if True, 16 images will be diffed in parallel and simultaneously

####################################################################################
ref_file_time = '2018-22' #### ADAPT THIS TIME!!! ####
print('difference imaging with reference image time:',ref_file_time)
if ref_file_time != '2014-22':
    diffimages_folder = 'diffimages_'+ref_file_time
else:
    diffimages_folder = 'diffimages'
####################################################################################

def min_adder(dire,file):
    data = fits.open(str(dire)+os.sep+str(file))
    mini = np.min(data[0].data)
    if mini < 0:
        data[0].data = data[0].data - mini
    data.writeto('diff_files_temp'+os.sep+str(file),overwrite=True)
    data.close()
       
def hotpants(imag,imag_err,ref,ref_err):
        date = str(imag)[-31:-9]
        hotpants = ("./hotpants" +
                    " -inim "+str(imag)+" -ini "+str(imag_err)+" -iu 83000.00 -iuk 83000.00" +
                    " -tmplim "+str(ref)+" -tni "+str(ref_err)+" -tu 83000.00 -tuk 83000.00" +
                    " -nsx 5 -nsy 5 -afssc 0 -ssf stamps_lco.txt -bgo 2 -ng 3 4 0.8 3 2.4 2 4.0 -c t -ko 2" +
                    " -outim "+diffimages_folder+os.sep+"diff_"+str(date)+".fits -oni "+diffimages_folder+os.sep+"noise_"+str(date)+".fits" +
                    " -savexy "+diffimages_folder+os.sep+"stamps_"+str(date)+".reg")
        os.system(hotpants)

def hotpants_parallel(imag,imag_err,ref,ref_err,simnum):
        date = str(imag)[-31:-9]
        hotpants = ("./hotpants" +
                    " -inim "+str(imag)+" -ini "+str(imag_err)+" -iu 83000.00 -iuk 83000.00" +
                    " -tmplim "+str(ref)+" -tni "+str(ref_err)+" -tu 83000.00 -tuk 83000.00" +
                    " -nsx 5 -nsy 5 -afssc 0 -ssf stamps_lco.txt -bgo 2 -ng 3 4 0.8 3 2.4 2 4.0 -c t -ko 2" +
                    " -outim "+diffimages_folder+os.sep+"diff_"+str(date)+".fits -oni "+diffimages_folder+os.sep+"noise_"+str(date)+".fits" +
                    " -savexy "+diffimages_folder+os.sep+"stamps_"+str(date)+".reg &")
        os.system(hotpants)
        os.system('P'+str(simnum)+'=$!')
    
def make_stamptxt(imag,imag_err,ref,ref_err):
    date = str(imag)[-31:-9]
    hotpants = ("./hotpants" +
                " -inim "+str(imag)+" -ini "+str(imag_err)+" -iu 83000.00 -iuk 83000.00" +
                " -tmplim "+str(ref)+" -tni "+str(ref_err)+" -tu 83000.00 -tuk 83000.00" +
                " -nsx 5 -nsy 5 -afssc 1 -bgo 2 -ng 3 4 0.8 3 2.4 2 4.0 -c t -ko 2" +
                " -outim "+diffimages_folder+os.sep+"diff_"+str(date)+".fits -oni "+diffimages_folder+os.sep+"noise_"+str(date)+".fits" +
                " -savexy "+diffimages_folder+os.sep+"stamps_"+str(date)+".reg")
    os.system(hotpants)

print(' ')

local = (str(os.path.abspath(os.getcwd())))
path = Path(local)
print('local path:',path)

combi_fits_files = list(path.glob('combined'+os.sep+'c_a_t_*.fits'))
error_fits_files = list(path.glob('combi_error'+os.sep+'e_a_t_*.fits'))
print('number of total combi-fits-files:',len(combi_fits_files))
print('number of total error-fits-files:',len(error_fits_files))
    
if len(combi_fits_files) == len(error_fits_files):
    N = len(combi_fits_files)
else:
    print('## WARNING ## non matching number of c_*.fits and e_*.fits files!')
    print('## WARNING ## --> interupting galfitting.py')
    N = 0
    
min_adder('refimages','ref_'+ref_file_time+'.fits')

print(' ')

####

if make_stamps == True:
    print('####')
    print('#### making stamps:')
    print('####')
    for i in range(0,8,2): # use good image! --> try multiple times, here 4 times...
    
        stamp_file = str(combi_fits_files[i])
        error_file = stamp_file.replace('combined'+os.sep+'c_a_t','combi_error'+os.sep+'e_a_t')
        
        c_a_t_file = str(stamp_file.replace(str(path)+os.sep+'combined'+os.sep,''))
        min_adder('combined',c_a_t_file)
            
        image_file = 'diff_files_temp'+os.sep+str(c_a_t_file)
        error_file = str(error_file)
        ref_file = 'diff_files_temp'+os.sep+'ref_'+ref_file_time+'.fits'
        ref_error_file = 'refimages'+os.sep+'error_ref_'+ref_file_time+'.fits'
            
        make_stamptxt(image_file,error_file,ref_file,ref_error_file)

elif make_stamps == False:
    print('####')
    print('#### difference imaging:')
    print('####')
    if parallel == True:
            simnum = 0
    
    for day in range(N):
        print('####')
        print('#### day',str(day+1)+':')
        print('####')
        current_combi_file = str(combi_fits_files[day])
        current_error_file = current_combi_file.replace('combined'+os.sep+'c_a_t','combi_error'+os.sep+'e_a_t')
        #print(current_combi_file,current_error_file)
        if len([file for file in error_fits_files if current_error_file in str(file)]) != 1:
            print('####')
            print('#### warning: not able to match c_*.fits and e_*.fits of the day!')
            print('####')
            break
        c_date = str(current_combi_file)[-17:-9]
        e_date = str(current_error_file)[-17:-9]
        if c_date != e_date:
            print(c_date,e_date)
            break
         
        c_a_t_file = str(str(combi_fits_files[day]).replace(str(path)+os.sep+'combined'+os.sep,''))
        min_adder('combined',c_a_t_file)
         
        image_file = 'diff_files_temp'+os.sep+str(c_a_t_file)
        error_file = str(current_error_file)
        ref_file = 'diff_files_temp'+os.sep+'ref_'+ref_file_time+'.fits'
        ref_error_file = 'refimages'+os.sep+'error_ref_'+ref_file_time+'.fits'
               
        if parallel == True and day+1 < N:
            if simnum == 16:
                os.system('wait $P1 $P2 $P3 $P4 $P5 $P6 $P7 $P8 $P9 $P10 $P11 $P12 $P13 $P14 $P15 $P16')
                simnum = 0
            simnum = simnum + 1
            hotpants_parallel(image_file,error_file,ref_file,ref_error_file,simnum)
            
        else:
            hotpants(image_file,error_file,ref_file,ref_error_file)
        
        if day+1 == N:
            time.sleep(60) # waiting for 1 miute after the last image, to be sure all processes are finished
        
        print('####')
        print('#### finished day',day+1,'of',N)
        print('####')
        
else:
    print('####')
    print('#### Error!')
    print('####')
####

os.system('rm diff_files_temp'+os.sep+'*.fits') # removing the temp.-files
print(' ')
time.sleep(7) # 7 more seconds just for fun...
print('diffimage_hotpants.py is finished')