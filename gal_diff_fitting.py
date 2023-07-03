import numpy as np
import matplotlib.pyplot as plt

import os
import sys
from pathlib import Path

from astropy.io import fits
from astropy.stats import sigma_clipped_stats

##############################################################################
##############################################################################
# ONLY THING THAT HAS TO BE ADAPTED HERE AFTER SAYING MAKE IN GALFIT FOLDER! #
# (and the switch use_better_positions in line 41 to use better positions!)  #
##############################################################################
# Which quasar is worked on? --> Set that one to True, all others to false!
he1104 = False
he2149 = False
q2237 = True
##############################################################################
##############################################################################

local = (str(os.path.abspath(os.getcwd())))
path = Path(local)
print('local path:',path)

diff_fits_files = list(path.glob('diff*.fits')) # in galfitting: combi_fits_files --> here: diff-images
error_fits_files = list(path.glob('noise*.fits')) # corresponding error images to the diff-images
print('number of total combi-fits-files:',len(diff_fits_files))
print('number of total error-fits-files:',len(error_fits_files))

if len(diff_fits_files) == len(error_fits_files):
    N = len(diff_fits_files)
else:
    print('## WARNING ## non matching numbers of diff*.fits and noise*.fits files!')
    print('## WARNING ## --> interupting galfitting.py')
    N = 0
    sys.exit('non matching numbers of diff*.fits and noise*.fits files')

#######################################################################
use_better_positions = True
# if True: Better positions from betterquasarposition.txt from
#          the median of qsoA positions, corrected with psf star shift 
#          known from its pm from GAIA, plus linear psf star shift
#          again from its pm from GAIA, are used as image A positions.
#######################################################################
use_best_positions = True
# if True: If this is also true we even use the best positions, which
#          come from improved initial guesses not clever masking.
#          SO FAR ONLY FOR Q2237!!!
#######################################################################

#######################################################################
# do NOT change...
if q2237:
    four_images = True # if True: 4 quasar images; if False: 2 images #
else:
    four_images = False
#######################################################################

plotting = False # if True, plots of psfs will be shown

# NEW: in make_galfit_input_file muss in position nach den positions noch angegeben werden,
#      ob die Positionen festgehalten werden sollen (+' 0 0  ') oder ob sie wie davor optimiert werden sollen (+ 1 1  ').
def make_galfit_input_file(data_file_name,error_file_name,object_region,position,fluxesti,second_fit):
    galfit_input_txt = open("galfit.input", "w")
    galfit_input_txt.write("A) "+str(data_file_name)+" \n"+  
                           "B) out.fits "+"\n"+
                           "C) "+str(error_file_name)+" \n"+
                           "D) psf.fits "+"\n"+
                           "E) 1 "+"\n"+
                           "F)   "+"\n"+                  
                           "G)   "+"\n"+                     
                           "H) "+str(object_region)+" \n"+
                           "I) 51    51 "+"\n"+           
                           "J) 25.000 "+"\n"+           
                           "K) 0.387  0.387 "+"\n"+        
                           "O) regular "+"\n"+            
                           "P) 0 "+"\n"+                   
                           "S) 0 "+"\n"+                  
                           "0)      q2237  "+"\n"+       
                           " 1) "+str(position)+"\n"+
                           " 3) "+str(fluxesti)+"     1 "+"\n"+       
                           " 4) "+str(second_fit)+"\n"+       
                           " 5)      0.0   0 "+"\n"+       
                           " 6)      0.0   0 "+"\n"+       
                           " 7) 0.0000     0 "+"\n"+       
                           " 8) 1.0000     0 "+"\n"+        
                           " 9) 0.0000     0 "+"\n"+       
                           "10) 0.0000     0 "+"\n"+      
                           " Z) 0            "+"\n"+       
                           " 0)        sky   "+"\n"+      
                           " 1) 10.0       1 "+"\n"+       
                           " 2) 0.0000     0 "+"\n"+      
                           " 3) 0.0000     0 "+"\n"+      
                           " Z) 0 ")
    galfit_input_txt.close()

def make_galfit_input_file_4images(data_file_name,error_file_name,object_region,position,fluxesti,second_fit,third_fit,fourth_fit):
    galfit_input_txt = open("galfit.input", "w")
    galfit_input_txt.write("A) "+str(data_file_name)+" \n"+  
                           "B) out.fits "+"\n"+
                           "C) "+str(error_file_name)+" \n"+
                           "D) psf.fits "+"\n"+
                           "E) 1 "+"\n"+
                           "F)   "+"\n"+                  
                           "G)   "+"\n"+                     
                           "H) "+str(object_region)+" \n"+
                           "I) 51    51 "+"\n"+           
                           "J) 25.000 "+"\n"+           
                           "K) 0.387  0.387 "+"\n"+        
                           "O) regular "+"\n"+            
                           "P) 0 "+"\n"+                   
                           "S) 0 "+"\n"+                  
                           "0)      q2237  "+"\n"+       
                           " 1) "+str(position)+"\n"+           
                           " 3) "+str(fluxesti)+"     1 "+"\n"+       
                           " 4) "+str(second_fit)+"\n"+       
                           " 5) "+str(third_fit)+"\n"+       
                           " 6) "+str(fourth_fit)+"\n"+       
                           " 7) 0.0000     0 "+"\n"+       
                           " 8) 1.0000     0 "+"\n"+        
                           " 9) 0.0000     0 "+"\n"+       
                           "10) 0.0000     0 "+"\n"+      
                           " Z) 0            "+"\n"+       
                           " 0)        sky   "+"\n"+      
                           " 1) 10.0       1 "+"\n"+       
                           " 2) 0.0000     0 "+"\n"+      
                           " 3) 0.0000     0 "+"\n"+      
                           " Z) 0 ")
    galfit_input_txt.close()

def make_galfit_input_file_4images_with_galaxy(data_file_name,error_file_name,object_region,position,fluxesti,second_fit,third_fit,fourth_fit,galaxy_position):
    galfit_input_txt = open("galfit.input", "w")
    galfit_input_txt.write("A) "+str(data_file_name)+" \n"+  
                           "B) out.fits "+"\n"+
                           "C) "+str(error_file_name)+" \n"+
                           "D) psf.fits "+"\n"+
                           "E) 1 "+"\n"+
                           "F)   "+"\n"+                  
                           "G)   "+"\n"+                     
                           "H) "+str(object_region)+" \n"+
                           "I) 51    51 "+"\n"+           
                           "J) 25.000 "+"\n"+           
                           "K) 0.387  0.387 "+"\n"+        
                           "O) regular "+"\n"+            
                           "P) 0 "+"\n"+                   
                           "S) 0 "+"\n"+                  
                           " 0)      q2237  "+"\n"+       
                           " 1) "+str(position)+"\n"+ 
                           " 3) "+str(fluxesti)+"     1 "+"\n"+       
                           " 4) "+str(second_fit)+"\n"+       
                           " 5) "+str(third_fit)+"\n"+       
                           " 6) "+str(fourth_fit)+"\n"+       
                           " 7) 0.0000     0 "+"\n"+       
                           " 8) 1.0000     0 "+"\n"+        
                           " 9) 0.0000     0 "+"\n"+       
                           "10) 0.0000     0 "+"\n"+      
                           " Z) 0            "+"\n"+
                           "0) sersic        "+"\n"+
                           "1) "+str(galaxy_position)+" 1 1  "+"\n"+ # position
                           "3) 10.000      1          "+"\n"+ # integrated magnitude
                           "4) 4.7000      0          "+"\n"+ # R_e (half-light radius) [pix]
                           "5) 4.0000      0          "+"\n"+ # Sersic index n (de Vaucouleurs n=4)
                           "8) 0.6500      0          "+"\n"+ # axis ratio (b/a)
                           "9) 58.36       0          "+"\n"+ # position angle (PA) [deg: Up=0, Left=90]
                           "Z) 0                      "+"\n"+
                           "0) expdisk       "+"\n"+
                           "1) "+str(galaxy_position)+" 1 1  "+"\n"+ # position
                           "3) 10.000      1          "+"\n"+ # central surface brightness
                           "4) 10.400      0          "+"\n"+ # R_0 (exponetial scale length) [pix]
                           "8) 0.6500      0          "+"\n"+ # axis ratio (b/a)
                           "9) 58.36       0          "+"\n"+ # position angle (PA) [deg: Up=0, Left=90]
                           "Z) 0                      "+"\n"+
                           " 0)        sky   "+"\n"+      
                           " 1) 10.0       1 "+"\n"+     
                           " 2) 0.0000     0 "+"\n"+      
                           " 3) 0.0000     0 "+"\n"+      
                           " Z) 1 ")
    galfit_input_txt.close()

######################################### NEW ###########################################
# quasar positions from galfitting.py of the combined images (i.e. the non-diff.-images)

quapos_file = np.genfromtxt('quasarpositions.txt',skip_header=2,usecols=(0),dtype='str')
quapos_time,quapos_xpos,quapos_xpos_err,quapos_ypos,quapos_ypos_err = \
np.loadtxt('quasarpositions.txt',skiprows=2,usecols=(1,2,3,4,5),unpack=True)

if use_better_positions == True:
    better_quapos_file = np.genfromtxt('betterquasarposition.txt',skip_header=2,usecols=(0),dtype='str')
    better_quapos_time,better_quapos_xpos,better_quapos_xpos_err,better_quapos_ypos,better_quapos_ypos_err = \
    np.loadtxt('betterquasarposition.txt',skiprows=2,usecols=(1,2,3,4,5),unpack=True)        
    if (np.abs(quapos_time - better_quapos_time) < 0.01).all() and (quapos_file == better_quapos_file).all():
        # overwrite everything with better values:
        print('#### using better GAIA positions for fixed quasar image A position')
        if use_best_positions == True:
            print('#### using even the BEST!!!')
        quapos_file,quapos_time,quapos_xpos,quapos_xpos_err,quapos_ypos,quapos_ypos_err = \
        better_quapos_file,better_quapos_time,better_quapos_xpos,better_quapos_xpos_err,better_quapos_ypos,better_quapos_ypos_err
    else:
        sys.exit('data from quasarpositions.txt and betterquasarposition.txt do not fit!')
        
if use_best_positions == True:
    quapos_file = np.genfromtxt('bestquasarposition.txt',skip_header=2,usecols=(0),dtype='str')
    quapos_time,quapos_xpos,quapos_xpos_err,quapos_ypos,quapos_ypos_err = \
    np.loadtxt('bestquasarposition.txt',skiprows=2,usecols=(1,2,3,4,5),unpack=True) 
    print('#### overwritten better with best position data!!!')
    
#print('Quasar positions:')
#for i in range(len(quapos_time)):
#    print(quapos_file[i],quapos_time[i],quapos_xpos[i],quapos_xpos_err[i],quapos_ypos[i],quapos_ypos_err[i])
#print('These positions are now fixed for galfitting the difference images!')

#########################################################################################

# initilize data arrays for all N and both data values and errors:
quasar_posi = np.zeros((N,2,2)) # position of image A at day, values and errors, both x and y
if four_images == False:
    quasar_flux = np.zeros((N,2,2)) # quasar fluxes at day, values and errors, of images A and B
else:
    quasar_flux = np.zeros((N,2,4)) # quasar fluxes at day, values and errors, of images A, B, C and D
backgr_flux = np.zeros((N,2))   # sky/background flux at day, value and error
chi2red_val = np.zeros(N)       # chi-squared-reduced-values of quasar-galfit
residuals = np.zeros(N)         # array for mean value of residuals
residual_errors = np.zeros(N)   # array for std of mean of residuals

qso_posi_devi = np.ones(N) # is initially set to ones to see in the end 
                           # that the positions where kept fixed, i.e. that the deviation is zero!

# time array amd file list:
time = np.zeros(N)
file_list = []

################################################################
# galfitting of reference image for the difference lightcurves #
################################################################
galffiting_of_refimage = True
reference_images = list(path.glob('ref_20*.fits'))
ref_image_errors = list(path.glob('error_ref_20*.fits'))
refimageresults = False # DO NOT CHANGE!!!
if len(reference_images) == 1 and len(ref_image_errors) == 1 and galffiting_of_refimage:
    print('galfitting of ref-image')
    refimageresults = True
    reference_image = str(reference_images[0])[-16:]
    ref_image_error = str(ref_image_errors[0])[-22:]
    print(reference_image,ref_image_error)
    print('make psf.fits for ref-image')
    data = fits.open(reference_image)
    # remove possible remaining background from ref-image:
    mean,median,std = sigma_clipped_stats(data[0].data[100:-100,100:-100],sigma = 3.0,std_ddof=1)
    background = median
    error_background = std # std i.e. typical fluctuation size of the individual background pixels
    print('ref-image background =',background,'+/-',error_background) # Der Fehler des Hintergrunds wurde im ursprünglichen combining schon im Fehlerbild berücksichtigt (nur in den diff-Versuchen manchmal wieder das negative minimum im Bild dazuaddiert und abgespeichert)
    refimage_backred = data[0].data - background
    if q2237:
        psf_cutout = refimage_backred[2120:2150,2213:2243]
    elif he2149:
        psf_cutout = refimage_backred[1703:1733,1786:1816]
    elif he1104:
        psf_cutout = refimage_backred[1770:1800,1700:1730]
    else:
        sys.exit('Quasar-bool-error!')
    psf_fits = fits.PrimaryHDU(psf_cutout)
    psf_fits = fits.HDUList([psf_fits])
    psf_fits[0].data = psf_cutout
    psf_fits.writeto('psf.fits',overwrite=True)
    # copy input-file:
    psf_fits_str = 'cp psf.fits savelog'+os.sep+'psf_refimage.fits'
    psffits = fits.open('psf.fits')
    psf_cutout = psffits[0].data
    if plotting:
        plt.imshow(np.log10(abs(psf_cutout)+1.0),origin='lower')
        plt.title('psf of reference image')
        plt.colorbar()
        plt.show()
        plt.clf()
         
    print('make galfit.input file for ref-image')
    if q2237:
        quasar_region = '1709 1739 1737 1767'
        quasar_position = '1726.9 1755.0  1 1' #'1724.3 1754.9  1 1'
    elif he2149:
        quasar_region = '1689 1719 1693 1723'
        quasar_position = '1703.7 1708.3  1 1'
    elif he1104:
        quasar_region = '1712 1742 1732 1762'
        quasar_position = '1723.8 1745.6  1 1' 
    else:
        sys.exit('Quasar-bool-error!')
    fluxestimate = '100000.0'
    second_fit = '50000.0     1'
    if four_images == True:
        third_fit = '20000.0     1'
        fourth_fit = '10000.0     1'
        galaxy_position = '1726.7 1752.6' #'1725.5 1751.5'
        make_galfit_input_file_4images_with_galaxy(reference_image,ref_image_error,quasar_region,quasar_position,fluxestimate,second_fit,third_fit,fourth_fit,galaxy_position)
    else:
        make_galfit_input_file(reference_image,ref_image_error,quasar_region,quasar_position,fluxestimate,second_fit)
   
    # copy input-file:
    galinput_str = 'cp galfit.input savelog'+os.sep+'galfit_refimage.input'
    os.system(galinput_str)
    print('run galfit on ref-image')
    galfit_str = './galfit galfit.input > savelog'+os.sep+'conout_refimage'
    os.system(galfit_str) # execute galfit with no console-output, which is instead saved in conout_*
    os.system('rm galfit.??') # remove extra files produced by galfit
    # rename output-files:
    outfits_str = 'mv out.fits out_refimage.fits'
    os.system(outfits_str)
    fitlog_str = 'mv fit.log fit_refimage.log'
    os.system(fitlog_str)
    psffits_str = 'mv psf.fits psf_refimage.fits'
    os.system(psffits_str)
    print('read fit results of ref-image')
    if four_images == False:
        with open('fit_refimage.log') as f:
            # extract relevant information:
            galfit_results = f.readlines()[7:13]
        quasar_line = galfit_results[0].split()
        qu_err_line = galfit_results[1].split()
        backgr_line = galfit_results[2].split()
        ba_err_line = galfit_results[3].split()
        chi2nu_line = galfit_results[5].split()
    elif four_images == True:
        with open('fit_refimage.log') as f:
            # extract relevant information:
            galfit_results = f.readlines()[7:17]
        quasar_line = galfit_results[0].split()
        qu_err_line = galfit_results[1].split()
        backgr_line = galfit_results[6].split()
        ba_err_line = galfit_results[7].split()
        chi2nu_line = galfit_results[9].split()
    else:
        print('##########################ERROR1#############################')
    #print(quasar_line,qu_err_line,backgr_line,ba_err_line,chi2nu_line)
    if quasar_line[0] != 'q2237' or backgr_line[0] != 'sky':
        print('## WARNING ## error in extracting fit values and results!')
        print('## WARNING ## --> no ref-image-results!')
    else:
        if four_images == False:
            quasar_refim_xpos = float(quasar_line[2][1:-1])
            quasar_refim_ypos = float(quasar_line[3][0:-1])
            quasar_refim_deviation = np.sqrt((quasar_refim_xpos-float(quasar_position.split()[0]))**2+(quasar_refim_ypos-float(quasar_position.split()[1]))**2)
            quasar_refim_fluxA = float(quasar_line[4])
            quasar_refim_fluxB = float(quasar_line[5])
            quasar_refim_xposerr = float(qu_err_line[0][1:-1])
            quasar_refim_yposerr = float(qu_err_line[1][0:-1])
            quasar_refim_fluxAerr = float(qu_err_line[2])
            quasar_refim_fluxBerr = float(qu_err_line[3])
            backgr_refim_flux = float(backgr_line[4])
            backgr_refim_fluxerr = float(ba_err_line[0])
            refim_chi2red = float(chi2nu_line[2])        
        elif four_images == True:
            quasar_refim_xpos = float(quasar_line[2][1:-1])
            quasar_refim_ypos = float(quasar_line[3][0:-1])
            quasar_refim_deviation = np.sqrt((quasar_refim_xpos-float(quasar_position.split()[0]))**2+(quasar_refim_ypos-float(quasar_position.split()[1]))**2)
            quasar_refim_fluxA = float(quasar_line[4])
            quasar_refim_fluxB = float(quasar_line[5])
            quasar_refim_fluxC = float(quasar_line[6])
            quasar_refim_fluxD = float(quasar_line[7])
            quasar_refim_xposerr = float(qu_err_line[0][1:-1])
            quasar_refim_yposerr = float(qu_err_line[1][0:-1])
            quasar_refim_fluxAerr = float(qu_err_line[2])
            quasar_refim_fluxBerr = float(qu_err_line[3])
            quasar_refim_fluxCerr = float(qu_err_line[4])
            quasar_refim_fluxDerr = float(qu_err_line[5])
            backgr_refim_flux = float(backgr_line[4])
            backgr_refim_fluxerr = float(ba_err_line[0])
            refim_chi2red = float(chi2nu_line[2])
        else:
            print('##########################ERROR2#############################')
    # read out residuals from out*.fits:
    residual_data = fits.open('out_refimage.fits')
    residual_array = residual_data[3].data
    if plotting:
        plt.imshow(residual_array,origin='lower')
        plt.title('galfit-residuals of ref-image')
        plt.colorbar()
        plt.show()
        plt.clf()
    trim = 8
    refim_residual = np.mean(residual_array[trim:-trim,trim:-trim])
    refim_resi_err = np.std(residual_array[trim:-trim,trim:-trim],ddof=1)/(30-2*trim)
    # move output-files:
    outfits_str = 'mv out_refimage.fits savelog'+os.sep+'out_refimage.fits'
    os.system(outfits_str)
    fitlog_str = 'mv fit_refimage.log savelog'+os.sep+'fit_refimage.log'
    os.system(fitlog_str)
    psffits_str = 'mv psf_refimage.fits savelog'+os.sep+'psf_refimage.fits'
    os.system(psffits_str)
# accessing the ref-image results:
if refimageresults:
    print('ref-image results:')
    print('x-position:',quasar_refim_xpos,'+/-',quasar_refim_xposerr)
    print('y-position:',quasar_refim_ypos,'+/-',quasar_refim_yposerr)
    print('position-deviation:',quasar_refim_deviation)
    print('A-flux:',quasar_refim_fluxA,'+/-',quasar_refim_fluxAerr)
    print('B-flux:',quasar_refim_fluxB,'+/-',quasar_refim_fluxBerr)
    if four_images == True:
        print('C-flux:',quasar_refim_fluxC,'+/-',quasar_refim_fluxCerr)
        print('D-flux:',quasar_refim_fluxD,'+/-',quasar_refim_fluxDerr)
    print('background:',backgr_refim_flux,'+/-',backgr_refim_fluxerr)
    print('residual:',refim_residual,'+/-',refim_resi_err)
    print('chi-squared-reduced:',refim_chi2red)
    print('galfitting of ref-image finished!')    

#########################################################
######  NOW: galfitting, etc. of all diffimages:  #######
#########################################################

for day in range(N):
    ####
    print('####')
    print('#### day',str(day+1),'of',str(N))
    print('####')
    cdata_of_day = str(diff_fits_files[day])
    cdata_of_day = cdata_of_day.replace(str(path)+'/','')
    print(cdata_of_day)
    error_of_day = cdata_of_day.replace('diff','noise')
    if len([file for file in error_fits_files if error_of_day in str(file)]) != 1:
        print('## WARNING ## not able to match c_*.fits and e_*.fits of the day!')
        print('## WARNING ## --> skip day')
        break
    print(error_of_day)
    telescope = cdata_of_day[5:19]
    year = cdata_of_day[19:23]
    month = cdata_of_day[23:25]
    days = cdata_of_day[25:27]
    date = cdata_of_day[19:27]
    #print(cdata_of_day,telescope,year,month,days,date)
    time[day] = int(year)+(int(month)-1)/12+(int(days)-1)/365
    print('year of image:',time[day])
    file_list.append(cdata_of_day)    
    ####
    
    #### NEW: finde korrespondierenden Index in quapos-tabelle:
    quapos_index = 0
    for originalfilenames in quapos_file:
        eightdigitsdate = originalfilenames[-8:]
        if eightdigitsdate == date:
            #print('quapos index check by date and time comparison:')
            #print(eightdigitsdate,date)
            #print(quapos_time[quapos_index],time[day])
            print('quapos_index =',quapos_index)
            break
        quapos_index = quapos_index + 1
    #### NEW: quapos_index liefert jetzt das richtige Element aus quapos-listen/arrays!
    
    ####
    # psf kommt vom ursprünglichen galfit:
    psf_fits_str = 'cp savelog'+os.sep+'psf_'+str(telescope)+str(date)+'.fits psf.fits'
    print('psf_'+str(telescope)+str(date)+'.fits')
    os.system(psf_fits_str)
    psffits = fits.open('psf.fits')
    psf_cutout = psffits[0].data
    if plotting:
        plt.imshow(np.log10(abs(psf_cutout)+1.0),origin='lower')
        plt.title('psf of image no. '+str(day+1))
        plt.colorbar()
        plt.show()
        plt.clf()
    ####
    
    ####
    print('make galfit.input file')
    if q2237:
        quasar_region = '1709 1739 1737 1767'
    elif he2149:
        quasar_region = '1689 1719 1693 1723'
    elif he1104:
        quasar_region = '1712 1742 1732 1762'
    else:
        sys.exit('Quasar-bool-error!')
    quasar_position = str(quapos_xpos[quapos_index])+' '+str(quapos_ypos[quapos_index])+' 0 0  ' 
    fluxestimate = '0.0'
    second_fit = '0.0     1'
    if four_images == True:
        third_fit = '0.0     1'
        fourth_fit = '0.0     1'
        make_galfit_input_file_4images(cdata_of_day,error_of_day,quasar_region,quasar_position,fluxestimate,second_fit,third_fit,fourth_fit)
    else:
        make_galfit_input_file(cdata_of_day,error_of_day,quasar_region,quasar_position,fluxestimate,second_fit)
    # copy input-file:
    galinput_str = 'cp galfit.input savelog'+os.sep+'galfit_'+str(day+1)+'.input'
    os.system(galinput_str)
    ####
    
    ####
    print('run galfit')
    galfit_str = './galfit galfit.input > savelog'+os.sep+'conout_'+str(day+1)
    os.system(galfit_str) # execute galfit with no console-output, which is instead saved in conout_*
    os.system('rm galfit.??') # remove extra files produced by galfit
    # rename output-files:
    outfits_str = 'mv out.fits out_'+str(day+1)+'.fits'
    os.system(outfits_str)
    fitlog_str = 'mv fit.log fit_'+str(day+1)+'.log'
    os.system(fitlog_str)
    ####
    
    ####
    print('read fit results')
    if four_images == False:
        with open('fit_'+str(day+1)+'.log') as f:
            # extract relevant information:
            galfit_results = f.readlines()[7:13]
        quasar_line = galfit_results[0].split()
        qu_err_line = galfit_results[1].split()
        backgr_line = galfit_results[2].split()
        ba_err_line = galfit_results[3].split()
        chi2nu_line = galfit_results[5].split()
    elif four_images == True:
        with open('fit_'+str(day+1)+'.log') as f:
            # extract relevant information:
            galfit_results = f.readlines()[7:17]
        quasar_line = galfit_results[0].split()
        qu_err_line = galfit_results[1].split()
        backgr_line = galfit_results[2].split()
        ba_err_line = galfit_results[3].split()
        chi2nu_line = galfit_results[5].split()
        # this is actually the same, as there are no galaxy-output-lines!!!
    else:
        print('##########################ERROR3#############################')
    #print(quasar_line,qu_err_line,backgr_line,ba_err_line,chi2nu_line)
    
    # read out residuals from out*.fits:
    residual_data = fits.open('out_'+str(day+1)+'.fits')
    residual_array = residual_data[3].data
    if plotting:
        plt.imshow(residual_array,origin='lower')
        plt.title('galfit-residuals no. '+str(day+1))
        plt.colorbar()
        plt.show()
        plt.clf()
    trim = 8
    residuals[day] = np.mean(residual_array[trim:-trim,trim:-trim])
    residual_errors[day] = np.std(residual_array[trim:-trim,trim:-trim],ddof=1)/(30-2*trim)
    
    # move output-files:
    outfits_str = 'mv out_'+str(day+1)+'.fits savelog'+os.sep+'out_'+str(day+1)+'.fits'
    os.system(outfits_str)
    fitlog_str = 'mv fit_'+str(day+1)+'.log savelog'+os.sep+'fit_'+str(day+1)+'.log'
    os.system(fitlog_str)
    ####
    
    ####
    if four_images == False:
        if quasar_line[0] != 'q2237' or backgr_line[0] != 'sky':
            print('## WARNING ## error in extracting fit values and results!')
            print('## WARNING ## --> interupting galfitting.py')
            break
        else:
            quasar_posi[day,0,0] = float(quasar_line[2][1:-1])
            quasar_posi[day,0,1] = float(quasar_line[3][0:-1])
            quasar_flux[day,0,0] = float(quasar_line[4])
            quasar_flux[day,0,1] = float(quasar_line[5])
            quasar_posi[day,1,0] = float(qu_err_line[0][1:-1])
            quasar_posi[day,1,1] = float(qu_err_line[1][0:-1])
            quasar_flux[day,1,0] = float(qu_err_line[2])
            quasar_flux[day,1,1] = float(qu_err_line[3])
            backgr_flux[day,0] = float(backgr_line[4])
            backgr_flux[day,1] = float(ba_err_line[0])
            chi2red_val[day] = float(chi2nu_line[2])
        print('position:',quasar_posi[day,0],'+/-',quasar_posi[day,1])
        print('fluxes:',quasar_flux[day,0],'+/-',quasar_flux[day,1])
        print('background:',backgr_flux[day,0],'+/-',backgr_flux[day,1])
        print('chi-squared-reduced:',chi2red_val[day])
    elif four_images == True:
        if quasar_line[0] != 'q2237' or backgr_line[0] != 'sky':
            print('## WARNING ## error in extracting fit values and results!')
            print('## WARNING ## --> no ref-image-results!')
            break
        else:
            quasar_posi[day,0,0] = float(quasar_line[2][1:-1])
            quasar_posi[day,0,1] = float(quasar_line[3][0:-1])
            quasar_flux[day,0,0] = float(quasar_line[4])
            quasar_flux[day,0,1] = float(quasar_line[5])
            quasar_flux[day,0,2] = float(quasar_line[6])
            quasar_flux[day,0,3] = float(quasar_line[7])
            quasar_posi[day,1,0] = float(qu_err_line[0][1:-1])
            quasar_posi[day,1,1] = float(qu_err_line[1][0:-1])
            quasar_flux[day,1,0] = float(qu_err_line[2])
            quasar_flux[day,1,1] = float(qu_err_line[3])
            quasar_flux[day,1,2] = float(qu_err_line[4])
            quasar_flux[day,1,3] = float(qu_err_line[5])
            backgr_flux[day,0] = float(backgr_line[4])
            backgr_flux[day,1] = float(ba_err_line[0])
            chi2red_val[day] = float(chi2nu_line[2])
        print('position:',quasar_posi[day,0],'+/-',quasar_posi[day,1])
        print('fluxes:',quasar_flux[day,0],'+/-',quasar_flux[day,1])
        print('background:',backgr_flux[day,0],'+/-',backgr_flux[day,1])
        print('chi-squared-reduced:',chi2red_val[day])
    else:
        print('##########################ERROR4#############################')
    ####
    
    #### check whether position was kept fixed:
    position_shift = np.sqrt((quapos_xpos[quapos_index]-quasar_posi[day,0,0])**2+(quapos_ypos[quapos_index]-quasar_posi[day,0,1])**2)
    #print('Position shift:',position_shift)
    if position_shift > 0.01: # possibility of rounding error from putting the exact value into GALFIT, which seams to round to two decimal points...
        print('WARNING: position shift despite fixed input!')
        print(position_shift)
        print(quapos_xpos[quapos_index],quasar_posi[day,0,0],quapos_ypos[quapos_index],quasar_posi[day,0,1])
    qso_posi_devi[day] = position_shift
    ####
    
    # end of loop
    
####
print('###')
print('### final steps')
print('###')

position_deviation = qso_posi_devi
#print(position_deviation) # true deviations are 0.0 as they where kept fixed from last time!

# Achtung: no scales, because no ref-stars!

image_A = quasar_flux[:,0,0]
image_B = quasar_flux[:,0,1]
error_A = quasar_flux[:,1,0]
error_B = quasar_flux[:,1,1]
if four_images == True:
    image_C = quasar_flux[:,0,2]
    image_D = quasar_flux[:,0,3]
    error_C = quasar_flux[:,1,2]
    error_D = quasar_flux[:,1,3]
backflux = backgr_flux[:,0]
backerror =backgr_flux[:,1]

# chi-reduced-squared and residual statistics:
print('chi-squared-reduced: minimum:',np.min(chi2red_val),', median:',np.median(chi2red_val),'and maximum:',np.max(chi2red_val))
print('mean central residuals: minimum:',round(np.min(residuals),3),', median:',round(np.median(residuals),3),'and maximum:',round(np.max(residuals),3))

# text-file with all (relative to first of day) scales, background stds and weights:
if four_images == False:
    if use_better_positions == False:
        txt = open("diffgalfitresults.txt", "w")
    elif use_best_positions == True:
        txt = open("diffgalfitresults_bestpos.txt", "w")
    else:
        txt = open("diffgalfitresults_betterpos.txt", "w")
    txt.write("file"+"\t"+"\t"+"\t"+"time"+"\t"+
              "x-pos"+"\t"+"x-pos-err"+"\t"+"y-pos"+"\t"+"y-pos-err"+"\t"+"pos-dev"+"\t"+
              "A-flux"+"\t"+"A-err"+"\t"+"B-flux"+"\t"+"B-err"+"\t"+
              "sky"+"\t"+"sky-err"+"\t"+"residuals"+"\t"+"chi2red"+"\n")
    txt.write("\n")
    if refimageresults:
        txt.write(str(reference_image)+"\t"+str(0.0)+"\t"+
                  str(quasar_refim_xpos)+"\t"+str(quasar_refim_xposerr)+"\t"+str(quasar_refim_ypos)+"\t"+str(quasar_refim_yposerr)+"\t"+str(quasar_refim_deviation)+"\t"+
                  str(quasar_refim_fluxA)+"\t"+str(quasar_refim_fluxAerr)+"\t"+str(quasar_refim_fluxB)+"\t"+str(quasar_refim_fluxBerr)+"\t"+
                  str(backgr_refim_flux)+"\t"+str(backgr_refim_fluxerr)+"\t"+str(refim_residual)+"\t"+str(refim_chi2red)+"\n")
        txt.write("\n")
    for i in range(N):
        txt.write(str(file_list[i])+"\t"+str(round(time[i],4))+"\t"+
                  str(quasar_posi[i,0,0])+"\t"+str(quasar_posi[i,1,0])+"\t"+str(quasar_posi[i,0,1])+"\t"+str(quasar_posi[i,1,1])+"\t"+str(round(position_deviation[i],4))+"\t"+
                  str(round(image_A[i],4))+"\t"+str(round(error_A[i],4))+"\t"+str(round(image_B[i],4))+"\t"+str(round(error_B[i],4))+"\t"+
                  str(round(backflux[i],3))+"\t"+str(round(backerror[i],3))+"\t"+str(round(residuals[i],3))+"\t"+str(chi2red_val[i])+"\n")
    txt.close()
else:
    if use_better_positions == False:
        txt = open("diffgalfitresults.txt", "w")
    elif use_best_positions == True:
        txt = open("diffgalfitresults_bestpos.txt", "w")
    else:
        txt = open("diffgalfitresults_betterpos.txt", "w")
    txt.write("file"+"\t"+"\t"+"\t"+"time"+"\t"+
              "x-pos"+"\t"+"x-pos-err"+"\t"+"y-pos"+"\t"+"y-pos-err"+"\t"+"pos-dev"+"\t"+
              "A-flux"+"\t"+"A-err"+"\t"+"B-flux"+"\t"+"B-err"+"\t"+"C-flux"+"\t"+"C-err"+"\t"+"D-flux"+"\t"+"D-err"+"\t"+
              "sky"+"\t"+"sky-err"+"\t"+"residuals"+"\t"+"chi2red"+"\n")
    txt.write("\n")
    if refimageresults:
        txt.write(str(reference_image)+"\t"+str(0.0)+"\t"+
                  str(quasar_refim_xpos)+"\t"+str(quasar_refim_xposerr)+"\t"+str(quasar_refim_ypos)+"\t"+str(quasar_refim_yposerr)+"\t"+str(quasar_refim_deviation)+"\t"+
                  str(quasar_refim_fluxA)+"\t"+str(quasar_refim_fluxAerr)+"\t"+str(quasar_refim_fluxB)+"\t"+str(quasar_refim_fluxBerr)+"\t"+
                  str(quasar_refim_fluxC)+"\t"+str(quasar_refim_fluxCerr)+"\t"+str(quasar_refim_fluxD)+"\t"+str(quasar_refim_fluxDerr)+"\t"+
                  str(backgr_refim_flux)+"\t"+str(backgr_refim_fluxerr)+"\t"+str(refim_residual)+"\t"+str(refim_chi2red)+"\n")
        txt.write("\n")
    for i in range(N):
        txt.write(str(file_list[i])+"\t"+str(round(time[i],4))+"\t"+
                  str(quasar_posi[i,0,0])+"\t"+str(quasar_posi[i,1,0])+"\t"+str(quasar_posi[i,0,1])+"\t"+str(quasar_posi[i,1,1])+"\t"+str(round(position_deviation[i],4))+"\t"+
                  str(round(image_A[i],4))+"\t"+str(round(error_A[i],4))+"\t"+str(round(image_B[i],4))+"\t"+str(round(error_B[i],4))+"\t"+
                  str(round(image_C[i],4))+"\t"+str(round(error_C[i],4))+"\t"+str(round(image_D[i],4))+"\t"+str(round(error_D[i],4))+"\t"+
                  str(round(backflux[i],3))+"\t"+str(round(backerror[i],3))+"\t"+str(round(residuals[i],3))+"\t"+str(chi2red_val[i])+"\n")
    txt.close()

# check all results and mask potentially bad ones regarding position-deviation, chi2red, residual and background:
exclude_counter = 0
exclude_mask = np.zeros(N)
for i in range(N):
    if position_deviation[i] > 1.0 or chi2red_val[i] > 10.0 or residuals[i] > 40.0 or abs(backflux[i]) > 50.0:
        # fit to initial distance larger than 1.0 pixel --> seems wrong
        print('## WARNING ## exclude data point from file:',file_list[i],'because:')
        print('## WARNING ## either the fitted position deviates strong from initial guess:',round(position_deviation[i],4))
        print('## WARNING ## or chi-squared-reduced is way too high:',chi2red_val[i])
        print('## WARNING ## or the background deviates unreasonably from zero:',round(backflux[i],4))
        print('## WARNING ## or the residual is too high:',round(residuals[i],4))
        exclude_mask[i] = 1
        exclude_counter = exclude_counter + 1
print(exclude_counter,'days excluded!')

# masking the arrays that are being plotted:
time = np.ma.array(time,mask=exclude_mask)
image_A = np.ma.array(image_A,mask=exclude_mask)
image_B = np.ma.array(image_B,mask=exclude_mask)
error_A = np.ma.array(error_A,mask=exclude_mask)
error_B = np.ma.array(error_B,mask=exclude_mask)
if four_images == True:
    image_C = np.ma.array(image_C,mask=exclude_mask)
    image_D = np.ma.array(image_D,mask=exclude_mask)
    error_C = np.ma.array(error_C,mask=exclude_mask)
    error_D = np.ma.array(error_D,mask=exclude_mask)
backflux = np.ma.array(backflux,mask=exclude_mask)
backerror = np.ma.array(backerror,mask=exclude_mask)
chi2red_val = np.ma.array(chi2red_val,mask=exclude_mask)
residuals = np.ma.array(residuals,mask=exclude_mask)
residual_errors = np.ma.array(residual_errors,mask=exclude_mask)
position_deviation = np.ma.array(position_deviation,mask=exclude_mask)
####

####
# plotting:

if q2237:
    title = 'lightcurves of quasar Q2237+0305'
elif he2149:
    title = 'lightcurves of quasar HE2149-2745'
elif he1104:
    title = 'lightcurves of quasar Q2237-1805'
else:
    sys.exit('Quasar-bool-error!')
    
# lightcurves plot with linear y-axis:
plt.errorbar(time,image_A,yerr=error_A,fmt='o',markersize=2,elinewidth=1,capsize=3,label='image A')
plt.errorbar(time,image_B,yerr=error_B,fmt='o',markersize=2,elinewidth=1,capsize=3,label='image B')
if four_images == True:
    plt.errorbar(time,image_C,yerr=error_C,fmt='o',markersize=2,elinewidth=1,capsize=3,label='image C')
    plt.errorbar(time,image_D,yerr=error_D,fmt='o',markersize=2,elinewidth=1,capsize=3,label='image D')
plt.title(title)
plt.legend()
if use_better_positions == False:
    plt.savefig('lightcurves_lin.jpg',dpi=800)
elif use_best_positions == True:
    plt.savefig('bestpos_lightcurves_lin.jpg',dpi=800)
else:
    plt.savefig('betterpos_lightcurves_lin.jpg',dpi=800)
if plotting:
    plt.show()
plt.clf()
# lightcurves plot with logarithmic y-axis:
plt.errorbar(time,image_A,yerr=error_A,fmt='o',markersize=2,elinewidth=1,capsize=3,label='image A')
plt.errorbar(time,image_B,yerr=error_B,fmt='o',markersize=2,elinewidth=1,capsize=3,label='image B')
if four_images == True:
    plt.errorbar(time,image_C,yerr=error_C,fmt='o',markersize=2,elinewidth=1,capsize=3,label='image C')
    plt.errorbar(time,image_D,yerr=error_D,fmt='o',markersize=2,elinewidth=1,capsize=3,label='image D')
plt.yscale('log')
plt.legend()
plt.title(title)
if use_better_positions == False:
    plt.savefig('lightcurves_log.jpg',dpi=800)
elif use_best_positions == True:
    plt.savefig('bestpos_lightcurves_log.jpg',dpi=800)
else:
    plt.savefig('betterpos_lightcurves_log.jpg',dpi=800)
if plotting:
    plt.show()
plt.clf()
# lightcurves plot of 2015-2018:
plt.errorbar(time,image_A,yerr=error_A,fmt='o',markersize=2,elinewidth=1,capsize=3,label='image A')
plt.errorbar(time,image_B,yerr=error_B,fmt='o',markersize=2,elinewidth=1,capsize=3,label='image B')
if four_images == True:
    plt.errorbar(time,image_C,yerr=error_C,fmt='o',markersize=2,elinewidth=1,capsize=3,label='image C')
    plt.errorbar(time,image_D,yerr=error_D,fmt='o',markersize=2,elinewidth=1,capsize=3,label='image D')
plt.xlim(2013.9,2016.4)
plt.title(title)
plt.legend()
if use_better_positions == False:
    plt.savefig('lightcurves_part1.jpg',dpi=800)
elif use_best_positions == True:
    plt.savefig('bestpos_lightcurves_part1.jpg',dpi=800)
else:
    plt.savefig('betterpos_lightcurves_part1.jpg',dpi=800)
if plotting:
    plt.show()
plt.clf()
# lightcurves plot of 2018-2022:
plt.errorbar(time,image_A,yerr=error_A,fmt='o',markersize=2,elinewidth=1,capsize=3,label='image A')
plt.errorbar(time,image_B,yerr=error_B,fmt='o',markersize=2,elinewidth=1,capsize=3,label='image B')
if four_images == True:
    plt.errorbar(time,image_C,yerr=error_C,fmt='o',markersize=2,elinewidth=1,capsize=3,label='image C')
    plt.errorbar(time,image_D,yerr=error_D,fmt='o',markersize=2,elinewidth=1,capsize=3,label='image D')
plt.xlim(2018.0,2023.2)
plt.title(title)
plt.legend()
if use_better_positions == False:
    plt.savefig('lightcurves_part2.jpg',dpi=800)
elif use_best_positions == True:
    plt.savefig('bestpos_lightcurves_part2.jpg',dpi=800)
else:
    plt.savefig('betterpos_lightcurves_part2.jpg',dpi=800)
if plotting:
    plt.show()
plt.clf()
# background, chi-reduced-squared, residual and position-deviation plot:
plt.errorbar(time,backflux,yerr=backerror,fmt='o',markersize=2,elinewidth=1,capsize=3,label='background flux')
plt.errorbar(time,residuals,yerr=residual_errors,fmt='o',markersize=2,elinewidth=1,capsize=3,label='mean central residual (fit vs. data)')
plt.scatter(time,position_deviation,marker='x',color='tab:red',label='position deviation (fit vs. intitial)')
plt.scatter(time,chi2red_val-1.0,marker='+',color='tab:purple',label='$\chi_{red}^2-1$')
plt.legend(fontsize=7)
plt.title('additional information')
if use_better_positions == False:
    plt.savefig('lightcurves_supplement.jpg',dpi=800)
elif use_best_positions == True:
    plt.savefig('bestpos_lightcurves_supplement.jpg',dpi=800)
else:
    plt.savefig('betterpos_lightcurves_supplement.jpg',dpi=800)
if plotting:
    plt.show()
plt.clf()
# NEW: difflightcurve with offset from refimage:
if refimageresults:
    plt.errorbar(time,image_A+quasar_refim_fluxA,yerr=error_A,fmt='o',markersize=2,elinewidth=1,capsize=3,label='image A')
    plt.errorbar(time,image_B+quasar_refim_fluxB,yerr=error_B,fmt='o',markersize=2,elinewidth=1,capsize=3,label='image B')
    if four_images == True:
        plt.errorbar(time,image_C+quasar_refim_fluxC,yerr=error_C,fmt='o',markersize=2,elinewidth=1,capsize=3,label='image C')
        plt.errorbar(time,image_D+quasar_refim_fluxD,yerr=error_D,fmt='o',markersize=2,elinewidth=1,capsize=3,label='image D')
    plt.title(title)
    plt.legend()
    if use_better_positions == False:
        plt.savefig('lightcurves_diffwithref.jpg',dpi=800)
    else:
        plt.savefig('betterpos_lightcurves_diffwithref.jpg',dpi=800)
    if plotting:
        plt.show()
    plt.clf() 
    # and now with logarithmic y-axis:
    plt.errorbar(time,image_A+quasar_refim_fluxA,yerr=error_A,fmt='o',markersize=2,elinewidth=1,capsize=3,label='image A')
    plt.errorbar(time,image_B+quasar_refim_fluxB,yerr=error_B,fmt='o',markersize=2,elinewidth=1,capsize=3,label='image B')
    if four_images == True:
        plt.errorbar(time,image_C+quasar_refim_fluxC,yerr=error_C,fmt='o',markersize=2,elinewidth=1,capsize=3,label='image C')
        plt.errorbar(time,image_D+quasar_refim_fluxD,yerr=error_D,fmt='o',markersize=2,elinewidth=1,capsize=3,label='image D')
    plt.yscale('log')
    plt.title(title)
    plt.legend()
    if use_better_positions == False:
        plt.savefig('lightcurves_diffwithref_log.jpg',dpi=800)
    elif use_best_positions == True:
        plt.savefig('bestpos_lightcurves_diffwithref_log.jpg',dpi=800)
    else:
        plt.savefig('betterpos_lightcurves_diffwithref_log.jpg',dpi=800)
    if plotting:
        plt.show()
    plt.clf() 
####

print('gal_diff_fitting.py is finished!')