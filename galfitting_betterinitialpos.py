import numpy as np
import matplotlib.pyplot as plt

import sys
import os
from pathlib import Path

from astropy.io import fits

local = (str(os.path.abspath(os.getcwd())))
path = Path(local)
print('local path:',path)

combi_fits_files = list(path.glob('c_a_t_*.fits'))
error_fits_files = list(path.glob('e_a_t_*.fits'))
print('number of total combi-fits-files:',len(combi_fits_files))
print('number of total error-fits-files:',len(error_fits_files))

if len(combi_fits_files) == len(error_fits_files):
    N = len(combi_fits_files)
else:
    print('## WARNING ## non matching number of c_*.fits and e_*.fits files!')
    print('## WARNING ## --> interupting galfitting.py')
    N = 0

#N = 2 # to test code with first few day(s)

four_images_with_galaxy = True # not only 2 images, but 4 plus galaxy, e.g.: q2237! THIS CODE WAS MADE FOR Q2237 --> SHOULD BE SET TO TRUE!!!!

refstar_position_output = True # if True, table with the refstarpostions is made for further checks

plotting = False # if True, plots of psfs, quasar regions and residuals will be shown

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
                           " 1) "+str(position)+" 1 1  "+"\n"+ 
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

# For looking at the model components, use Z: 'Skip this model in output image? (yes=1, no=0)'!
# Set P to 1 to create only the model as output!
# Then look at model: out_1.fits[2] but NOT in savelog but in the galfit-directory itself!
    
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
                           " 1) "+str(position)+" 1 1  "+"\n"+ 
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

###################################################################### NEW!!!!!!! ######################################################################

# use better positions from improved positions (line from GAIA-slope and masked median offset) 
# as new intial position guesses to improve the galfit-fit to improve the positions, 
# to use them to get better positions with no masking and to therefor use these best positions for plotting and galdifffitting!

use_better_positions = True #### MUST BE TRUE!!! --> THIS IS THE IDEA OF THIS CODE!!! DOESNT WORK OTHERWISE...
if use_better_positions == False:
    print('WARNING!!!')
    print('WARNING!!!')
    print('WARNING!!!')
    print('WARNING!!!')
    print('WARNING!!!')
    print('WARNING!!!')
    print('WARNING!!!')
if use_better_positions == True:
    better_quapos_file = np.genfromtxt('betterquasarposition.txt',skip_header=2,usecols=(0),dtype='str')
    better_quapos_time,better_quapos_xpos,better_quapos_xpos_err,better_quapos_ypos,better_quapos_ypos_err = \
    np.loadtxt('betterquasarposition.txt',skiprows=2,usecols=(1,2,3,4,5),unpack=True)
    print('#### using the better quasar positions as initial guesses for the quasar image A position')
    print('#### where we shift the positions randomly a small bit to avoid problems with galfit...')
    quapos_file,quapos_time,quapos_xpos,quapos_xpos_err,quapos_ypos,quapos_ypos_err = \
    better_quapos_file,better_quapos_time,better_quapos_xpos+(2.0*np.random.random(len(better_quapos_xpos))-1.0)*0.1,better_quapos_xpos_err,better_quapos_ypos+(2.0*np.random.random(len(better_quapos_ypos))-1.0)*0.1,better_quapos_ypos_err
    if False:
        print(better_quapos_xpos[0:5],better_quapos_ypos[0:5])
        print(quapos_xpos[0:5],quapos_ypos[0:5])
        print(quapos_file[0:5],quapos_time[0:5])
        sys.exit('test')
    
###################################################################### NEW!!!!!!! ######################################################################

# initilize data arrays for all N and both data values and errors:
quasar_posi = np.zeros((N,2,2)) # position of image A at day, values and errors, both x and y
if four_images_with_galaxy == False:
    quasar_flux = np.zeros((N,2,2)) # quasar fluxes at day, values and errors, of images A and B
else:
    quasar_flux = np.zeros((N,2,4)) # quasar fluxes at day, values and errors, of images A, B, C and D
    galaxy_results = np.zeros((N,2,6)) # results and errors of galaxy fits: 4 positions and 2 magnitudes
backgr_flux = np.zeros((N,2))   # sky/background flux at day, value and error
chi2red_val = np.zeros(N)       # chi-squared-reduced-values of quasar-galfit
residuals = np.zeros(N)         # array for mean value of residuals
residual_errors = np.zeros(N)   # array for std of mean of residuals

# time array amd file list:
time = np.zeros(N)
file_list = []

# ref-stars for image scaling:
# initial guesses for star positions for aperture photometry for scale factors:
# q2237:   
initial_positions = np.array([[1578.0,1825.0],[1600.0,1586.0],[2066.0,1498.0],[2203.0,1689.0],[1890.0,1904.0],
                              [1837.0,2049.0],[1677.0,1991.0],[1271.0,1603.0],[2191.0,2213.0],[2421.0,2327.0],
                              [1721.0,2276.0],[1505.0,2288.0],[1339.0,2281.0],[1052.0,1921.0],[ 812.0,1835.0],
                              [ 874.0,1565.0],[1123.0,1287.0],[1408.0,1368.0],[2055.0,1346.0],[2732.0,1166.0],
                              [3266.0,2570.0]])
# he2149:
#initial_positions = np.array([[2013.0,1854.0],[1274.0,1823.0],[1280.0,1650.0],[2150.0,1518.0],[2239.0,1552.0],
#                              [2199.0,1889.0],[2482.0,1831.0],[1851.0,1891.0],[1086.0,1854.0],[1043.0,1618.0],
#                              [ 993.0,1574.0],[ 959.0,1523.0],[1190.0,1986.0],[ 507.0,1928.0],[ 656.0,1384.0],
#                              [1602.0,1381.0],[1823.0,1359.0],[2387.0,1375.0],[2898.0,2169.0],[2733.0,2107.0],
#                              [2318.0,2246.0],[1793.0,2304.0],[1595.0,1049.0],[1504.0, 763.0]])
# he1104:
#initial_positions = np.array([[1816.0,1649.0],[1847.0,1657.0],[2080.0,1895.0],[2159.0,1977.0],[1354.0,1910.0],
#                              [1210.0,1983.0],[1662.0, 882.0],[ 977.0, 900.0],[2814.0,1040.0],[3018.0,1638.0],
#                              [2849.0,2156.0],[2228.0,2063.0],[1742.0,2317.0],[1304.0,2398.0],[ 364.0,2267.0],
#                              [ 992.0,1694.0],[1087.0,1443.0],[1826.0,1707.0],[1805.0,1760.0],[2342.0, 729.0],
#                              [ 782.0,1648.0],[3082.0, 737.0]])

ref_positions = [str(np.array(x,dtype=int))[1:-1] for x in initial_positions.tolist()] #converts to ['abcd wxyz',...] where abcd and wxyz are the ds9 pixel positions (in the order given by ds9):
                 
#ref_positions = ref_positions[:1] # to test quickly
#print(ref_positions)

#ref_regions:
halfside = 15 # size (half side) of region around ref-stars
ref_regions = []
for position in ref_positions:
    xy_pos = position.split()
    center = np.floor(np.array([float(xy_pos[0]),float(xy_pos[1])])).astype(int)
    edge1,edge2,edge3,edge4 = center[0]-halfside,center[0]+halfside,center[1]-halfside,center[1]+halfside
    region = str(edge1)+' '+str(edge2)+' '+str(edge3)+' '+str(edge4)
    ref_regions.append(region)
#print(ref_regions)

ref_fluxes = np.zeros((N,len(ref_regions)))
refstarpos = np.zeros((N,len(ref_regions),4))

for day in range(N):
    ####
    print('####')
    print('#### day',str(day+1),'of',str(N))
    print('####')
    cdata_of_day = str(combi_fits_files[day])
    cdata_of_day = cdata_of_day.replace(str(path)+'/','')
    print(cdata_of_day)
    error_of_day = cdata_of_day.replace('c_a_t','e_a_t')
    if len([file for file in error_fits_files if error_of_day in str(file)]) != 1:
        print('## WARNING ## not able to match c_*.fits and e_*.fits of the day!')
        print('## WARNING ## --> skip day')
        break
    print(error_of_day)
    year = cdata_of_day[20:24]
    month = cdata_of_day[24:26]
    days = cdata_of_day[26:28]
    date = cdata_of_day[20:28]
    time[day] = int(year)+(int(month)-1)/12+(int(days)-1)/365
    print('year of image:',time[day])
    file_list.append(cdata_of_day[6:28])    
    ####
    
    ########################################################################################
    #### NEW: finde korrespondierenden Index in quapos-tabelle:
    quapos_index = 0
    for originalfilenames in quapos_file:
        eightdigitsdate = originalfilenames[-8:]
        if eightdigitsdate == date:
            break
        quapos_index = quapos_index + 1
    print('## quapos index date and time check (betterpos-table vs. c_a_t-file):')
    print('##',eightdigitsdate,date)
    print('##',quapos_time[quapos_index],time[day])
    print('## quapos_index =',quapos_index)
    #### NEW: quapos_index liefert jetzt das richtige Element aus quapos-listen/arrays!    
    ########################################################################################
    
    ####
    print('make psf.fits')
    data = fits.open(cdata_of_day)
    ######################################################################################################
    # q2237:
    psf_cutout = data[0].data[2120:2150,2213:2243]
    # he2149:
    #psf_cutout = data[0].data[1703:1733,1786:1816]
    # he1104: 
    #psf_cutout = data[0].data[1770:1800,1700:1730]
    ######################################################################################################
    if plotting:
        plt.imshow(psf_cutout,origin='lower')
        plt.title('psf of image no. '+str(day+1))
        plt.colorbar()
        plt.show()
        plt.clf()
    psf_fits = fits.PrimaryHDU(psf_cutout)
    psf_fits = fits.HDUList([psf_fits])
    psf_fits[0].data = psf_cutout
    psf_fits.writeto('psf.fits',overwrite=True)
    # copy input-file:
    psf_fits_str = 'cp psf.fits savelog'+os.sep+'psf_'+str(day+1)+'.fits'
    os.system(psf_fits_str)
    ####
    
    ####
    print('make galfit.input file')
    ######################################################################################################
    # q2237:
    quasar_region = '1709 1739 1737 1767'
    quasar_position = '1724.3 1754.9'
    # he2149:
    #quasar_region = '1689 1719 1693 1723'
    #quasar_position = '1703.7 1708.3'
    # he1104:
    #quasar_region = '1712 1742 1732 1762' 
    #quasar_position = '1723.8 1745.6' 
    ######################################################################################################
    if plotting:
        qrsplit = [int(x) for x in quasar_region.split(' ')]
        plt.imshow(data[0].data[qrsplit[2]:qrsplit[3],qrsplit[0]:qrsplit[1]],origin='lower')
        plt.title('quasar cutout of image no. '+str(day+1))
        plt.colorbar()
        plt.show()
        plt.clf()
        
    fluxestimate = '50000.0'
    second_fit = '25000.0     1'
    
    if use_better_positions == True:
        # override quasar positions with better positions:
        print('## quapos before:',quasar_position)
        quasar_position = str(quapos_xpos[quapos_index])+' '+str(quapos_ypos[quapos_index])
        print('## quapos after:',quasar_position)
    
    if four_images_with_galaxy == True:
        third_fit = '20000.0     1'
        fourth_fit = '10000.0     1'
        galaxy_position = '1724.5 1752.5' # q2237
        
        if use_better_positions == True:
            # override galaxy position with better positions + HST quasar to galaxy distance :
            print('## galpos before:',galaxy_position)
            galaxy_position = str(quapos_xpos[quapos_index]-0.19)+' '+str(quapos_ypos[quapos_index]-2.43)
            print('## galpos after:',galaxy_position)
            
            print('## initial guesses for quasar and galaxy positions improved!')
            
        make_galfit_input_file_4images_with_galaxy(cdata_of_day,error_of_day,quasar_region,quasar_position,fluxestimate,second_fit,third_fit,fourth_fit,galaxy_position)
    
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
    if four_images_with_galaxy == False:
        with open('fit_'+str(day+1)+'.log') as f:
            # extract relevant information:
            galfit_results = f.readlines()[7:13]
        quasar_line = galfit_results[0].split()
        qu_err_line = galfit_results[1].split()
        backgr_line = galfit_results[2].split()
        ba_err_line = galfit_results[3].split()
        chi2nu_line = galfit_results[5].split()
    elif four_images_with_galaxy == True:
        with open('fit_'+str(day+1)+'.log') as f:
            # extract relevant information:
            galfit_results = f.readlines()[7:17]
        quasar_line = galfit_results[0].split()
        qu_err_line = galfit_results[1].split()
        backgr_line = galfit_results[6].split()
        ba_err_line = galfit_results[7].split()
        chi2nu_line = galfit_results[9].split()
        # additional galaxy-fit-results:
        galaxy_lines = []
        galaxy_lines.append(galfit_results[2].split())
        galaxy_lines.append(galfit_results[3].split())
        galaxy_lines.append(galfit_results[4].split())
        galaxy_lines.append(galfit_results[5].split())
        #print(galaxy_lines)
    else:
        break
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
    
    #### REFERENCE STARS ####
    print('reference star measurements')
    for ref_star in range(len(ref_regions)):
        refstar_region = ref_regions[ref_star]
        refstar_position = ref_positions[ref_star]
        no_second_fit = '     0.0   0'
        fluxestimate = '500000.0'
        make_galfit_input_file(cdata_of_day,error_of_day,refstar_region,refstar_position,fluxestimate,no_second_fit)
        os.system('./galfit galfit.input > /dev/null') # execute galfit with no console-output
        os.system('rm galfit.??') # remove extra files produced by galfit
        os.system('rm out.fits') # remove output fits file
        try:
            with open('fit.log') as f:
                galfit_ref_results = f.readlines()[7:9]
                #print(galfit_ref_results)
            os.system('rm fit.log') # remove fit.log
        except:
            print('refstar-galfit-error with refstar',str(ref_star+1),'! --> skipped!')
            continue
        ref_star_flux = galfit_ref_results[0].split()[4]
        if ')' in ref_star_flux:
            ref_fluxes[day,ref_star] = galfit_ref_results[0].split()[5]
        elif ',' in ref_star_flux:
            ref_fluxes[day,ref_star] = galfit_ref_results[0].split()[6]
        else:
            ref_fluxes[day,ref_star] = ref_star_flux
        # now refstar positions if needed:
        if refstar_position_output == True:
            #print(galfit_ref_results[0].split(),galfit_ref_results[1].split())
            if galfit_ref_results[0].split()[2] == '(':
                refstarpos[day,ref_star,0] = galfit_ref_results[0].split()[3].replace(',','') # x ref pos
                refstarpos[day,ref_star,1] = galfit_ref_results[1].split()[0].replace('(','').replace(',','') # x pos err
                refstarpos[day,ref_star,2] = galfit_ref_results[0].split()[4].replace(')','') # y ref pos
                refstarpos[day,ref_star,3] = galfit_ref_results[1].split()[1].replace(')','') # y pos err
            else:
                refstarpos[day,ref_star,0] = galfit_ref_results[0].split()[2].replace('(','').replace(',','') # x ref pos
                refstarpos[day,ref_star,1] = galfit_ref_results[1].split()[0].replace('(','').replace(',','') # x pos err
                refstarpos[day,ref_star,2] = galfit_ref_results[0].split()[3].replace(')','') # y ref pos
                refstarpos[day,ref_star,3] = galfit_ref_results[1].split()[1].replace(')','') # y pos err
            #print(refstarpos[day,ref_star])
    #### REFERENCE STARS ####
    
    ####
    if quasar_line[0] != 'q2237' or backgr_line[0] != 'sky':
        print('## WARNING ## error in extracting fit values and results!')
        print('## WARNING ## --> interupting galfitting.py')
        break
    else:
        if four_images_with_galaxy == False:
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
        elif four_images_with_galaxy == True:
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
            # additional galaxy-fit-results:
            galaxy_results[day,0,0]=float(galaxy_lines[0][2][1:-1]) # sersic x position
            galaxy_results[day,0,1]=float(galaxy_lines[0][3][0:-1]) # sersic y position
            galaxy_results[day,0,2]=float(galaxy_lines[0][4])       # sersic magnitude
            galaxy_results[day,0,3]=float(galaxy_lines[2][2][1:-1]) # expdisc x position
            galaxy_results[day,0,4]=float(galaxy_lines[2][3][0:-1]) # expdisc y position
            galaxy_results[day,0,5]=float(galaxy_lines[2][4])       # expdisc magnitude
            galaxy_results[day,1,0]=float(galaxy_lines[1][0][1:-1]) # error of sersic x position
            galaxy_results[day,1,1]=float(galaxy_lines[1][1][0:-1]) # error of sersic y position
            galaxy_results[day,1,2]=float(galaxy_lines[1][2])       # error of sersic magnitude
            galaxy_results[day,1,3]=float(galaxy_lines[3][0][1:-1]) # error of expdisc x position
            galaxy_results[day,1,4]=float(galaxy_lines[3][1][0:-1]) # error of expdisc y position
            galaxy_results[day,1,5]=float(galaxy_lines[3][2])       # error of expdisc magnitude
            #print(galaxy_results)
        else:
            break
    print('position:',quasar_posi[day,0],'+/-',quasar_posi[day,1])
    print('fluxes:',quasar_flux[day,0],'+/-',quasar_flux[day,1])
    print('background:',backgr_flux[day,0],'+/-',backgr_flux[day,1])
    print('chi-squared-reduced:',chi2red_val[day])
    ####

####
print('###')
print('### final steps')
print('###')

# important check preparation: position-deviation from initial guess:
position_deviation = np.sqrt((quasar_posi[:,0,0]-float(quasar_position.split(' ')[0]))**2+(quasar_posi[:,0,1]-float(quasar_position.split(' ')[1]))**2)
print('position deviations from initial input:')
print(position_deviation)

# scales and plotting preparation:
rel_ref_fluxes = ref_fluxes[0]/ref_fluxes
scales = np.median(rel_ref_fluxes,axis=1)
#scale_errors = np.maximum(abs(scales-np.quantile(rel_ref_fluxes,0.25,axis=1)),abs(np.quantile(rel_ref_fluxes,0.75,axis=1)-scales))
scale_errors = abs(np.quantile(rel_ref_fluxes,0.75,axis=1)-np.quantile(rel_ref_fluxes,0.25,axis=1))/2.0
#scale_errors = np.std(rel_ref_fluxes,ddof=1,axis=1)*np.sqrt(np.pi/(2.0*len(ref_positions)))
scale_errors[0] = np.median(scale_errors)
#scale_errors = np.zeros(N)
print('scales:')
print(scales,'+/-',scale_errors)

if four_images_with_galaxy == False:
    image_A = quasar_flux[:,0,0]*scales
    image_B = quasar_flux[:,0,1]*scales
    error_A = np.sqrt(quasar_flux[:,1,0]*quasar_flux[:,1,0]*scales*scales+quasar_flux[:,0,0]*quasar_flux[:,0,0]*scale_errors*scale_errors)
    error_B = np.sqrt(quasar_flux[:,1,1]*quasar_flux[:,1,1]*scales*scales+quasar_flux[:,0,1]*quasar_flux[:,0,1]*scale_errors*scale_errors)
elif four_images_with_galaxy == True:
    image_A = quasar_flux[:,0,0]*scales
    image_B = quasar_flux[:,0,1]*scales
    image_C = quasar_flux[:,0,2]*scales
    image_D = quasar_flux[:,0,3]*scales
    error_A = np.sqrt(quasar_flux[:,1,0]*quasar_flux[:,1,0]*scales*scales+quasar_flux[:,0,0]*quasar_flux[:,0,0]*scale_errors*scale_errors)
    error_B = np.sqrt(quasar_flux[:,1,1]*quasar_flux[:,1,1]*scales*scales+quasar_flux[:,0,1]*quasar_flux[:,0,1]*scale_errors*scale_errors)
    error_C = np.sqrt(quasar_flux[:,1,2]*quasar_flux[:,1,2]*scales*scales+quasar_flux[:,0,2]*quasar_flux[:,0,2]*scale_errors*scale_errors)
    error_D = np.sqrt(quasar_flux[:,1,3]*quasar_flux[:,1,3]*scales*scales+quasar_flux[:,0,3]*quasar_flux[:,0,3]*scale_errors*scale_errors)
else:
    print('ERROR!!! four_images_with_galaxy-switch not set correctly!')
backflux = backgr_flux[:,0]*scales
backerror = np.sqrt(backgr_flux[:,1]*backgr_flux[:,1]*scales*scales+backgr_flux[:,0]*backgr_flux[:,0]*scale_errors*scale_errors)

# chi-reduced-squared and residual statistics:
print('chi-squared-reduced: minimum:',np.min(chi2red_val),', median:',np.median(chi2red_val),'and maximum:',np.max(chi2red_val))
print('mean central residuals: minimum:',round(np.min(residuals),3),', median:',round(np.median(residuals),3),'and maximum:',round(np.max(residuals),3))

# text-file with all (relative to first of day) scales, background stds and weights:
txt = open("galfit-output-infos.txt", "w")
if four_images_with_galaxy == False:
    txt.write("file"+"\t"+"time"+"\t"+
              "x-pos"+"\t"+"x-pos-err"+"\t"+"y-pos"+"\t"+"y-pos-err"+"\t"+"pos-dev"+"\t"+
              "A-flux"+"\t"+"A-err"+"\t"+"B-flux"+"\t"+"B-err"+"\t"+
              "scales"+"\t"+"scale_err"+"\t"+"sky"+"\t"+"sky-err"+"\t"+"residuals"+"\t"+"chi2red"+"\n")
    txt.write("\n")
    for i in range(N):
        txt.write(str(file_list[i])+"\t"+str(round(time[i],4))+"\t"+
                  str(quasar_posi[i,0,0])+"\t"+str(quasar_posi[i,1,0])+"\t"+str(quasar_posi[i,0,1])+"\t"+str(quasar_posi[i,1,1])+"\t"+str(round(position_deviation[i],4))+"\t"+
                  str(round(image_A[i],4))+"\t"+str(round(error_A[i],4))+"\t"+str(round(image_B[i],4))+"\t"+str(round(error_B[i],4))+"\t"+
                  str(round(scales[i],3))+"\t"+str(round(scale_errors[i],3))+"\t"+str(round(backflux[i],3))+"\t"+str(round(backerror[i],3))+"\t"+str(round(residuals[i],3))+"\t"+str(chi2red_val[i])+"\n")
else:
    txt.write("file"+"\t"+"time"+"\t"+
              "x-pos"+"\t"+"x-pos-err"+"\t"+"y-pos"+"\t"+"y-pos-err"+"\t"+"pos-dev"+"\t"+
              "A-flux"+"\t"+"A-err"+"\t"+"B-flux"+"\t"+"B-err"+"\t"+"C-flux"+"\t"+"C-err"+"\t"+"D-flux"+"\t"+"D-err"+"\t"+
              "scales"+"\t"+"scale_err"+"\t"+"sky"+"\t"+"sky-err"+"\t"+"residuals"+"\t"+"chi2red"+"\n")
    txt.write("\n")
    for i in range(N):
        txt.write(str(file_list[i])+"\t"+str(round(time[i],4))+"\t"+
                  str(quasar_posi[i,0,0])+"\t"+str(quasar_posi[i,1,0])+"\t"+str(quasar_posi[i,0,1])+"\t"+str(quasar_posi[i,1,1])+"\t"+str(round(position_deviation[i],4))+"\t"+
                  str(round(image_A[i],4))+"\t"+str(round(error_A[i],4))+"\t"+str(round(image_B[i],4))+"\t"+str(round(error_B[i],4))+"\t"+str(round(image_C[i],4))+"\t"+str(round(error_C[i],4))+"\t"+str(round(image_D[i],4))+"\t"+str(round(error_D[i],4))+"\t"+
                  str(round(scales[i],3))+"\t"+str(round(scale_errors[i],3))+"\t"+str(round(backflux[i],3))+"\t"+str(round(backerror[i],3))+"\t"+str(round(residuals[i],3))+"\t"+str(chi2red_val[i])+"\n")
txt.close()

if four_images_with_galaxy == True:
    #print(galaxy_results)
    txt2 = open("galaxy_results.txt", "w")
    txt2.write("file"+"\t"+"time"+"\t"+
               "quapos x"+"\t"+"qua-err x"+"\t"+"quapos y"+"\t"+"qua-err y"+"\t"+
               "sersic x"+"\t"+"ser-err x"+"\t"+"sersic y"+"\t"+"ser-err y"+"\t"+
               "expdisk x"+"\t"+"exp-err x"+"\t"+"expdisk y"+"\t"+"exp-err y"+"\t"+
               "sersic galmag"+"\t"+"ser-err galmag"+"\t"+"expdisk galmag"+"\t"+"exp-err galmag"+"\t"+              
               "scales"+"\t"+"scale_err"+"\n")
    txt2.write("\n")
    for i in range(N):
        txt2.write(str(file_list[i])+"\t"+str(round(time[i],5))+"\t"+
                  str(quasar_posi[i,0,0])+"\t"+str(quasar_posi[i,1,0])+"\t"+str(quasar_posi[i,0,1])+"\t"+str(quasar_posi[i,1,1])+"\t"+
                  str(galaxy_results[i,0,0])+"\t"+str(galaxy_results[i,1,0])+"\t"+str(galaxy_results[i,0,1])+"\t"+str(galaxy_results[i,1,1])+"\t"+
                  str(galaxy_results[i,0,3])+"\t"+str(galaxy_results[i,1,3])+"\t"+str(galaxy_results[i,0,4])+"\t"+str(galaxy_results[i,1,4])+"\t"+
                  str(galaxy_results[i,0,2])+"\t"+str(galaxy_results[i,1,2])+"\t"+str(galaxy_results[i,0,5])+"\t"+str(galaxy_results[i,1,5])+"\t"+
                  str(round(scales[i],5))+"\t"+str(round(scale_errors[i],5))+"\n")
    txt2.close()

if refstar_position_output == True:
    for j in range(len(ref_positions)):
        if j < 9:
            refstar = '0'+str(int(j+1))
        else:
            refstar = str(int(j+1))
        txt3_name = "refstar"+refstar+"positions.txt"
        txt3 = open(txt3_name, "w")
        txt3.write("time"+"\t"+"\t"+"refxpos"+"\t"+"xposerr"+"\t"+"refypos"+"\t"+"yposerr"+"\n")
        for i in range(N):
            txt3.write(str(round(time[i],5))+
                       "\t"+str(round(refstarpos[i,j,0],5))+"\t"+str(round(refstarpos[i,j,1],5))+
                       "\t"+str(round(refstarpos[i,j,2],5))+"\t"+str(round(refstarpos[i,j,3],5))+"\n")
        txt3.close()

# check all results and mask potentially bad ones regarding position-deviation, chi2red, residual and background:
exclude_counter = 0
exclude_mask = np.zeros(N)
for i in range(N):
    if ( position_deviation[i] > 5.0 or chi2red_val[i] > 10.0 ) and ( residuals[i] > 40.0 or abs(backflux[i]) > 50.0 ) :
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
if four_images_with_galaxy == True:
    image_C = np.ma.array(image_C,mask=exclude_mask)
    image_D = np.ma.array(image_D,mask=exclude_mask)
    error_C = np.ma.array(error_C,mask=exclude_mask)
    error_D = np.ma.array(error_D,mask=exclude_mask)
backflux = np.ma.array(backflux,mask=exclude_mask)
backerror = np.ma.array(backerror,mask=exclude_mask)
scales = np.ma.array(scales,mask=exclude_mask)
scale_errors =  np.ma.array(scale_errors,mask=exclude_mask)
chi2red_val = np.ma.array(chi2red_val,mask=exclude_mask)
residuals = np.ma.array(residuals,mask=exclude_mask)
residual_errors = np.ma.array(residual_errors,mask=exclude_mask)
position_deviation = np.ma.array(position_deviation,mask=exclude_mask)
####

####
# plotting:
    
# lightcurves plot with linear y-axis:
plt.errorbar(time,image_A,yerr=error_A,fmt='o',markersize=2,elinewidth=1,capsize=3,label='image A')
plt.errorbar(time,image_B,yerr=error_B,fmt='o',markersize=2,elinewidth=1,capsize=3,label='image B')
if four_images_with_galaxy == True:
    plt.errorbar(time,image_C,yerr=error_C,fmt='o',markersize=2,elinewidth=1,capsize=3,label='image C')
    plt.errorbar(time,image_D,yerr=error_D,fmt='o',markersize=2,elinewidth=1,capsize=3,label='image D')   
plt.title('lightcurves of quasar images')
plt.legend()
plt.savefig('lightcurves_lin.jpg',dpi=800)
if plotting:
    plt.show()
plt.clf()
# lightcurves plot with logarithmic y-axis:
plt.errorbar(time,image_A,yerr=error_A,fmt='o',markersize=2,elinewidth=1,capsize=3,label='image A')
plt.errorbar(time,image_B,yerr=error_B,fmt='o',markersize=2,elinewidth=1,capsize=3,label='image B')
if four_images_with_galaxy == True:
    plt.errorbar(time,image_C,yerr=error_C,fmt='o',markersize=2,elinewidth=1,capsize=3,label='image C')
    plt.errorbar(time,image_D,yerr=error_D,fmt='o',markersize=2,elinewidth=1,capsize=3,label='image D')
plt.yscale('log')
plt.legend()
plt.title('lightcurves of quasar images')
plt.savefig('lightcurves_log.jpg',dpi=800)
if plotting:
    plt.show()
plt.clf()
# background, scale, chi-reduced-squared, residual and position-deviation plot:
plt.errorbar(time,backflux,yerr=backerror,fmt='o',markersize=2,elinewidth=1,capsize=3,label='background flux')
plt.errorbar(time,scales,yerr=scale_errors,fmt='o',markersize=2,elinewidth=1,capsize=3,label='image scale')
plt.errorbar(time,residuals,yerr=residual_errors,fmt='o',markersize=2,elinewidth=1,capsize=3,label='mean central residual (fit vs. data)')
plt.scatter(time,position_deviation,marker='x',color='tab:red',label='position deviation (fit vs. intitial)')
plt.scatter(time,chi2red_val-1.0,marker='o',s=4,color='tab:purple',label='$\chi_{red}^2-1$')
plt.legend(fontsize=7)
plt.title('additional info for lightcurves')
plt.savefig('lightcurves_supplement.jpg',dpi=800)
if plotting:
    plt.show()
plt.clf()
####

print('galfitting_betterinitialpos.py is finished!')
