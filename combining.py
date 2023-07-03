import numpy as np
import matplotlib.pyplot as plt

import os
from pathlib import Path

from astropy.io import fits
from astropy.stats import sigma_clipped_stats, gaussian_sigma_to_fwhm
from photutils.aperture import CircularAperture, aperture_photometry

from photutils.centroids import gaussian1d_moments
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.modeling.models import Const1D, Gaussian1D, Moffat1D
from astropy.utils.exceptions import AstropyUserWarning
import warnings

print('starting combine.py!')

local = (str(os.path.abspath(os.getcwd())))
local_path = Path(local)
data_path = Path(local+os.sep+'aligned')
combined_path = Path(local+os.sep+'combined')
combi_error_path = Path(local+os.sep+'combi_error')
print('local path:      ',local_path)
print('data path:       ',data_path)
print('combined path:   ',combined_path)
print('combi-error path:',combi_error_path)
fits_files = list(data_path.glob('a_t_*.fits'))
print('number of total fits-files:',len(fits_files))

N_pix = 3700 # number of pixels in both directions of the data-files ('a_t_*.fits')
             # i.e. the shape of the data in the input and output files is always N_pix*N_pix!
             # Must be N_pix = 3700 for a_t_*.fits-files of he2149, q2237, ...!

dates = [] # dates, aber inklusive Telsekopdetails im string davor!
for x in fits_files:
    y = str(x)
    z = y.replace(str(data_path),'') # remove data path
    fitsname = z[1:]                 # remove "/" infront of a_t_...
    date = fitsname[4:26]            # remove "a_t_" at beginning and stop with date
    if not (date in dates):
        dates.append(date)
print('all dates with telescope details:',dates)
print('total number of dates i.e. number of output-c_*.fits-files:',len(dates))

# .txt-file for scales, weights and seeing per image:
txt = open("imageinfos_combined.txt", "w")
txt.write("filename"+"\t\t\t\t\t"+"scale"+"\t"+"sca-err"+"\t"+"weight"+"\t"+
          "wei_err"+"\t"+"backgr"+"\t"+"bac-err"+"\t"+"seeing"+"\t"+"see-err"+"\n")

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

#print(initial_positions,'initial positions of ref-stars for scales, weights and seeing')

iternumber = 3 # number of iterations of centroiding and FWHM-estimating

plotting = False # if True, plots of scale-stars, their apertures and the final combined images are shown

############################################
### switch for using a min-max-rejection ###
############################################
use_MinMaxRej = True
############################################
    
###############################################################################
###############################################################################
###############################################################################
###############################################################################

# modified centroid function for ref-star positions and FWHM:

# Source code for photutils.centroids.gaussian to be found under
# https://photutils.readthedocs.io/en/stable/_modules/photutils/centroids/gaussian.html
# Licensed under a 3-clause BSD style license - see LICENSE.rst
# Here only centroid_1dg, but modified and called centroid_1dgm
# Original function: from photutils.centroids import centroid_1dg 
# Modified by Christian Sorgenfrei in September 2021

def centroid_1dgm(data, error=None, mask=None):
    
    #########################################################################
    # Calculate the centroid of a 2D array by fitting 1D Gaussians to the   #
    # marginal ``x`` and ``y`` distributions of the array.                  #
    #                                                                       #
    # Non-finite values (e.g., NaN or inf) in the ``data`` or ``error``     #
    # arrays are automatically masked. These masks are combined.            #
    #                                                                       #
    # Parameters                                                            #
    # ----------                                                            #
    # data : array_like                                                     #
    #     The 2D data array.                                                #
    #                                                                       #
    # error : array_like, optional                                          #
    #     The 2D array of the 1-sigma errors of the input ``data``.         #
    #                                                                       #
    # mask : array_like (bool), optional                                    #
    #     A boolean mask, with the same shape as ``data``, where a `True`   #
    #     value indicates the corresponding element of ``data`` is masked.  #
    #                                                                       #
    # Returns                                                               #
    # -------                                                               #
    # centroid : `~numpy.ndarray`                                           #        
    #     The ``x, y`` coordinates of the centroid.                         #                    
    #                                                                       #
    # and NEW : arrays with FWHM estimates around the centroids             #
    #           also for 1D-Moffat-fits, not only 1D-Gaussian-fits          #
    #########################################################################
    
    data = np.ma.asanyarray(data)
    
    if mask is not None and mask is not np.ma.nomask:
        mask = np.asanyarray(mask)
        if data.shape != mask.shape:
            raise ValueError('data and mask must have the same shape.')
        data.mask |= mask

    if np.any(~np.isfinite(data)):
        data = np.ma.masked_invalid(data)
        warnings.warn('Input data contains non-finite values (e.g., NaN or '
                      'inf) that were automatically masked.',
                      AstropyUserWarning)

    if error is not None:
        error = np.ma.masked_invalid(error)
        if data.shape != error.shape:
            raise ValueError('data and error must have the same shape.')
        data.mask |= error.mask
        error.mask = data.mask

        xy_error = [np.sqrt(np.ma.sum(error**2, axis=i)) for i in (0, 1)]
        xy_weights = [(1.0 / xy_error[i].clip(min=1.e-30)) for i in (0, 1)]
    else:
        xy_weights = [np.ones(data.shape[i]) for i in (1, 0)]

    # assign zero weight where an entire row or column is masked
    if np.any(data.mask):
        bad_idx = [np.all(data.mask, axis=i) for i in (0, 1)]
        for i in (0, 1):
            xy_weights[i][bad_idx[i]] = 0.

    xy_data = [np.ma.sum(data, axis=i).data for i in (0, 1)]

    constant_init = np.ma.min(data)
    g_centroid = [] # gaussian centroid positions
    g_cenwidth = [] # gaussian fwhm of centroid
    m_centroid = [] # moffat centroid positions
    m_cenwidth = [] # moffat fwhm of centroid
    for (data_i, weights_i) in zip(xy_data, xy_weights):
        params_init = gaussian1d_moments(data_i)
        #print(params_init)
        g_init = Const1D(constant_init) + Gaussian1D(*params_init)
        fitter = LevMarLSQFitter()
        x = np.arange(data_i.size)
        g_fit = fitter(g_init, x, data_i, weights=weights_i)
        #print(g_fit)
        
        g_centroid.append(g_fit.mean_1.value)
        g_FWHM = gaussian_sigma_to_fwhm*g_fit.stddev_1.value
        g_cenwidth.append(g_FWHM)
        
        # now with moffat, not gaussian:
            
        moffat_alpha = 2.0, # include alpha for moffat in params_init tuple
        params_init = params_init + moffat_alpha # use amplitude, mean and gamma~sigma from above
        #print(params_init)
        m_init = Const1D(constant_init) + Moffat1D(*params_init)
        #m_init.alpha_1.fixed = True # to have the possibility to keep alpha fixed
        m_fit = fitter(m_init, x, data_i, weights=weights_i)
        #print(m_fit)
        
        m_centroid.append(m_fit.x_0_1.value)
        gamma = abs(m_fit.gamma_1.value)
        alpha = abs(m_fit.alpha_1.value)
        #print('moffat gamma and alpha:',gamma,alpha)
        m_FWHM = 2.0*gamma*np.sqrt(2.0**(1.0/alpha)-1)
        m_cenwidth.append(m_FWHM)

    return np.array(g_centroid), np.array(g_cenwidth), np.array(m_centroid), np.array(m_cenwidth)

###############################################################################
###############################################################################
###############################################################################
###############################################################################

# special warning counters and warn-day-lists:
zero_file_count = 0
zero_file_list = []
type_warn_count = 0
type_warn_list = []
data_shape_count = 0
data_shape_list = []
negative_data_count = 0
negative_data_list = []
back_medi_count = 0
back_medi_list = []
fit_error_count = 0
fit_error_list = []
centroid_warn_count = 0
centroid_warn_list = []
aperture_warn_count = 0
aperture_warn_list = []
fatal_scale_count = 0
fatal_scale_list = []
scale_warn_count = 0
scale_warn_list = []
one_file_count = 0
one_file_list = []
errcorr_factor_count = 0
errcorr_factor_list = []
date_error_count = 0
date_error_list = []

#### beginning of looping through all dates:

day_counter = 0
for day in dates:
    day_counter = day_counter + 1
    print('###########')
    print('###########')
    print('########### combine images of date',day_counter,'of',len(dates))
    print('###########')
    print('###########')
    
    #### finding the current fits-files:
    current_fits_files = []
    for file in fits_files:
        if day in str(file):
            current_fits_files.append(file)
    N_current_fits_files = len(current_fits_files)
    print('number of combined images (at',day + '):',N_current_fits_files)
    #print('list of current fits files:',current_fits_files)
    if N_current_fits_files < 1:
        # cannot happen!
        zero_file_count = zero_file_count + 1
        zero_file_list.append(day)
        print('#### warning: number of files under 1; impossible error!')
        continue # continue with next day
    
    same_type_check = []
    for file in current_fits_files:
        #print(file)
        same_type_check.append(str(file)[-40:-14]+str(file)[-9:])
    #print(same_type_check)
    if (all(x == same_type_check[0] for x in same_type_check)) != True:
        # should not be possible to happen anymore as we include the tel.-typ in the dates!
        type_warn_count = type_warn_count + 1
        type_warn_list.append(day)
        print('#### warning: multiple instruments combined at date',day)       

    #### combining starts:    
    print('start combining!')
    
    # arrays for aperture phot. stars and fwhms:
    star_flux = np.zeros((N_current_fits_files,len(initial_positions))) # needed for calculating the scales
    star_fwhm = np.zeros((N_current_fits_files,len(initial_positions))) # needed for estimating the seeing
    
    # arrays for background estimates and their errors:
    backgr_est = np.zeros(N_current_fits_files)
    backgr_est_error = np.zeros(N_current_fits_files)

    # initialising the data_cube and the error_cube:
    data_cube = np.zeros((N_current_fits_files,N_pix,N_pix))
    error_cube = np.zeros((N_current_fits_files,N_pix,N_pix))
    
    # counter for massive scaling problems:
    fatal_scaling_error_counter = 0
    
    # extracting the data from the files, preparing the image-data-cube and ref-star centroiding, FWHMs and aperture photometry:
    current_counter = 0
    for file in current_fits_files:
        #print(file)
        data = fits.open(file)
        if data[0].data.shape != (N_pix,N_pix):
            data_shape_count = data_shape_count + 1
            data_shape_list.append(str(file)[-40:])
            print('#### warning: file has wrong shape!')
            break # interrupt for this day due to wrong data shape - really should not happen!!!
        
        trim = 100 # just for sigma-clipping-background-estimation
        mean,median,std = sigma_clipped_stats(data[0].data[trim:-trim,trim:-trim],sigma = 3.0,std_ddof=1)
        background = median
        error_background = std # std i.e. typical fluctuation size of the individual background pixels
        print('background =',background,'+/-',error_background)
        image = data[0].data - background
        
        # save background estimates:
        backgr_est[current_counter] = background
        backgr_est_error[current_counter] = error_background
        
        # error image preparation and calculation:
        
        if len(data[0].data[data[0].data<0]) != 0:
            negative_data_count = negative_data_count + 1
            negative_data_list.append(str(file)[-40:])
            print('#### warning: negative data values in file!')
            
        N_back_pix = np.size(image[np.absolute(image) < 3.0*error_background]) # number of not clipped (i.e. background) pixels
        err_median_background = np.sqrt(np.pi/2.0)*error_background/np.sqrt(N_back_pix) # error of median of background, that got subtracted
        if abs(err_median_background) > 1.0:
            back_medi_count = back_medi_count + 1
            back_medi_list.append(str(file)[-40:])
            print('#### warning: unlikely large background median error of:',err_median_background,'--> set to zero!')
            err_median_background = 0.0
        
        image_error = np.sqrt(data[0].data + err_median_background**2) # first term: poisson error of data; second term: error of median of background
        
        ####################################################################
        # save data (i.e. image and image_error) in data-cubes:
        data_cube[current_counter] = image
        error_cube[current_counter] = image_error
        ####################################################################
        
        # centroids for apertures and fwhm for seeing:
        positions = np.zeros(np.shape(initial_positions))
        positions = positions + initial_positions # always (for all files) use initial_positions as first position guesses
        #print(positions)
        cuthalfsize = 20 # +/- around position for first centroids and FWHM iteration
        
        centroid_warnings = 0 # counts the number of warnings from the ref-stars per file
        
        for i in range(len(positions)):
            # for all 'i' ref-stars:
            for iteration in range(iternumber):
                #print(cuthalfsize,'is the half edge length of iteration',iteration)
            
                edge1 = np.floor(positions[i,1]-cuthalfsize).astype(int)
                edge2 = np.ceil(positions[i,1]+cuthalfsize).astype(int)
                edge3 = np.floor(positions[i,0]-cuthalfsize).astype(int)
                edge4 = np.ceil(positions[i,0]+cuthalfsize).astype(int)
                #print(edge1,edge2,edge3,edge4)
                
                try:
                    g_pos,g_fwhm,m_pos,m_fwhm = centroid_1dgm(image[edge1:edge2,edge3:edge4],
                                                          error=image_error[edge1:edge2,edge3:edge4])
                except:
                    print('WARNING: centroid_1dgm failed --> skipping this iteration')
                    fatal_scaling_error_counter = fatal_scaling_error_counter + 1
                    continue
                g_pos = g_pos + np.array([edge3,edge1])
                m_pos = m_pos + np.array([edge3,edge1])
                #print(g_pos,m_pos)
                
                # geometric means:
                g_fwhm_geomean = np.sqrt(g_fwhm[0]*g_fwhm[1])
                m_fwhm_geomean = np.sqrt(m_fwhm[0]*m_fwhm[1])
                #print(g_fwhm_geomean,m_fwhm_geomean)
                
                # use moffat results for new positions and edge cutting:
                positions[i,0] = m_pos[0]
                positions[i,1] = m_pos[1]
                cuthalfsize = 3.0*m_fwhm_geomean # new half edge is 3.0*fwhm
                
                # check after last iteration whether fit values seem plausible and if gaussian and moffat fit got similar results:
                if iteration == iternumber-1:                 
                    #print(iteration)
                    
                    ini_pos_diff = np.absolute(positions[i]-initial_positions[i])
                    if np.any(g_fwhm > 20.0) or np.any(m_fwhm > 20.0) or np.any(ini_pos_diff > 10.0):
                        fit_error_count = fit_error_count + 1
                        fit_error_list.append(str(file)[-40:])
                        print('#### warning: implausible fit values for g/m-fwhm and centroid - initial position',g_fwhm,m_fwhm,ini_pos_diff)
                    
                    pos_dev = np.sqrt((g_pos[0]-m_pos[0])**2+(g_pos[1]-m_pos[1])**2)
                    fwhm_dev = abs(g_fwhm_geomean-m_fwhm_geomean)
                    if pos_dev > 1.0 and fwhm_dev > 2.0:
                        # warning if large position- and fwhm-deviations:
                        centroid_warnings = centroid_warnings + 1
                        if centroid_warnings == 2:
                            # only warn at the end if it happend at least twice in the file:
                            centroid_warn_count = centroid_warn_count + 1
                            centroid_warn_list.append(str(file)[-40:])
                        print('#### warning: potentially large gaussian and moffat fit deviation in ref star '+str(i)+':')
                        print('              position deviation is:',round(pos_dev,3),'and FWHM deviation is:',round(fwhm_dev,3))
    
            # use moffat results of final iteration as final positions and fwhm
            #print('final pos:',positions[i])
            #print('final fwhm:',m_fwhm_geomean)
            star_fwhm[current_counter,i] = m_fwhm_geomean
         
            # apertur photometrie for image scaling:
            try:
                aperture_radius = 1.0*star_fwhm[current_counter,i] # set circular aperture radius too once the estimated FWHM
                aperture = CircularAperture(positions[i],r=aperture_radius)
                #print(aperture)
                phot_table = aperture_photometry(image,aperture)
                star_flux[current_counter,i] = phot_table['aperture_sum'][0]
            except:
                print('WARNING: aperture flux failed')
                aperture_warn_count = aperture_warn_count + 1 
                aperture_warn_list.append(str(file)[-40:])
                if i < 0:
                    star_flux[current_counter,i] = np.mean(star_flux[current_counter,0:i-1])  
                else:
                    star_flux[current_counter,i] = np.nan
            
            if plotting:
                edge1 = np.floor(positions[i,1]-20).astype(int)
                edge2 = np.ceil(positions[i,1]+20).astype(int)
                edge3 = np.floor(positions[i,0]-20).astype(int)
                edge4 = np.ceil(positions[i,0]+20).astype(int)
                plt.imshow(np.log10(abs(image[edge1:edge2,edge3:edge4])+0.1),origin='lower')
                plt.title('log-image of scaling star '+str(i+1))
                #plt.imshow(image[edge1:edge2,edge3:edge4],origin='lower')
                #plt.title('lin-image of scaling star '+str(i+1))
                aperture.plot(color='r',origin=(edge3,edge1))
                plt.colorbar()
                plt.show()
                plt.clf()
                
        current_counter = current_counter + 1
        
    #scale calculation:
    #print(star_flux)
    rel_flux = star_flux[0]/star_flux # inverse relative fluxes
    #print(rel_flux)
    scales = np.nanmedian(rel_flux,axis=1)
    if N_current_fits_files == 1:
        scales_error = np.array([0.0]) # if only one picture scale/weight error is 0 not nan!
    else:
        scales_error = np.nanstd(rel_flux,axis=1,ddof=1)*np.sqrt(np.pi/(2.0*len(positions)))
        scales_error[0] = np.mean(scales_error[1:]) # error of the first scale should not be zero!
    
    # set scales to 1.0 if a major problem occured!
    for i in range(N_current_fits_files):
        if np.isnan(scales[i]) == True or np.isinf(scales[i]) == True or scales[i] >= 20.0 or scales[i] <= 0.05 or np.isnan(scales_error[i]) == True or np.isinf(scales_error[i]) == True:
            print('WARNING: set scales and errors to 1, because scale-determination failed!:',scales[i],scales_error[i])
            scales[i] = 1.0
            scales_error[i] = 1.0
            fatal_scaling_error_counter = fatal_scaling_error_counter + 1
                
    print('scales: ',np.around(scales,4),'+/-',np.around(scales_error,4))
    #print('scales:',scales,'+/-',scales_error,'i.e. rel error:',scales_error/scales*100.0,'%')
    
    if fatal_scaling_error_counter > 2:
        fatal_scale_count = fatal_scale_count + 1
        fatal_scale_list.append(day)
    
    # weights calculation:
    weights = 1.0/scales # comes from variance weighting
    weights_error = weights**2*scales_error
    print('weights:',np.around(weights,4),'+/-',np.around(weights_error,4))
    #print('weights:',weights,'+/-',weights_error,'i.e. rel error:',weights_error/weights*100.0,'%')
    
    if (np.any(abs(scales-1.0) >= 0.05)) or (np.any(abs(weights-1.0) >= 0.05)):
        scale_warn_count = scale_warn_count + 1
        scale_warn_list.append(day)
        print('#### warning: possible scaling/weighting-error!')
    
    #seeing calculation:
    #print(star_fwhm)
    seeing = np.median(star_fwhm,axis=1)
    seeing_error = np.std(star_fwhm,axis=1,ddof=1)*np.sqrt(np.pi/(2.0*len(positions)))
    print('seeing: ',np.around(seeing,4),'+/-',np.around(seeing_error,4),'pixels')
    #print('seeing:',seeing,'+/-',seeing_error,'pixels i.e. rel error:',seeing_error/seeing*100.0,'%')
    
    # background estimates:
    #print('background',backgr_est,'+/-',backgr_est_error,'pixels i.e. rel error:',backgr_est_error/backgr_est*100.0,'%')
    
    # text-file with all (relative to first of day) scales, background stds and weights:
    txt.write("\n") # one empty line between the days
    txt_count = 0
    for file in current_fits_files:
        txt.write(str(file)[-40:]+
                  "\t"+str(round(scales[txt_count],4))+"\t"+str(round(scales_error[txt_count],4))+
                  "\t"+str(round(weights[txt_count],4))+"\t"+str(round(weights_error[txt_count],4))+
                  "\t"+str(round(backgr_est[txt_count],1))+"\t"+str(round(backgr_est_error[txt_count],1))+
                  "\t"+str(round(seeing[txt_count],4))+"\t"+str(round(seeing_error[txt_count],4))+"\n")
        txt_count = txt_count + 1                          
        
    if N_current_fits_files == 1:
        one_file_count = one_file_count + 1
        one_file_list.append(day)
        print('#### warning: only one file to combine!')
    
    # number of symmetric N_rej:
    if use_MinMaxRej == True:
        if N_current_fits_files < 3:
            N_rej = 0
        elif N_current_fits_files >= 3 and N_current_fits_files <= 5:
            N_rej = 1
        else:
            N_rej = 2
    else:
        N_rej = 0 # ATTENTION: this manually overwrites the minmax-rejection, i.e. no minmaxrej!!!                                                             
    
    print('Will reject',2*N_rej,'of',N_current_fits_files,'image-values per pixel.')
    
    #### Now: weighted combining of the scaled images including the minmax-rejection via np.argsort:
    
    # image scaling:
    for i in range(N_current_fits_files):
        data_cube[i] = data_cube[i]*scales[i]
    
    # sorting (of scaled images, image errors, scales, scale errors weights and weight errors) via np.argsort:
    sorting = np.argsort(data_cube,axis=0)
    #print('sorting:',sorting,sorting.shape)
    data_cube = np.take_along_axis(data_cube,sorting,axis=0)
    error_cube = np.take_along_axis(error_cube,sorting,axis=0)
    #print(data_cube,data_cube.shape)
    scale_cube =  np.ones(np.shape(data_cube))
    scale_cube_error = np.ones(np.shape(data_cube))
    for i in range(N_current_fits_files):
        scale_cube[i] = scale_cube[i]*scales[i] 
        scale_cube_error[i] = scale_cube_error[i]*scales_error[i] 
    scale_cube = np.take_along_axis(scale_cube,sorting,axis=0)
    scale_cube_error = np.take_along_axis(scale_cube_error,sorting,axis=0)
    weight_cube = 1.0/scale_cube
    weigth_cube_error = scale_cube_error/scale_cube**2
    #print(scale_cube,scale_cube_error,weight_cube,weigth_cube_error)
    
    # minmax-rejection only for N_rej>0, 
    # or else one would slice arrays with [0:0] to empty arrays where no rejecting should have happend:
    if N_rej > 0:
        data_cube = data_cube[N_rej:-N_rej]
        error_cube = error_cube[N_rej:-N_rej]
        #scale_cube = scale_cube[N_rej:-N_rej] # not needed later
        #scale_cube_error = scale_cube_error[N_rej:-N_rej] # not needed later
        weight_cube = weight_cube[N_rej:-N_rej]
        weigth_cube_error = weigth_cube_error[N_rej:-N_rej]
    #print(data_cube.shape,error_cube.shape,weight_cube.shape,weigth_cube_error.shape)
    
    ##########
    
    #### actually combining:
    combined_data = np.average(data_cube,weights=weight_cube,axis=0)
    
    #### error image:
    combined_error = (1.0/np.sum(weight_cube,axis=0)) * np.sqrt(np.sum(error_cube**2 + combined_data**2*weigth_cube_error**2,axis=0))
    
    ##########
    
    # final back-stats and error-scaling-correction:
        
    final_combi_stat = sigma_clipped_stats(combined_data,sigma = 3.0,std_ddof=1)
    print('back-stats in combined data:',[round(stat_val,4) for stat_val in final_combi_stat])
    # first two values should be ~ 0 i.e. vanishing background in combined image
    # and the third is its variation (std of the background pixel; how the background varies)
    
    #quick_back_estim = np.sqrt(np.mean(backgr_est)/((N_current_fits_files-2.0*N_rej)*np.mean(weights)))
    #print('expected background error:',round(quick_back_estim,4))
    # is approx. equal to median of error pic, because that is an approximation of error formula for background
    # this should be equal to the std from before, but often isn't, because of minmaxrej, median worse than mean, small number of images, etc.
    
    prefinal_error_stat = sigma_clipped_stats(combined_error,sigma = 3.0,std_ddof=1)
    print('prefinal back-stats in combined error:',[round(stat_val,4) for stat_val in prefinal_error_stat])
    # second two should now approx. equal to std of combined background, because it is expectation value of background errors, and the third value is its variation
    
    # correction factor to correct the scale of errors:
    #error_correction_factor = final_combi_stat[2]/quick_back_estim + 0.001 # add to make the errors save and slight bigger
    #print('error-image-correction-factor:',round(error_correction_factor,4))
    
    # correction factor to correct the scale of errors:
    error_correction_factor = final_combi_stat[2]/prefinal_error_stat[1] # true combined std / 'wrong' median error
    print('error-image-correction-factor:',round(error_correction_factor,4))
    
    
    if error_correction_factor < 0.5 or error_correction_factor > 2.0:
        errcorr_factor_count =  errcorr_factor_count + 1
        errcorr_factor_list.append(day)
        print('#### warning: error-image-correction-factor scales errors more than 100%')
    
    # correcting the error image by including the factor, i.e. scaling it acording to the background variations:
    combined_error = error_correction_factor*combined_error
    
    #final_error_stat = sigma_clipped_stats(combined_error,sigma = 3.0,std_ddof=1)
    #print('back-stats in combined error:',[round(stat_val,4) for stat_val in final_error_stat])
    # second two should now approx. equal to std of combined background, because it is expectation value of background errors, and the third value is its variation
    
    ##########
    
    #### The combined data is in the array 'combined_data'
    #### and the error data of it is in 'combined_error'.
    
    if plotting:
        plt.imshow(np.log10(abs(combined_data[1550:1851,1550:1851])+1.0),origin='lower')
        #plt.imshow(np.log10(abs(combined_data)+1.0), origin='lower')
        plt.title('combined image no. '+str(day_counter))
        plt.colorbar()
        plt.show()
        plt.clf()
    
    print('Make two new fits containing combined_data and combined_error:')
    
    # make new fits file with combined data:
    data = fits.open(current_fits_files[0])
    data[0].data = combined_data # overwrite with combined_data and then soon write to new folder and file  
    u = str(current_fits_files[0])
    v = u.replace(str(data_path),'')
    fitsname_without_number = v[1:27]+v[32:]
    #print(fitsname_without_number)
    if not day in fitsname_without_number:
        #should not be possible to happen anymore!
        date_error_count = date_error_count + 1
        date_error_list.append(day)
        print('### warning: date error in',day,'###')
    
    # combined data output:
    print('combined output',day_counter,'of',len(dates),': ...'+os.sep+'combined'+os.sep+'c_'+fitsname_without_number)
    data.writeto(str(combined_path)+os.sep+'c_'+str(fitsname_without_number),overwrite=True)
    data.close()
    
    # same for error data:
    err_data = fits.open(current_fits_files[0])
    err_data[0].data = combined_error # overwrite with combined_error and then soon write to new folder and file  
    #print(fitsname_without_number)
    
    print('combi-error output',day_counter,'of',len(dates),': ...'+os.sep+'combi_error'+os.sep+'e_'+fitsname_without_number)
    err_data.writeto(str(combi_error_path)+os.sep+'e_'+str(fitsname_without_number),overwrite=True)
    err_data.close()

txt.close()

# end of program: show the warnings and then finish:
print('###########')
print('###########')
if type_warn_count != 0:
    print('WARNING!: multiple instruments combined at',type_warn_count,'dates!')
    print(type_warn_list)
if data_shape_count != 0:
    print('WARNING!: wrong data shape in',data_shape_count,'files!')
    print(data_shape_list)
if negative_data_count != 0:
    print('WARNING!: negative data found in',negative_data_count,'files!')
    print(negative_data_list)
if back_medi_count != 0:
    print('WARNING!: unlikely large background median error in',back_medi_count,'files!')
    print(back_medi_list)
if fit_error_count != 0:
    print('WARNING!: implausible final fit values encountered',fit_error_count,'times!')
    compact_fit_error_list = list(set([x[4:26] for x in fit_error_list]))
    print(compact_fit_error_list)
else:
    compact_fit_error_list = []
if centroid_warn_count != 0:
    print('WARNING!: potentially erroneous centroiding in',centroid_warn_count,'files!')
    print(centroid_warn_list)
if aperture_warn_count != 0:
    print('WARNING!: aperture flux estimations failed',aperture_warn_count,'times!')
    print(aperture_warn_list)
if fatal_scale_count != 0:
    print('WARNING!: potentially fatal scaling problems at',fatal_scale_count,'dates!')
    print(fatal_scale_list)
if scale_warn_count != 0:
    print('WARNING!: potentially erroneous scale/weight at',scale_warn_count,'dates!')
    print(scale_warn_list)
if one_file_count != 0:
    print('WARNING!: only one image per day in',one_file_count,'dates!')
    print(one_file_list)
if zero_file_count != 0:
    print('WARNING!: impossibly less then one image per day in',zero_file_count,'dates!')
    print(zero_file_list)
if errcorr_factor_count != 0:
    print('WARNING!: large error changes through correction factor in',errcorr_factor_count,'dates!')
    print(errcorr_factor_list)
if date_error_count != 0:
    print('WARNING!: wrong date in',date_error_count,'dates!')
    print(date_error_list)
print('###########')
print('###########')
alldatelist = type_warn_list+fatal_scale_list+scale_warn_list+one_file_list+zero_file_list+errcorr_factor_list+date_error_list+compact_fit_error_list
allfilelist = data_shape_list+negative_data_list+back_medi_list+centroid_warn_list+aperture_warn_list
thelist = alldatelist + [x[4:26] for x in allfilelist]
seen1 = set()
dupl1 = [x for x in thelist if x in seen1 or seen1.add(x)]
seen2 = set()
dupl2 = [x for x in dupl1 if x in seen2 or seen2.add(x)]
seen3 = set()
dupl3 = [x for x in dupl2 if x in seen3 or seen3.add(x)]
dupllist_dateatleast4times = list(set(dupl3))
print('WARNING!: potentially bad combined images with at least 4 problems at',len(dupllist_dateatleast4times),'dates:')
print(dupllist_dateatleast4times)
print('###########')
print('###########')
print('WARNING!: potentially bad combined images with at least 4 problems to open in ds9:')
print('ds9 *'+'* *'.join(list(set([x[-8:] for x in dupllist_dateatleast4times])))+'*')
print('combining.py is finished!')