import os
import sys
import numpy as np
from astropy.time import Time

#### 
#### INPUT THE DATE of the astrometry_ref.fits file as string of 8 numbers, e.g.: '20140526' for q2237!!!!
##########################
ref_date = '20140526' ######## q2237: '20140526', he2149: '20211005,'
##########################
#### 

######################################################################################################
#### DO NOT CHANGE THIS FILE FROM HERE ON - IT WILL BE ACCESSED BY aligning.py VIA interp.csh!!!! ####
######################################################################################################

# switch set by aligning.py: decides whether this correction will be applied
txt = open("gaiastarpositioncorrection_bool.txt", "r")
gaiastarpositioncorrection_bool = txt.read()
txt.close()

#print(sys.argv) # gaiaposcorrectionforaligning.py arguments

# check if this program got zero arguments:
# --> then refalignimage.data stars will be matched with gaia stars and pm's for good.data
if len(sys.argv) == 1: # sys.argv counts pythonfilename also as argument
    if gaiastarpositioncorrection_bool == 'True':
        # print statements can be read in alignpy_consollog.txt afterwards
        print('#########################################################')
        print('# reference and gaia star and pm matching for good.data #')
        print('#########################################################')
        
        # star matching:
        
        # load isis and gaia values:
        refaligndataRA,refaligndataDEC = np.loadtxt('refalignimage.data',usecols=(0,1),unpack=True)
        #print(refaligndataRA,refaligndataDEC)
        gaiaRA,gaiaDEC,pmRA,pmDEC,pmSN = np.loadtxt('gaia_pixelpositionandpm_list.txt',skiprows=1,usecols=(0,1,2,3,4),unpack=True)
        #print(gaiaRA,gaiaDEC,pmRA,pmDEC,pmSN)
        
        # open refalignimage.data to read it and 
        # open goodpm.data to write into it, i.e. append lines with pm:
        refalignimagedata = open('refalignimage.data','r')
        datalines = refalignimagedata.readlines()
        refalignimagedata.close()
        gooddata = open('goodpm.data','w')
        
        # match stars:
        match_counter = 0
        max_distance = 15.0 # distance in pixels must be below this to have a match
        min_pmSN = 5.0 # pmSN must be above this, to correct a star for pm (is already at at least 5 from the gaialist)
        #print(len(refaligndataRA),len(gaiaRA))
        for i in range(len(refaligndataRA)):
            # starpmRA/DEC_forgooddata will be written to good.data and:
            # if no star match: pmRA = 0.0 and pmDEC = 0.0 will be added to good.data of this star.
            starpmRA_forgooddata = 0.0
            starpmDEC_forgooddata = 0.0
            distance = np.sqrt((refaligndataRA[i]-gaiaRA)**2+(refaligndataDEC[i]-gaiaDEC)**2)
            mindex = np.argmin(distance)
            if distance[mindex] < max_distance and pmSN[mindex] > min_pmSN:                
                #print(distance[mindex],refaligndataRA[i],gaiaRA[mindex],refaligndataDEC[i],gaiaDEC[mindex])
                match_counter = match_counter + 1
                starpmRA_forgooddata = pmRA[mindex]
                starpmDEC_forgooddata = pmDEC[mindex]
            #print(starpmRA_forgooddata,starpmDEC_forgooddata)
            # write pm into the end of the line of the .data-file if it does not have a pm in there yet:
            #print(datalines[i])
            if str(datalines[i])[-4:-1] == 'nan':
                datalines[i] = datalines[i].strip('\n')+' '+str(starpmRA_forgooddata)+' '+str(starpmDEC_forgooddata)+' \n'
            #print(datalines[i])
            gooddata.writelines(datalines[i])
            
        gooddata.close()
            
        print('matched',match_counter,'of',len(refaligndataRA),'stars with a distance smaller than',max_distance,'pixels,')
        print('a pmSN of above',min_pmSN,'and written their pmRA, pmDEC and good.data to goodpm.data')
        
    else:
        os.system('cp refalignimage.data good.data')
        # nothing changed --> registering as usual:
        # the ref star positions are unchanged and copied to good data, where they are used by isis as usual


# check if this program got one argument:
# --> then refalignimage.data star positions will be matched with the gaia pm's for good.data
elif len(sys.argv) == 2:
    date = str(sys.argv[1])[24:32]

    if gaiastarpositioncorrection_bool == 'True':
        print('#########################################################')
        print('# daily star data correction with gaia pm for good.data #')
        print('#########################################################')
    
        #### initial test: just remove some stars from the list (the test is the third line with ####)
        #### do not activate anymore, as the first 100 stars are the 'best' ones!!!
        #### os.system('tail -n +101 <refalignimage.data >good.data') # remove the first 100 stars from good.data
        #### end of test
        
        # time of daily image:
        year = str(date[0:4])
        month = str(date[4:6])
        day = str(date[6:8])
        print('good.data positions will be corrected for:',day+'.'+month+'.'+year)
        utc_string = year+'-'+month+'-'+day+'T00:00:00.000' # UTC 00:00:00.000am on day-month-year
        julian_date = Time(utc_string,format='isot',scale='utc').jd # julian date
        #print(utc_string,julian_date)
        
        # time of reference image:
        ref_date_utc = str(ref_date[0:4])+'-'+str(ref_date[4:6])+'-'+str(ref_date[6:8])+'T00:00:00.000'
        julian_ref_date = Time(ref_date_utc,format='isot',scale='utc').jd
        #print(ref_date_utc,julian_ref_date)
        
        # time differnce to refernce image in days:
        time_difference = julian_date - julian_ref_date # time difference in days
        print('time difference to astrometry_ref.fits:',time_difference,'days')
    
        # position correction:
        initial_RA,initial_DEC,pm_RA_masperyear,pm_DEC_masperyear = np.loadtxt('goodpm.data',usecols=(0,1,5,6),unpack=True)
        pixel_scale = 0.387 # in arcsec/pixel
        pm_RA_pixelperday = pm_RA_masperyear/(365.25*1000*pixel_scale)
        pm_DEC_pixelperday = pm_DEC_masperyear/(365.25*1000*pixel_scale)
        RA_shift = pm_RA_pixelperday * time_difference
        DEC_shift = pm_DEC_pixelperday * time_difference
        print('maximum calculated absolute pixel shifts:',round(np.max(np.abs(RA_shift)),4),'and',round(np.max(np.abs(DEC_shift)),4))
        corrected_RA = initial_RA + RA_shift
        corrected_DEC = initial_DEC - DEC_shift
        
        # write good.data with corrected values for isis to align properly
        # open refalignimage.data to read it and 
        # open goodpm.data to write into it, i.e. append lines with pm:
        refalignimagedata_day = open('refalignimage.data','r')
        datalines_day = refalignimagedata_day.readlines()
        refalignimagedata_day.close()
        gooddata_day = open('good.data','w')
        if len(datalines_day) != len(corrected_RA):
            print('fatal error')
            sys.exit('number of lines to write and number of corrected lines are not equal')
        for i in range(len(datalines_day)):
            newgooddataline = ' '.join([str(corrected_RA[i]),str(corrected_DEC[i])]+datalines_day[i].split()[2:])+'\n'
            gooddata_day.writelines(newgooddataline)
        gooddata_day.close()
        
        print('#########################################################') # can be read in alignpy_consollog.txt afterwards
        print('# daily star data correction of good.data has finished! #')
        print('#########################################################')
        
    else:
        # the daily correction is not happening and isis registers unchanged as before without any correction
        print('#########################################################')
        print('# NO daily star data correction and good.data unchanged #')
        print('#########################################################')
        
    

# if refalignimage.data did not get zero or one argument interupt everything as it must get zero or one argument!
else:
    print('argument error: wrong number of arguments given to gaiaposcorrectionforaligning.py!')
    sys.exit('argument error: wrong number of arguments given to gaiaposcorrectionforaligning.py!')    
