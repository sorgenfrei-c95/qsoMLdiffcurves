import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import iqr

# Which quasar is worked on? --> Set that one to True, all others to false
he1104 = False
he2149 = False
q2237 = True

# masking switch for better median 'true' position for q2237 because of galaxy model problems especially in V discussed with Robert on Dec. 21, 2022:
q2237_masking = True # only effects which data points are used for median determination with q2237 data!

# should be False for console use:
show_plots = False 

# loading quasarpositions.txt data:
file = np.genfromtxt('quasarpositions.txt',skip_header=2,usecols=(0),dtype='str') 
time,quapos_x,qua_err_x,quapos_y,qua_err_y=np.loadtxt('quasarpositions.txt',skiprows=2,usecols=(1,2,3,4,5),unpack=True) 

# color map for plotting:
cmap = plt.get_cmap("tab10")

# x position shift for nicer plotting:
x_shift = np.floor(np.median(quapos_y)-np.median(quapos_x))+3.0
#print(x_shift)
if x_shift > 0.0:
    x_shift_sign = '$+$'
else:
    x_shift_sign = '$-$'

# linear function:
def lin(t,a,b): 
    return a*(t-2014)+b

# linear fit to data for first and third plot:
x_popt,x_pcov = curve_fit(lin,time,quapos_x+x_shift,p0=(0.0,0.0),sigma=qua_err_x) 
x_a = x_popt[0]
x_a_err = np.sqrt(x_pcov[0,0])
x_b = x_popt[1]
y_popt,y_pcov = curve_fit(lin,time,quapos_y,p0=(0.0,0.0),sigma=qua_err_y) 
y_a = y_popt[0]
y_a_err = np.sqrt(y_pcov[0,0])
y_b = y_popt[1]

# load gaia_pixelpositionandpm_list.txt to extract psf star data by comparing the pixel position of the psf star:
PSF_RA,PSF_DEC,PSF_RApm,PSF_DECpm,PSF_SN=np.loadtxt('../gaia_pixelpositionandpm_list.txt',skiprows=1,usecols=(0,1,2,3,4),unpack=True) # RA[pix],DEC[pix],RApm[mas/year],DECpm[mas/year],SN[totalpm/error]

# psf star position box (DEC;RA not the other way around, because of ds9) from galfitting.py:
quasar_counter = 0
if he1104 == True:
    psfpixpos = [1770,1800,1700,1730] # he1104
    quasar_counter = quasar_counter + 1
elif he2149 == True:
    psfpixpos = [1703,1733,1786,1816] # he2149
    quasar_counter = quasar_counter + 1
elif q2237 == True:
    psfpixpos = [2120,2150,2213,2243] # q2237
    quasar_counter = quasar_counter + 1
else:
    sys.exit('NO QUASAR SELECTED!')
if quasar_counter != 1:
    sys.exit('NOT EXACTLY ONE QUASAR SELECTED!')

# find psf star index in gaia_pixelpositionandpm_list.txt:
post_register_trimming = 100 # after registering with ISIS for which the GAIA table was written the images where trimmed by 100 pixels!
PSF_index_x = np.where((PSF_RA-post_register_trimming > psfpixpos[2]) & (PSF_RA-post_register_trimming < psfpixpos[3]))
PSF_index_y = np.where((PSF_DEC-post_register_trimming > psfpixpos[0]) & (PSF_DEC-post_register_trimming < psfpixpos[1]))
#print(PSF_index_x,PSF_index_y,np.intersect1d(PSF_index_x,PSF_index_y))
if len(np.intersect1d(PSF_index_x,PSF_index_y)) == 1:
    PSF_index = np.intersect1d(PSF_index_x,PSF_index_y)
else:
    sys.exit('PSF STAR GAIA VALUES NOT FOUND!')
#print(PSF_index,PSF_RA[PSF_index],PSF_DEC[PSF_index],PSF_RApm[PSF_index],PSF_DECpm[PSF_index],PSF_SN[PSF_index])

# psf proper motion from GAIA for betterpositions and from second plot on:
psf_RApm = PSF_RApm[PSF_index][0]   # in -mas/year from GAIA        (example q2237 approx. -24.485)
psf_DECpm = PSF_DECpm[PSF_index][0] # in mas/year from GAIA         (example q2237 approx. -13.302)
psf_SNpm = PSF_SN[PSF_index][0]     # totalpm/error from GAIA       (example q2237 approx. 337.6)  
psf_pm_x = -psf_RApm/1000.0         # in arcsec/year
psf_pm_y = psf_DECpm/1000.0         # in arcsec/year
pixscale = 0.387                    # in arcsec/pixel
psfshift_x = psf_pm_x/pixscale      # in pixel/year
psfshift_y = psf_pm_y/pixscale      # in pixel/year
#print(psf_RApm,psf_DECpm,psfshift_x,psfshift_y)

#'true' quasar position estimation:
if q2237_masking == True and q2237 == True:
    # 'true' quasar position estimate with masked data for better median in q2237 --> manuel determined with Robert on Dec. 21, 2022: 
    # mask all data points for median determination where x is smaller than 1757-x_shift, this changes a lot for V and nearly nothing for R:
    mask = np.where(quapos_x+x_shift >= 1757.0)[0] # indices of 'good' values to use only them...
    #print(mask,quapos_x,quapos_x[mask])
    #print(np.median(quapos_x-lin(time,psfshift_x,0.0))+x_shift,np.median(quapos_x[mask]-lin(time[mask],psfshift_x,0.0))+x_shift,len(mask))
    quasar_true_position_x = np.median(quapos_x[mask]-lin(time[mask],psfshift_x,0.0))
    quasar_true_position_x_err = 0.5*iqr(quapos_x[mask]-lin(time[mask],psfshift_x,0.0))/np.sqrt(len(time[mask]))
    quasar_true_position_y = np.median(quapos_y[mask]-lin(time[mask],psfshift_y,0.0))
    quasar_true_position_y_err = 0.5*iqr(quapos_y[mask]-lin(time[mask],psfshift_y,0.0))/np.sqrt(len(time[mask]))
    #print('masking used for better median --> medians changed:')
    quasar_unmasked_pos_x = np.median(quapos_x-lin(time,psfshift_x,0.0))
    quasar_unmasked_pos_y = np.median(quapos_y-lin(time,psfshift_y,0.0))
    N_total = len(time)
    N_masked = N_total - len(time[mask])
    #print('x:',quasar_unmasked_pos_x+x_shift,'-->',quasar_true_position_x+x_shift)
    #print('y:',quasar_unmasked_pos_y,'-->',quasar_true_position_y)
else:
    # original 'true' quasar position estimate by median of psf shift corrected quasar positions with error without masking for the other quasars:
    quasar_true_position_x = np.median(quapos_x-lin(time,psfshift_x,0.0))
    quasar_true_position_x_err = 0.5*iqr(quapos_x-lin(time,psfshift_x,0.0))/np.sqrt(len(time))
    quasar_true_position_y = np.median(quapos_y-lin(time,psfshift_y,0.0))
    quasar_true_position_y_err = 0.5*iqr(quapos_y-lin(time,psfshift_y,0.0))/np.sqrt(len(time))
#print(quasar_true_position_x,quasar_true_position_x_err,quasar_true_position_y,quasar_true_position_y_err)

# final better quasar image A positions for the fifth plot, betterquasarposition.txt and with that the diffgalfitting for better lightcurves!
# the better pos. are linear functions with GAIA-psf-star-slope and median 'true' quasar position as offset, as developed in the 5 plots:
better_qsopos_x = lin(time,psfshift_x,quasar_true_position_x)
better_qsopos_y = lin(time,psfshift_y,quasar_true_position_y)

# first plot - position data and linear fit:
plt.figure(figsize=(10,7),dpi=400)
plt.title('quasar position shift from psf proper motion with linear fit',size=17)
plt.errorbar(time,quapos_x+x_shift,yerr=qua_err_x,fmt='o',markersize=2,elinewidth=1,capsize=3,label='$x$ position')
plt.errorbar(time,quapos_y,yerr=qua_err_y,fmt='o',markersize=2,elinewidth=1,capsize=3,label='$y$ position')
plt.plot(time,lin(time,x_a,x_b),label='lin. fit with slope = '+str(round(x_a,3))+' +/- '+str(round(x_a_err,3)))
plt.plot(time,lin(time,y_a,y_b),label='lin. fit with slope = '+str(round(y_a,3))+' +/- '+str(round(y_a_err,3)))
plt.xlabel('time in years',size=12)
plt.ylabel('positions $x$'+x_shift_sign+str(abs(x_shift))+' and $y$ in pixels',size=12)
plt.legend(fontsize=9,ncol=2,loc=3)
plt.savefig('betterquasarpositions_1.png')
if show_plots == True:
    plt.show()
plt.clf()

# second plot - position data and expectation from GAIA psf proper motion:
plt.figure(figsize=(10,7),dpi=400)
plt.title('quasar position shift and psf star pos. from GAIA proper motions',size=17)
plt.errorbar(time,quapos_x+x_shift,yerr=qua_err_x,fmt='o',markersize=2,elinewidth=1,capsize=3,label='$x$ position')
plt.errorbar(time,quapos_y,yerr=qua_err_y,fmt='o',markersize=2,elinewidth=1,capsize=3,label='$y$ position')
plt.plot(time,lin(time,psfshift_x,quasar_true_position_x+x_shift),label='psf star $x$ pos. (GAIA pm, shifted, NOT a fit!), slope = '+str(round(psfshift_x,3))+' pixel/year',color=cmap(4))
plt.plot(time,lin(time,psfshift_y,quasar_true_position_y),label='psf star $y$ pos. (GAIA pm, shifted, NOT a fit!), slope = '+str(round(psfshift_y,3))+' pixel/year',color=cmap(5))
plt.text(0.005,0.97,'GAIA pm data for psf star (pixel scale = '+str(round(pixscale,3))+\
         ' arcsec/pixel; pmSN[totalpm/error] = '+str(round(psf_SNpm,1))+'; rel. sign between RA and x):'+\
         '\n'+'RApm$='+str(round(psf_RApm,3))+'$ mas/year $='+str(round(psfshift_x,3))+\
         '$ pixel/year and DECpm$='+str(round(psf_DECpm,3))+'$ mas/year $='+str(round(psfshift_y,3))+\
         '$ pixel/year',horizontalalignment='left',verticalalignment='center',transform=plt.gca().transAxes,fontsize=8)
plt.xlabel('time in years',size=12)
plt.ylabel('positions $x$'+x_shift_sign+str(abs(x_shift))+' and $y$ in pixels',size=12)
plt.legend(fontsize=9,ncol=2,loc=3)
plt.savefig('betterquasarpositions_2.png')
if show_plots == True:
    plt.show()
plt.clf()

# third plot - first and second together for comparison:
plt.figure(figsize=(10,7),dpi=400)
plt.title('quasar pos. shift with lin. fit and psf star pos. from GAIA pm',size=17)
plt.errorbar(time,quapos_x+x_shift,yerr=qua_err_x,fmt='o',markersize=2,elinewidth=1,capsize=3,label='$x$ pos.')
plt.errorbar(time,quapos_y,yerr=qua_err_y,fmt='o',markersize=2,elinewidth=1,capsize=3,label='$y$ pos.')
plt.plot(time,lin(time,x_a,x_b),label='lin. fit with slope = '+str(round(x_a,3))+' +/- '+str(round(x_a_err,3)))
plt.plot(time,lin(time,y_a,y_b),label='lin. fit with slope = '+str(round(y_a,3))+' +/- '+str(round(y_a_err,3)))
plt.plot(time,lin(time,psfshift_x,quasar_true_position_x+x_shift),label='$x$ pos. from GAIA psf pm (NOT a fit!), '+str(round(psfshift_x,3))+' pixel/year',color=cmap(4))
plt.plot(time,lin(time,psfshift_y,quasar_true_position_y),label='$y$ pos. from GAIA psf pm (NOT a fit!), '+str(round(psfshift_y,3))+' pixel/year',color=cmap(5))
plt.xlabel('time in years',size=12)
plt.ylabel('positions $x$'+x_shift_sign+str(abs(x_shift))+' and $y$ in pixels',size=12)
plt.legend(fontsize=9,ncol=3,loc=3)
plt.savefig('betterquasarpositions_3.png')
if show_plots == True:
    plt.show()
plt.clf()

# fourth plot - quasar position minus the GAIA psf motion to get the median (with error) of the 'true' unshifted and constant quasar position:
plt.figure(figsize=(10,7),dpi=400)
plt.title('quasar pos. shift minus GAIA psf shift and its median',size=17)
plt.errorbar(time,quapos_x-lin(time,psfshift_x,0.0)+x_shift,yerr=qua_err_x,fmt='o',markersize=2,elinewidth=1,capsize=3,label='$x$ position psf shift corrected')
plt.errorbar(time,quapos_y-lin(time,psfshift_y,0.0),yerr=qua_err_y,fmt='o',markersize=2,elinewidth=1,capsize=3,label='$y$ position psf shift corrected')
plt.plot([np.min(time),np.max(time)],[quasar_true_position_x+x_shift,quasar_true_position_x+x_shift],label='median true $x$ position $\pm$ half iqr$/\sqrt{n}$',color=cmap(4))
plt.plot([np.min(time),np.max(time)],[quasar_true_position_x+x_shift+quasar_true_position_x_err,quasar_true_position_x+x_shift+quasar_true_position_x_err],linestyle='dashed',color=cmap(4))
plt.plot([np.min(time),np.max(time)],[quasar_true_position_x+x_shift-quasar_true_position_x_err,quasar_true_position_x+x_shift-quasar_true_position_x_err],linestyle='dashed',color=cmap(4))
plt.plot([np.min(time),np.max(time)],[quasar_true_position_y,quasar_true_position_y],label='median true $y$ position $\pm$ half iqr$/\sqrt{n}$',color=cmap(5))
plt.plot([np.min(time),np.max(time)],[quasar_true_position_y+quasar_true_position_y_err,quasar_true_position_y+quasar_true_position_y_err],linestyle='dashed',color=cmap(5))
plt.plot([np.min(time),np.max(time)],[quasar_true_position_y-quasar_true_position_y_err,quasar_true_position_y-quasar_true_position_y_err],linestyle='dashed',color=cmap(5))
if q2237_masking == True and q2237 == True:
    plt.text(0.005,0.97,'Excluded '+str(N_masked)+' of '+str(N_total)+' data points from bad seeing images, where the GALFIT-q2237-and-galaxy-model failed, \n'+\
             'which improved the median positions, for x: '+str(round(quasar_unmasked_pos_x+x_shift,2))+' --> '+str(round(quasar_true_position_x+x_shift,2))+\
             ' and y: '+str(round(quasar_unmasked_pos_y,2))+' --> '+str(round(quasar_true_position_y,2))+'!',horizontalalignment='left',verticalalignment='center',transform=plt.gca().transAxes,fontsize=8)
plt.xlabel('time in years',size=12)
plt.ylabel('positions $x$'+x_shift_sign+str(abs(x_shift))+' and $y$ in pixels',size=12)
plt.legend(fontsize=9,ncol=2,loc=3)
plt.savefig('betterquasarpositions_4.png')
if show_plots == True:
    plt.show()
plt.clf()

# fifth plot - propagating the 'true' quasar position with the GAIA psf pm for having the best possible postions for GALFIT on the difference images:
plt.figure(figsize=(10,7),dpi=400)
plt.title('better quasar image A position for difference imaging with GALFIT',size=17)
plt.errorbar(time,quapos_x+x_shift,yerr=qua_err_x,fmt='o',markersize=2,elinewidth=1,capsize=3,alpha=0.5,label='original $x$ position')
plt.errorbar(time,quapos_y,yerr=qua_err_y,fmt='o',markersize=2,elinewidth=1,capsize=3,alpha=0.5,label='original $y$ position')
plt.errorbar(time,better_qsopos_x+x_shift,yerr=quasar_true_position_x_err,fmt='o',markersize=2,elinewidth=1,capsize=3,label='better $x$ position for diff-imag-GALFIT')
plt.errorbar(time,better_qsopos_y,yerr=quasar_true_position_y_err,fmt='o',markersize=2,elinewidth=1,capsize=3,label='better $y$ position for diff-imag-GALFIT')
plt.xlabel('time in years',size=12)
plt.ylabel('positions $x$'+x_shift_sign+str(abs(x_shift))+' and $y$ in pixels',size=12)
plt.legend(fontsize=9,ncol=2,loc=3)
plt.savefig('betterquasarpositions_5.png')
if show_plots == True:
    plt.show()
plt.clf()
    
# save better quasar positions in file to use them with GALFIT on the difference images for getting better lightcurves:
txt = open("betterquasarposition.txt","w")
txt.write("file"+"\t\t\t\t"+"time"+"\t"+"x-pos"+"\t"+"x-pos-err"+"\t"+"y-pos"+"\t"+"y-pos-err"+"\n")
txt.write("\n")
for i in range(len(time)):
    txt.write(str(file[i])+"\t"+str(time[i])+"\t"+str(better_qsopos_x[i])+"\t"+str(quasar_true_position_x_err)+"\t"+str(better_qsopos_y[i])+"\t"+str(quasar_true_position_y_err)+"\n")
txt.close()