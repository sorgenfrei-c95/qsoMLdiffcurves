import numpy as np
import os
from pathlib import Path
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from scipy.stats import iqr
import sys

local = (str(os.path.abspath(os.getcwd())))
path = Path(local)
print('local path:', path)

reftime = str(input('Enter time intervall of refimage (e.g. "2014-22"): '))

reffile = list(path.glob('ref*'+reftime+'.fits'))
reffile_error = list(path.glob('error_ref*'+reftime+'.fits'))
if len(reffile) != 1 or len(reffile_error) != 1:
    print('#### REFERENCE FILE ERROR ####')
    sys.exit('Interrupting galref_apparentmagoffset.py')

def make_galfit_input_file(data_file_name, error_file_name, object_region, position, fluxesti, second_fit):
    galfit_input_txt = open("galfit.input", "w")
    galfit_input_txt.write("A) "+str(data_file_name)+" \n" +
                           "B) out.fits "+"\n" +
                           "C) "+str(error_file_name)+" \n" +
                           "D) psf.fits "+"\n" +
                           "E) 1 "+"\n" +
                           "F)   "+"\n" +
                           "G)   "+"\n" +
                           "H) "+str(object_region)+" \n" +
                           "I) 51    51 "+"\n" +
                           "J) 25.000 "+"\n" +
                           "K) 0.387  0.387 "+"\n" +
                           "O) regular "+"\n" +
                           "P) 0 "+"\n" +
                           "S) 0 "+"\n" +
                           "0)      q2237  "+"\n" +
                           " 1) "+str(position)+" 1 1  "+"\n" +
                           " 3) "+str(fluxesti)+"     1 "+"\n" +
                           " 4) "+str(second_fit)+"\n" +
                           " 5)      0.0   0 "+"\n" +
                           " 6)      0.0   0 "+"\n" +
                           " 7) 0.0000     0 "+"\n" +
                           " 8) 1.0000     0 "+"\n" +
                           " 9) 0.0000     0 "+"\n" +
                           "10) 0.0000     0 "+"\n" +
                           " Z) 0            "+"\n" +
                           " 0)        sky   "+"\n" +
                           " 1) 10.0       1 "+"\n" +
                           " 2) 0.0000     0 "+"\n" +
                           " 3) 0.0000     0 "+"\n" +
                           " Z) 0 ")
    galfit_input_txt.close()

# ref-stars:
# initial guesses for star positions for aperture photometry for scale factors:
# q2237:


if str(path).find('q2237') != -1:
    print('refstars of quasar: q2237')
    initial_positions = np.array([[1578.0, 1825.0], [1600.0, 1586.0], [2066.0, 1498.0], [2203.0, 1689.0], [1890.0, 1904.0],
                                  [1837.0, 2049.0], [1677.0, 1991.0], [1271.0, 1603.0], [2191.0, 2213.0], [2421.0, 2327.0],
                                  [1721.0, 2276.0], [1505.0, 2288.0], [1339.0, 2281.0], [1052.0, 1921.0], [ 812.0, 1835.0],
                                  [ 874.0, 1565.0], [1123.0, 1287.0], [1408.0, 1368.0], [2055.0, 1346.0], [2732.0, 1166.0],
                                  [3266.0, 2570.0]])
# he2149:
if str(path).find('he2149') != -1:
    print('refstars of quasar: he2149')
    initial_positions = np.array([[2013.0, 1854.0], [1274.0, 1823.0], [1280.0, 1650.0], [2150.0, 1518.0], [2239.0, 1552.0],
                                  [2199.0, 1889.0], [2482.0, 1831.0], [1851.0, 1891.0], [1086.0, 1854.0], [1043.0, 1618.0],
                                  [ 993.0, 1574.0], [ 959.0, 1523.0], [1190.0, 1986.0], [ 507.0, 1928.0], [ 656.0, 1384.0],
                                  [1602.0, 1381.0], [1823.0, 1359.0], [2387.0, 1375.0], [2898.0, 2169.0], [2733.0, 2107.0],
                                  [2318.0, 2246.0], [1793.0, 2304.0], [1595.0, 1049.0], [1504.0,  763.0]])
# he1104:
if str(path).find('he1104') != -1:
    print('refstars of quasar: he1104')
    initial_positions = np.array([[1816.0, 1649.0], [1847.0, 1657.0], [2080.0, 1895.0], [2159.0, 1977.0], [1354.0, 1910.0],
                                  [1210.0, 1983.0], [1662.0,  882.0], [ 977.0,  900.0], [2814.0, 1040.0], [3018.0, 1638.0],
                                  [2849.0, 2156.0], [2228.0, 2063.0], [1742.0, 2317.0], [1304.0, 2398.0], [ 364.0, 2267.0],
                                  [ 992.0, 1694.0], [1087.0, 1443.0], [1826.0, 1707.0], [1805.0, 1760.0], [2342.0,  729.0],
                                  [ 782.0, 1648.0], [3082.0,  737.0]])

# converts to ['abcd wxyz',...] where abcd and wxyz are the ds9 pixel positions (in the order given by ds9):
ref_positions = [str(np.array(x, dtype=int))[1:-1]
                 for x in initial_positions.tolist()]

# print(ref_positions)

# ref_regions:
halfside = 15  # size (half side) of region around ref-stars
ref_regions = []
for position in ref_positions:
    xy_pos = position.split()
    center = np.floor(
        np.array([float(xy_pos[0]), float(xy_pos[1])])).astype(int)
    edge1, edge2, edge3, edge4 = center[0]-halfside, center[0] + \
        halfside, center[1]-halfside, center[1]+halfside
    region = str(edge1)+' '+str(edge2)+' '+str(edge3)+' '+str(edge4)
    ref_regions.append(region)
# print(ref_regions)

reffile = str(reffile[0]).replace(str(path)+'/', '')
reffile_error = str(reffile_error[0]).replace(str(path)+'/', '')
print('reffile:', reffile)
data = fits.open(reffile)
# remove possible remaining background from ref-image:
mean, median, std = sigma_clipped_stats(
    data[0].data[100:-100, 100:-100], sigma=3.0, std_ddof=1)
background = median
# std i.e. typical fluctuation size of the individual background pixels
error_background = std
# Der Fehler des Hintergrunds wurde im ursprünglichen combining schon im Fehlerbild berücksichtigt (nur in den diff-Versuchen manchmal wieder das negative minimum im Bild dazuaddiert und abgespeichert)
print('ref-image background =', background, '+/-', error_background)
refimage_backred = data[0].data - background
print('make psf.fits')
if str(path).find('q2237') != -1:
    psf_cutout = refimage_backred[2120:2150, 2213:2243]
if str(path).find('he2149') != -1:
    psf_cutout = refimage_backred[1703:1733, 1786:1816]
if str(path).find('he1104') != -1:
    psf_cutout = refimage_backred[1770:1800, 1700:1730]
psf_fits = fits.PrimaryHDU(psf_cutout)
psf_fits = fits.HDUList([psf_fits])
psf_fits[0].data = psf_cutout
psf_fits.writeto('psf.fits', overwrite=True)

ref_fluxes = np.zeros(len(ref_regions))
print('reference star measurements:')
for ref_star in range(len(ref_regions)):
    print('refstar no.', str(ref_star+1))
    refstar_region = ref_regions[ref_star]
    refstar_position = ref_positions[ref_star]
    no_second_fit = '     0.0   0'
    fluxestimate = '50000.0'
    make_galfit_input_file(reffile, reffile_error, refstar_region,
                           refstar_position, fluxestimate, no_second_fit)
    # execute galfit with no console-output
    os.system('./galfit galfit.input > /dev/null')
    os.system('rm galfit.??')  # remove extra files produced by galfit
    os.system('rm out.fits')  # remove output fits file
    try:
        with open('fit.log') as f:
            galfit_ref_results = f.readlines()[7:9]
            # print(galfit_ref_results)
        os.system('rm fit.log')  # remove fit.log
    except:
        print('refstar-galfit-error with refstar',
              str(ref_star+1), '! --> skipped!')
        continue
    ref_star_flux = galfit_ref_results[0].split()[4]
    if ')' in ref_star_flux:
        ref_fluxes[ref_star] = galfit_ref_results[0].split()[5]
    elif ',' in ref_star_flux:
        ref_fluxes[ref_star] = galfit_ref_results[0].split()[6]
    else:
        ref_fluxes[ref_star] = ref_star_flux

print('refstar fluxes:',ref_fluxes)
ref_mag = -2.5*np.log10(ref_fluxes)
print('calculated refstar magnitudes:',ref_mag)

# Compare refstar magnitudes with true values from GAIA data:


def GAIA_GminusV(Gbp, Grp):
    x = Gbp-Grp
    return -0.02704+0.01424*x-0.2156*x**2+0.01426*x**3


def GAIA_GminusR(Gbp, Grp):
    x = Gbp-Grp
    return -0.02275+0.3961*x-0.1243*x**2-0.01396*x**3+0.003775*x**4


def GAIA_V(GminusV, G):
    return G-GminusV


def GAIA_R(GminusR, G):
    return G-GminusR

if str(path).find('q2237') != -1:
    GAIAdata_table = '../q2237gaiaDR3table.csv'
if str(path).find('he2149') != -1:
    GAIAdata_table = '../he2149gaiaDR3table.csv'
if str(path).find('he1104') != -1:
    GAIAdata_table = '../he1104gaiaDR3table.csv'

G,Gbp,Grp,compareRApm,compareDECpm = np.genfromtxt(GAIAdata_table,dtype=float,skip_header=1,delimiter=',',usecols=(69,74,79,13,15),unpack=True)
# print(G,Gbp,Grp)

if str(path).find('/V/') != -1:
    GAIA_ref_mag = GAIA_V(GAIA_GminusV(Gbp, Grp), G)

if str(path).find('/R/') != -1:
    GAIA_ref_mag = GAIA_R(GAIA_GminusR(Gbp, Grp), G)

# print(GAIA_ref_mag)

RA,DEC,RApm,DECpm = np.loadtxt('../gaia_pixelpositionandpm_list.txt',skiprows=1,usecols=(0,1,2,3),unpack=True)
# print(RA,DEC)
# print(len(RA),len(GAIA_ref_mag))

counter = 0
for i in range(len(compareRApm)):
    pmdiff = np.sqrt((compareRApm[i]-RApm)**2+(compareDECpm[i]-DECpm)**2)
    mindex = np.argmin(pmdiff)
    if pmdiff[mindex] == 0.0:
        # pmdiff[mindex] should be zero
        #print(pmdiff[mindex],compareRApm[i],RApm[mindex],compareDECpm[i],DECpm[mindex])
        counter = counter + 1
GAIA_RA = np.zeros(counter)
GAIA_DEC = np.zeros(counter)
GAIA_MAG = np.zeros(counter)
counter2 = 0
for i in range(len(compareRApm)):
    pmdiff = np.sqrt((compareRApm[i]-RApm)**2+(compareDECpm[i]-DECpm)**2)
    mindex = np.argmin(pmdiff)
    if pmdiff[mindex] == 0.0:
        GAIA_RA[counter2] = RA[mindex]
        GAIA_DEC[counter2] = DEC[mindex]
        GAIA_MAG[counter2] = GAIA_ref_mag[i]
        counter2 = counter2 + 1

#print(GAIA_RA,GAIA_DEC,GAIA_MAG)

# match GAIA data to refstars:
ref_mag_GAIA = np.zeros(len(ref_mag))
match_counter = 10.0 # distance in pixels must be below this to have a match
trimmedpixel = 100 # between aligning and combining/galfitting
for i in range(len(ref_mag)):
    # starpmRA/DEC_forgooddata will be written to good.data and:
    distance = np.sqrt((initial_positions[i][0]-(GAIA_RA-trimmedpixel))**2+(initial_positions[i][1]-(GAIA_DEC-trimmedpixel))**2)
    mindex = np.argmin(distance)
    print('distance of best match:',distance[mindex])
    if distance[mindex] < 15.0:                
        #print(distance[mindex],initial_positions[i][0],GAIA_RA[mindex]-trimmedpixel,initial_positions[i][1],GAIA_DEC[mindex]-trimmedpixel)
        ref_mag_GAIA[i] = GAIA_MAG[mindex]
    else:
        print('no match --> mag set to NAN')
        ref_mag_GAIA[i] = np.nan

print('refstar magnitudes from GAIA:',ref_mag_GAIA)
diff_mag = ref_mag_GAIA-ref_mag
print('difference of refstar magnitudes GAIA-GALFIT:',diff_mag)
print('####')
offset = np.nanmedian(diff_mag)
offerr = iqr(diff_mag,nan_policy='omit')/2.0
print('magnitude offset, i.e. median of the differences (with iqr error):')
print('offset =',offset,'+/-',offerr)
print('####')

txt = open("refimage_"+reftime+"_magnitude_offset.txt", "w")
txt.write("offset"+"\t\t\t"+"error"+"\n")
txt.write(str(offset)+"\t"+str(offerr))
txt.close()
print('offset written in txt-file:',"refimage_"+reftime+"_magnitude_offset.txt") 
print('galref_apparentmagoffset.py is finished!')