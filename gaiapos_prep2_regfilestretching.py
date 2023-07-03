import numpy as np
import sys
import glob
from astropy.io import fits

"""
this program is used to manually improve the radial matching of stars in q2237gaiaregion_pixel.reg ontop of astrometry_ref.fits iteratively!
when finished, goto gaiapos_prep3_....py!
"""

# Change this value and check afterwards, whether most of the stars in astronomy_ref.fits are inside the green circles!
# The regionfile elements (green objects) will be squeezed together radially for POSITIVE values as typically needed!
linear_radial_correction_factor = float(input('Enter radial squeeze-pix-"factor" e.g. 10.0: '))
# At the edge of the image this correspondes approximately to a shift inwards of linear_radial_correction_factor pixels!

#### NEW: Small rotations around the center (in the mathematical positive sence) also possible:
small_rotation_angle = float(input('Enter small rotation angle in deg. e.g. 0.2: ')) # small angle in degrees!!!
if abs(small_rotation_angle) > 1.0:
    print('small_rotation_angle to big! --> set to zero!')
    small_rotation_angle = 0.0
####

# read in pixel <quasar>gaiaregion_pixel.reg file:

gaiaregion_pixelreg_list = glob.glob("*gaiaregion_pixel.reg")
if len(gaiaregion_pixelreg_list) != 1:
    sys.exit('ERROR: not one unambiguously correct *gaiaregion_pixel.reg-file!')
regfilename = gaiaregion_pixelreg_list[0]
    
old_reg = open(regfilename,"r")
lines = old_reg.readlines()
header = ''.join(lines[0:3])

N = int((len(lines)-3)/2)

ra = np.zeros(N)
dec = np.zeros(N)
ring = np.zeros(N)
pmtotal = np.zeros(N)
angle = np.zeros(N)

for i in range(N):
    circle = lines[2*i+3]
    vector = lines[2*i+4]
    circle_list = circle[7:-2].split(",")
    vector_list = vector[9:-11].split(",")
    
    ra[i] = float(circle_list[0])
    dec[i] = float(circle_list[1])
    ring[i] = float(circle_list[2])
    if vector_list[0] != circle_list[0] or vector_list[1] != circle_list[1]:
        sys.exit('ERROR: impossible ra and dec circle vector difference!')
    pmtotal[i] = float(vector_list[2])
    angle[i] = float(vector_list[3])

# change ra and dec pixel positions:
    
reffits = fits.open('astrometry_ref.fits')
shape = reffits[0].data.shape
reffits.close()
ra_center = shape[0]/2.0
dec_center = shape[1]/2.0
center = (ra_center + dec_center)/2.0
if center != ra_center or center != dec_center:
    sys.exit('ERROR: non-quadratic-shape warning for astrometry-ref.fits')

print('Stretching '+str(regfilename)+' by a factor of '+str(round((-1.0)*linear_radial_correction_factor/center,4))+' and rotating it by '+str(round(small_rotation_angle,4))+' degrees.')
# stretching:
ra = ra - linear_radial_correction_factor * (ra-ra_center) / center 
dec = dec - linear_radial_correction_factor * (dec-dec_center) / center
# rotating:
ra = ra - (dec-dec_center) * small_rotation_angle*np.pi/180
dec = dec + (ra-ra_center) * small_rotation_angle*np.pi/180

# write corrected <quasar>gaiaregion_pixel.reg file:

new_reg = open(regfilename,"w")
new_reg.write(header)
for i in range(N):
    new_reg.write('circle('+str(ra[i])+','+str(dec[i])+','+str(ring[i])+')'+'\n'+
              '# vector('+str(ra[i])+','+str(dec[i])+','+str(pmtotal[i])+','+str(angle[i])+') vector=1'+'\n')
new_reg.close()