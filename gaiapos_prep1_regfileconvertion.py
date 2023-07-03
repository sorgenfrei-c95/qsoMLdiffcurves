import numpy as np

#### BEFORE STARTING gaiapospm_regionfileconvertion.py ####

"""
Newer information follows this from line 21 on!
"""

# input ADQL Query at https://dc.zah.uni-heidelberg.de/__system__/adql/query/form
# to get the needed .csv-tabel of star-positions and propermotions close to the quasar; example input for q2237:
# 
# SELECT *
# FROM gaia.edr3lite
# WHERE 1=CONTAINS(
#   POINT('ICRS',340.126,3.358),
#   CIRCLE('ICRS',ra, dec, 0.4))
# AND ruwe < 1.4
#
# and to be save set timeout to 10s and limit to 10000 (at least).

"""
August 2022: switching from EDR3 to DR3!!! Access with ESAs GAIA archiv at https://gea.esac.esa.int/archive/
--> go to 'Search' and then 'Advanced(ADQL)': there choose job name, set download format to CSV and input query:

ADQL-example for q2237:

SELECT *
FROM gaiadr3.gaia_source
WHERE 1=CONTAINS(
  POINT('ICRS',340.126,3.358),
  CIRCLE('ICRS',ra, dec, 0.4))
AND ruwe < 1.4

--> then move csv into the two band-folders of the quasar-folder and rename it <quasar>gaiaDR3table.csv!
"""

#### read in .csv-tabel with gaia-data and mask bad SN stars ####

"""
first row should start with:
solution_id,designation,source_id,random_index,ref_epoch,ra,ra_error,dec,dec_error,parallax,parallax_error,parallax_over_error,pm,pmra,pmra_error,pmdec,pmdec_error,...
--> usecols=(5,7,6,8,13,15,14,16,12)
units: ra & dec: [deg] ; ra_error & dec_error: [mas] ; everything pm-related: [mas/yr]
"""

ra,dec,ra_error,dec_error,pmra,pmdec,pmra_error,pmdec_error,totalpmcheck = \
np.loadtxt('q2237gaiaDR3table.csv',skiprows=1,delimiter=',',usecols=(5,7,6,8,13,15,14,16,12),unpack=True) # for old ARI-DataCenter-EDR3lite-q2237gaiatable.csv: usecols=(1,2,3,4,5,6,7,8) without totalpmcheck
#print(ra,dec,ra_error,dec_error,pmra,pmdec,pmra_error,pmdec_error,totalpmcheck)

angle = (180-np.arctan2(pmdec,pmra)*180/np.pi+360)%360
pmtotal = np.sqrt(pmra**2 + pmdec**2)
if np.max(pmtotal-totalpmcheck) > 0.001:
    print('WARNING!!! POTENTIALLY FATAL ERORR!!!')
pmtotal_error = np.sqrt((pmra_error*pmra/pmtotal)**2+(pmdec_error*pmdec/pmtotal)**2)
SN = pmtotal/pmtotal_error
mask = np.where(SN<5,1,0) # mask stars with SN < 5

ra = np.ma.array(ra,mask=mask).compressed()
dec = np.ma.array(dec,mask=mask).compressed()
ra_error = np.ma.array(ra_error,mask=mask).compressed()
dec_error = np.ma.array(dec_error,mask=mask).compressed()
pmra = np.ma.array(pmra,mask=mask).compressed()
pmdec = np.ma.array(pmdec,mask=mask).compressed()
pmra_error = np.ma.array(pmra_error,mask=mask).compressed()
pmdec_error = np.ma.array(pmdec_error,mask=mask).compressed()
pmtotal = np.ma.array(pmtotal,mask=mask).compressed()
angle = np.ma.array(angle,mask=mask).compressed()
SN = np.ma.array(SN,mask=mask).compressed()

#### write .reg-file, i.e. a region file to be opend with ds9 ontop ####

reg = open("q2237gaiaregion.reg","w")
reg.write('# Region file format: DS9 version 4.1'+'\n'+
'global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'+'\n'+
'fk5'+'\n')
for i in range(len(SN)):
    reg.write('circle('+str(ra[i])+','+str(dec[i])+',3.0")'+'\n'+
              '# vector('+str(ra[i])+','+str(dec[i])+','+str(pmtotal[i])+'",'+str(angle[i])+') vector=1'+'\n')
reg.close()

# reminder: unit of pm from GAIA is mas/year

"""
#### AFTER EXECUTING gaiapos_prep1_regfileconvertion.py ####

# open the q2237gaiaregion.reg-file ontop of astronomy_ref.fits in ds9:
    
ds9 astrometry_ref.fits -region q2237gaiaregion.reg

# and if necessary shift the region (strg A while in edit region, then use arrows) to fit the stars as good as possible
# because they will be matched with isis positions and then save the regions but in PIXEL UNITS!!!!!!!!!
---> save it (in pixel units!!! coord. sys: IMAGE) as <quasar>gaiaregion_pixel.reg and go to gaiapos_prep2_....py
"""