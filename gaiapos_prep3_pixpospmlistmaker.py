import numpy as np

with open('q2237gaiaregion_pixel.reg') as file:
    lines = file.readlines()
gaialist = [x[9:-11] for x in lines if '# vector(' in x]
#print(gaialist)

ra_pix = np.zeros(len(gaialist))
dec_pix = np.zeros(len(gaialist))
# pm will be taken from original gaia-table, but verified by angle!
angle_check = np.zeros(len(gaialist))
i = 0
for starline in gaialist:
    qwert = starline.split(',')
    ra_pix[i] = qwert[0]
    dec_pix[i] = qwert[1]
    angle_check[i] = qwert[3]
    i = i + 1
#print(ra_pix,dec_pix,angle_check)

pmra,pmdec,pmra_error,pmdec_error = \
np.loadtxt('q2237gaiaDR3table.csv',skiprows=1,delimiter=',',usecols=(13,15,14,16),unpack=True)
#print(pmra,pmdec,pmra_error,pmdec_error)

angle = (180-np.arctan2(pmdec,pmra)*180/np.pi+360)%360
pmtotal = np.sqrt(pmra**2 + pmdec**2)
pmtotal_error = np.sqrt((pmra_error*pmra/pmtotal)**2+(pmdec_error*pmdec/pmtotal)**2)
SN = pmtotal/pmtotal_error
mask = np.where(SN<5,1,0) # mask stars with SN < 5

pmra = np.ma.array(pmra,mask=mask).compressed()
pmdec = np.ma.array(pmdec,mask=mask).compressed()
#pmra_error = np.ma.array(pmra_error,mask=mask).compressed()
#pmdec_error = np.ma.array(pmdec_error,mask=mask).compressed()
#pmtotal = np.ma.array(pmtotal,mask=mask).compressed()
angle = np.ma.array(angle,mask=mask).compressed()
SN = np.ma.array(SN,mask=mask).compressed()
#print(pmra,pmdec,pmra_error,pmdec_error,pmtotal,angle,SN)
print('GAIA position and pm data for',len(SN),'stars!')

# compare two lists: same angles? (ds9 changes it appearently by 180 degrees)
angle_diff = (angle_check-angle+360.0)%360-180.0
#print(angle_diff)
#mean_angle_diff = np.mean(angle_diff)
#print('mean angle difference:',mean_angle_diff)
error_count = 0
for i in range(len(angle_diff)):
    if angle_diff[i] > 0.001:
        error_count = error_count + 1
        print('WARNING: fail at star',i,'with angle difference of',angle_diff[i])
if error_count > 0:
    print('WARNING: potential star mismatching in',error_count,'cases!')
else:
    print(error_count,'angle difference errors! --> angle check successful!')
    
# combine the data and write to final gaia_pixelpositionandpm_list.txt
txt = open("gaia_pixelpositionandpm_list.txt","w")
txt.write('RA [pix]'+'\t'+'\t'+'DEC [pix]'+'\t'+'\t'+'RApm [mas/year]'+'\t'+'\t'+'DECpm [mas/year]'+'\t'+'SN [totalpm/error]'+'\n')
for i in range(len(SN)):
    txt.write(str(ra_pix[i])+'\t'+str(dec_pix[i])+'\t'+str(pmra[i])+'\t'+str(pmdec[i])+'\t'+str(round(SN[i],2))+'\n')
txt.close()

print('Successfully made gaia_pixelpositionandpm_list.txt for gaiaposcorrectionforaligning.py and alining.py!')