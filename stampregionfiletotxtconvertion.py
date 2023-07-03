import numpy as np

reg = open("stamps_lco.reg","r")
reglines = reg.readlines()
reg.close()

stamp_x = []
stamp_y = []

for i in range(len(reglines)):
    if reglines[i][0:4] == 'box(':
        upto = reglines[i].find(')')
        stamppixel = reglines[i][4:upto].split(',')
        stamp_x.append(stamppixel[0])
        stamp_y.append(stamppixel[1])

Nx = len(stamp_x)
Ny = len(stamp_y)
if Nx == Ny:
    N = Nx
    #print(stamp_x,stamp_y,N)
else:
    print('#### FATAL ERROR!!! ####')

txt = open("stamps_lco.txt","w")

for j in range(N):
    txt.write(stamp_x[j]+' '+stamp_y[j]+'\n')
    
txt.close()

unique_stamps = np.unique(np.loadtxt("stamps_lco.txt"),axis=0)

txt = open("stamps_lco.txt","w")

for j in range(int(unique_stamps.size/2)):
    txt.write(str(int(unique_stamps[j][0]))+' '+str(int(unique_stamps[j][1]))+'\n')

txt.close()

print('stamps_lco.txt written!')