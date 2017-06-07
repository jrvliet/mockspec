import numpy as np
from commoncmds import *
import os
import time

timestart = time.time()

kindir = '/home/matrix2/nnielsen/Kinematics/'
#kindir = '../'

dirloc = np.genfromtxt(kindir+'galaxies/isogals.dat',usecols=10,dtype=None)


# Galaxy subsamples
zgal,BK = np.genfromtxt(kindir+'Tables/isoprop.dat',skip_header=2,
                        usecols=(1,7),unpack=True)

Vcirc = np.genfromtxt(kindir+'Tables/massprop.dat',skip_header=2,usecols=5,
                      unpack=True)

# Only use galaxies which have a B-K color
nocolor = []
dirloc2,zgal2,BK2,Vcirc2 = [],[],[],[]
for i in range(len(dirloc)):
    if BK[i] != 9999:
        BK2.append(BK[i])
        dirloc2.append(dirloc[i])
        zgal2.append(zgal[i])
        Vcirc2.append(Vcirc[i])

# Median values of B-K color and zgal
BKmed = np.median(BK2)
zmed = np.median(zgal)

bluelozdir,bluehizdir = [],[]
bluelozgal,bluehizgal = [],[]
bluelozVcirc,bluehizVcirc = [],[]
redlozdir,redhizdir = [],[]
redlozgal,redhizgal = [],[]
redlozVcirc,redhizVcirc = [],[]

for i in range(len(dirloc2)):
    if BK2[i] < BKmed: # Blue galaxies
        if zgal2[i] < zmed: # at low redshift
            bluelozdir.append(dirloc2[i])
            bluelozgal.append(0)
            bluelozVcirc.append(0)
        elif zgal2[i] >= zmed: # at high redshift
            bluehizdir.append(dirloc2[i])
            bluehizgal.append(0)
            bluehizVcirc.append(0)
    elif BK2[i] >= BKmed: # Red galaxies
        if zgal2[i] < zmed: # at low redshift
            redlozdir.append(dirloc2[i])
            redlozgal.append(0)
            redlozVcirc.append(0)
        elif zgal2[i] >= zmed: # at high redshift
            redhizdir.append(dirloc2[i])
            redhizgal.append(0)
            redhizVcirc.append(0)


# Make region dictionary with velocity separations

## blueloz = makedictionary(bluelozdir,bluelozgal,bluelozVcirc)

## bluehiz = makedictionary(bluehizdir,bluehizgal,bluehizVcirc)

redloz  = makedictionary(redlozdir,redlozgal,redlozVcirc)

redhiz  = makedictionary(redhizdir,redhizgal,redhizVcirc)


## time1 = time.time()

## os.system('echo "Blue Low z" > BKz_wrtzabs.log')
## bluelozbinc,bluelozbinnum,bluelozdVi,bluelozmean,bluelozstd = \
##                           bootstrap(blueloz,bluelozgal,bluelozVcirc)

## time2 = time.time()

## os.system('echo "Blue High z" >> BKz_wrtzabs.log')
## bluehizbinc,bluehizbinnum,bluehizdVi,bluehizmean,bluehizstd = \
##                           bootstrap(bluehiz,bluehizgal,bluehizVcirc)


time3 = time.time()
## os.system('echo "Red Low z" >> BKz_wrtzabs.log')
redlozbinc,redlozbinnum,redlozdVi,redlozmean,redlozstd = \
                        bootstrap(redloz,redlozgal,redlozVcirc)


time4 = time.time()
## os.system('echo "Red High z" >> BKz_wrtzabs.log')
redhizbinc,redhizbinnum,redhizdVi,redhizmean,redhizstd = \
                        bootstrap(redhiz,redhizgal,redhizVcirc)

time5 = time.time()


#==============================================================================
# Save the TPCF for gaussian fitting, chi square testing

## blzdn,blzup = updownerrors(bluelozmean,bluelozstd)
## bhzdn,bhzup = updownerrors(bluehizmean,bluehizstd)

rlzdn,rlzup = updownerrors(redlozmean,redlozstd)
rhzdn,rhzup = updownerrors(redhizmean,redhizstd)

## blzsave = []
## for i in range(len(blzdn)):
##     blzsave.append([bluelozbinc[i],bluelozbinnum[i],9999,blzdn[i],blzup[i]])
## np.savetxt('2blueloz_wrtzabs.linbin',blzsave,fmt='%11.5e')
## write_header('\n\n\n','2blueloz_wrtzabs.linbin')

## bhzsave = []
## for i in range(len(bhzdn)):
##     bhzsave.append([bluehizbinc[i],bluehizbinnum[i],9999,bhzdn[i],bhzup[i]])
## np.savetxt('2bluehiz_wrtzabs.linbin',bhzsave,fmt='%11.5e')
## write_header('\n\n\n','2bluehiz_wrtzabs.linbin')

rlzsave = []
for i in range(len(rlzdn)):
    rlzsave.append([redlozbinc[i],redlozbinnum[i],9999,rlzdn[i],rlzup[i]])
np.savetxt('2redloz_wrtzabs.linbin',rlzsave,fmt='%11.5e')
write_header('\n\n\n','2redloz_wrtzabs.linbin')

rhzsave = []
for i in range(len(rhzdn)):
    rhzsave.append([redhizbinc[i],redhizbinnum[i],9999,rhzdn[i],rhzup[i]])
np.savetxt('2redhiz_wrtzabs.linbin',rhzsave,fmt='%11.5e')
write_header('\n\n\n','2redhiz_wrtzabs.linbin')


#==============================================================================
# Save the TPCF as I need it for plotting

## blzbinc,blzfx = shiftarr(bluelozbinc,bluelozbinnum)
## bhzbinc,bhzfx = shiftarr(bluehizbinc,bluehizbinnum)

## dumbinc,blzmean = shiftarr(bluelozbinc,bluelozmean)
## dumbinc,blzstd = shiftarr(bluelozbinc,bluelozstd)

## dumbinc,bhzmean = shiftarr(bluehizbinc,bluehizmean)
## dumbinc,bhzstd = shiftarr(bluehizbinc,bluehizstd)


rlzbinc,rlzfx = shiftarr(redlozbinc,redlozbinnum)
rhzbinc,rhzfx = shiftarr(redhizbinc,redhizbinnum)

dumbinc,rlzmean = shiftarr(redlozbinc,redlozmean)
dumbinc,rlzstd = shiftarr(redlozbinc,redlozstd)

dumbinc,rhzmean = shiftarr(redhizbinc,redhizmean)
dumbinc,rhzstd = shiftarr(redhizbinc,redhizstd)


## blzstddn,blzstdup = updownerrors(blzmean,blzstd)
## bhzstddn,bhzstdup = updownerrors(bhzmean,bhzstd)
rlzstddn,rlzstdup = updownerrors(rlzmean,rlzstd)
rhzstddn,rhzstdup = updownerrors(rhzmean,rhzstd)


## bluelozsave = []
## for i in range(len(blzbinc)):
##     bluelozsave.append([blzbinc[i],blzfx[i],blzstddn[i],blzstdup[i]])

## bluehizsave = []
## for i in range(len(bhzbinc)):
##     bluehizsave.append([bhzbinc[i],bhzfx[i],bhzstddn[i],bhzstdup[i]])

## np.savetxt('2blueloz_wrtzabs.tpcf',bluelozsave,fmt='%9.5e')
## np.savetxt('2bluehiz_wrtzabs.tpcf',bluehizsave,fmt='%9.5e')


redlozsave = []
for i in range(len(rlzbinc)):
    redlozsave.append([rlzbinc[i],rlzfx[i],rlzstddn[i],rlzstdup[i]])

redhizsave = []
for i in range(len(rhzbinc)):
    redhizsave.append([rhzbinc[i],rhzfx[i],rhzstddn[i],rhzstdup[i]])

np.savetxt('2redloz_wrtzabs.tpcf',redlozsave,fmt='%9.5e')
np.savetxt('2redhiz_wrtzabs.tpcf',redhizsave,fmt='%9.5e')

#==============================================================================

timeend = time.time()

print (timeend-timestart)/60.0,'Total time (all in minutes)'
## print (time2-time1)/60.0
## print (time3-time2)/60.0
print (time3-timestart)/60.0,'Time to create dictionaries'
print (time4-time3)/60.0,'Time to bootstrap red, low z'
print (time5-time4)/60.0,'Time to bootstrap red, high z'
print (timeend-time5)/60.0,'Time to save *.linbin and *.tpcf files'

print 'Median redshift, B-K: ',zmed,BKmed

## print len(blueloz),len(bluehiz),len(redloz),len(redhiz)
print 'Number of absorbers in red low z, high z samples: ',len(redloz),\
      len(redhiz)

## print len(bluelozdVi),len(bluehizdVi),len(redlozdVi),len(redhizdVi)
print 'Number of velocity splittings in red low z, high z samples: ',\
      len(redlozdVi),len(redhizdVi)
