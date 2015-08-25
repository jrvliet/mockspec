
# Functions to generate lines of sight

from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from numpy import linalg as LA
import math
import sys
from random import random, seed
from os import popen
from numpy import linspace as ls
from numba import jit
import subprocess as sp

def inBox(x,y,z,boxsize):
    size = boxsize
    if x<size and x>-1*size and y<size and y>-1*size and z<size and z>-1*size:
        return 1
    else:
        return 0


@jit
def findEnds(px, py, pz, dx, dy, dz, boxsize):
    size = boxsize/2.0
    xen, yen, zen, ten, xex, yex, zex, tex = 0,0,0,0,0,0,0,0
    t = ls(-5e3, 5e3, 5e5)
    
    for i in range(0,len(t)-1):
        x1 = px + dx*t[i]
        y1 = py + dy*t[i]
        z1 = pz + dz*t[i]
        x2 = px + dx*t[i+1]
        y2 = py + dy*t[1+i]
        z2 = pz + dz*t[i+1]
        
        if inBox(x1, y1, z1, size)==0 and inBox(x2, y2, z2, size)==1:
            xen = (x1+x2)/2.0
            yen = (y1+y2)/2.0
            zen = (z1+z2)/2.0
            ten = (t[i]+t[i+1])/2.0
        if inBox(x1, y1, z1, size)==1 and inBox(x2, y2, z2, size)==0:
            xex = (x1+x2)/2.0
            yex = (y1+y2)/2.0
            zex = (z1+z2)/2.0
            tex = (t[i]+t[i+1])/2.0
    return xen, yen, zen, ten, xex, yex, zex, tex

def func(x, y):
    return 0.0

def read_control_file(filename):

    f = open(filename)
    
    for i in range(0,5):
        line = f.readline()

    gasfile = line.split()[0]
    gasfile = gasfile.replace('"', '')
    galID = gasfile.split('_')[0]
    
    line = f.readline()
    rootname = line.split()[0]
    
    line = f.readline()
    aexpn = float(line.split()[0])
    
    for i in range(0,3):
        line = f.readline()
    summaryLoc = (line.split()[0]).replace('"','')

    f.close()

    return gasfile, galID, rootname, aexpn, summaryLoc


def read_summary(galID, aexpn, summaryLoc):

    import sys

    aexpn = float(aexpn)

    # Location of summaries
    f = open(summaryLoc + galID + '.dat')

    f.readline()
    f.readline()

    found = 0

    for line in f:

        l = line.split()
        a = float(l[0])

        if a==aexpn:
            
            found = 1

            # Read in data
            mvir = float(l[2])
            rvir = float(l[3])
            a11 = float(l[4])
            a12 = float(l[5])
            a13 = float(l[6])
            a21 = float(l[7])
            a22 = float(l[8])
            a23 = float(l[9])
            a31 = float(l[10])
            a32 = float(l[11])
            a33 = float(l[12])
            vpec_x = float(l[13])
            vpec_y = float(l[14])
            vpec_z = float(l[15])
            
    f.close()
    if found==1:
        return mvir, rvir, a11, a12, a13, a21, a22, a23, a31, a32, a33, vpec_x, vpec_y, vpec_z
    else:
        print 'ERROR: Could not find the expansion factor {0} in file: \n\t'.format(aexpn) + summaryLoc
        sys.exit()






def genLines(galID, gasfile, summaryLoc, expn, inc, nLOS, maximpact, ncores):

    tol = 1e-5
    seed(25525)
   
    # Read in galaxy properties
    mvir, rvir, a11, a12, a13, a21, a22, a23, a31, a32, a33, vpec_x, vpec_y, vpec_z = read_summary(galID, expn, summaryLoc)
    maximpact_kpc = maximpact * rvir
    
    # Initialize variables
    Xcom, Ycom, Zcom = 0., 0., 0.
    VXcom, VYcom, VZcom = 0., 0., 0.
    l, b = 0.0, 90.0
    vx_obs, vy_obs, vz_obs = 0., 0., 0.
     
    boxsize = 4.0*rvir

    # Open output files
    try:
        fout = open('lines.dat','w')
        finfo = open('lines.info','w')
    except IOError:
        print 'Error opening files lines.dat, lines.info for writing.'
        print 'Check permissions'
        print 'Exitting...'
        sys.exit()

    # Write headers
    fout.write('#        Enter points                      Exit Points\n')
    fout.write('# xen            yen         zen         xex         yex         zex\n')
    finfo.write('# More details on each LOS\n')
    finfo.write('# LOS num     b(kpc)      phi     Inclination\n')

    # Define the rotation matrix to convert from the box frame to the galaxy frame
    a_btg = np.matrix([[a11, a12, a13], [a21, a22, a23], [a31, a32, a33]])

    dxs = []
    dys = []
    dzs = []
    xens = []
    yens = []
    zens = []
    xexs = []
    yexs = []
    zexs = []


    # Generate random impact parameters:
    impacts = []
    for j in range(0,nLOS-1):
        b = maximpact_kpc * random()
        impacts.append(b)

    # Sort impact parameters
    impacts.sort()



    for j in range(0,nLOS):


        b = impacts[j]
        phi = random()*2*math.pi
        
        finfo.write('{0:>6}       {1:>6,.1f}    {2:>6,.1f}   {3:>6,.1f}\n'.format(j+1, b, math.degrees(phi), math.degrees(inc)))
           
        # Define the LOS directional vector in the sky frame:
        ds = np.matrix([0.0, 0.0, -1.0]).T
        
        # Define the impact point in the sky frame:
        ps = np.matrix([b*math.cos(phi), b*math.sin(phi), 0.0]).T


        # Define the rotation matrix to convert from the sky frame to the galaxy frame
        a_stg = np.matrix([[1,0,0], [0, math.cos(inc), -1*math.sin(inc)], [0, math.sin(inc), math.cos(inc)]])

        # Rotate the sky's directional vector and impact point into the galaxy's frame
        dg = a_stg*ds
        pg = a_stg*ps


        # To rotate from the sky's frame to the box frame, need a_gtb, which is the 
        # inverse of a_btg
        a_gtb = LA.inv(a_btg) 
        
        # Rotate the galaxy's directional vector and impact point into the box's frame
        db = a_gtb*dg
        pb = a_gtb*pg


        # Get the enter and exit points of the LOS
        xen, yen, zen, ten, xex, yex, zex, tex = findEnds(pb[0,0], pb[1,0], pb[2,0], db[0,0], db[1,0], db[2,0], boxsize)
        

        fout.write('{0:>12,.5f}{1:>12,.5f}{2:>12,.5f}{3:>12,.5f}{4:>12,.5f}{5:>12,.5f}\n'.format(xen,yen,zen,xex,yex,zex))

        xens.append(xen)
        yens.append(yen)
        zens.append(zen)
        xexs.append(xex)
        yexs.append(yex)
        zexs.append(zex)

    # Determine the inclination of the galaxy
    # Use the dot product of the normal vector and the LOS
    norm_g = np.matrix([0,0,10])
    norm_b = a_gtb*norm_g.T
    los_len = math.sqrt(db[0,0]*db[0,0] + db[1,0]*db[1,0] + db[2,0]*db[2,0])
    norm_len = math.sqrt(norm_b[0,0]*norm_b[0,0] + norm_b[1,0]*norm_b[1,0] + norm_b[2,0]*norm_b[2,0])
    dotprod = db[0,0]*norm_b[0,0] + db[1,0]*norm_b[1,0] + db[2,0]*norm_b[2,0]
    interior = dotprod / (los_len*norm_len)

    # The following conditional statements are needed to account 
    # for rounding errors that occur
    if interior<-1.0 and interior>-1.0-tol:
        interior=-1.0
    elif interior>1.0 and interior<1.0+tol:
        interior=1.0
    angle = math.acos(interior)
    if angle>math.pi/2.0 and angle<math.pi:
        angle = angle - math.pi/2.0

    print 'Calculated inclination: {0:.3f} degrees'.format(math.degrees(angle))

        
    fout.close()
    finfo.close()
   


 
def runCellfinder(codeLoc):
    command = codeLoc+'/funcs/cellfinder/cellfinder'
    try: 
        sp.check_call(command, shell=True)
    except:
        print '\n\nCould not run cellfinder with:\n\t{0:s}'.format(command)
        print 'Exiting...'
        sys.exit()
