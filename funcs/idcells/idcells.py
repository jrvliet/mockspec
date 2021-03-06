#!/usr/bin/python

# Reads in list of cells found by FBcellfinder.c 
# Finds the cells in the ion boxes
# Creates the los list for each ion

# Usage:
#  idcells.py <gal_props file>

import sys
import numpy as np
import subprocess as sp
import pandas as pd
import os

class losProps(object):

    def __init__ (self):

        self.xen = 0.0
        self.yen = 0.0
        self.zen = 0.0
        self.losx = 0.0 
        self.losy = 0.0 
        self.losz = 0.0 
        self.a11 = 0.0
        self.a12 = 0.0
        self.a13 = 0.0
        self.a21 = 0.0
        self.a22 = 0.0
        self.a23 = 0.0
        self.a31 = 0.0
        self.a32 = 0.0
        self.a33 = 0.0
        self.Xcom = 0.0
        self.Ycom = 0.0
        self.Zcom = 0.0
        self.VXcom = 0.0
        self.VYcom = 0.0
        self.VZcom = 0.0
        self.x0 = 0.0
        self.y0 = 0.0
        self.z0 = 0.0
        self.vx_obs = 0.0
        self.vy_obs = 0.0
        self.vz_obs = 0.0


def read_control_file(filename):

    f = open(filename)
    
    for i in range(0,4):
        line = f.readline()

    gasfile = line.split()[0]
    galID = gasfile.split('_')[0]
    aexpn_tmp = (gasfile.split('GZa')[1])
    aexpn = aexpn_tmp.split('.')[0] + '.' + aexpn_tmp.split('.')[1]
    
    line = f.readline()
    summaryLoc = (line.split()[0])
    
    line = f.readline()
    ion_list = line.split('#')[0].rstrip().split()

    ion_num = len(ion_list)

    f.close()

    return galID, aexpn, ion_list, ion_num


def read_lines(filename):
    
    f = open(filename)

    los_num, los_b, los_phi = [], [], []

    for line in f:

        if '#' not in line:
            
            l = line.split()
            los_num.append(int(l[0]))
            los_b.append(float(l[1]))
            los_phi.append(float(l[2]))

    return los_num, los_b, los_phi

    
#def write_OutfileHdr(outfile, aexpn, R0, phi, l, b, xen, yen, zen, losx, losy, losz, 
#                        a11, a12, a13, a21, a22, a23, a31, a32, a33, Xcom, Ycom, Zcom, 
#                        VXcom, VYcom, VZcom, x0, y0, z0, vx_obs, vy_obs, vz_obs):
def write_OutfileHdr(outfile, aexpn, R0, phi, l, b, los):

    # Write header information to output file

    from math import sqrt

    f = open(outfile,'w')
    
    s = ('aexpn   {0:.3f} Galaxy center in kpc:  {1:.6f}  {2:.6f}  {3:.6f} '
         'Peculiar velocity center in km/s:    {4:.3f}   {5:.3f} '
         '{6:.3f}\n'.format(float(aexpn), los.Xcom, los.Ycom, los.Zcom,
                            los.VXcom, los.VYcom, los.VZcom))
    f.write(s)
    
    s = ('Galactic Coordinates (l:b) in degrees:     {0:.1f}    {1:.1f}  '
        'Observer Position (R(kpc):phi(deg)):     {2:.6f} '
        '{3:.6f}\n'.format(l,b,R0,phi))
    f.write(s)
    
    s = ('Observer Velocity (vx vy vz) in km/s in Box Frame:        '
        '{0:.3E}      {1:.3E}      {2:.3E}'.format(los.vx_obs,
                                            los.vy_obs,los.vz_obs))
    f.write(s)
    f.write('\n')

    s = ('Box Entry Point (x0:y0:z0) in kpc:      {0:.6f}     {1:.6f}       '
        '{2:.6f}  LOS Direction Cosines (l:m:n):    {3:.6E}   {4:.6E}  '
        '{5:.6E}'.format(los.xen,los.yen,los.zen,
                         los.losx,los.losy,los.losz))
    f.write(s)
    f.write('\n')

    s = ('Z prime axis and Lvec/L (a31 a32 a33): '
        '{0: .7E} {1: .7E} {2: .7E}\n'.format(los.a31,los.a32,los.a33))
    f.write(s)

    s = ('X prime axis            (a11 a12 a13): '
        '{0: .7E} {1: .7E} {2: .7E}\n'.format(los.a11,los.a12,los.a13))
    f.write(s)

    s = ('Y prime axis            (a21 a22 a23): '
        '{0: .7E} {1: .7E} {2: .7E}\n'.format(los.a21,los.a22,los.a23))
    f.write(s)
    
    width = 15

    s = ('1'.center(width) + '2'.center(width) + '3'.center(width) + 
           '4'.center(width) + '5'.center(width) + '6'.center(width) + 
           '7'.center(width) + '8'.center(width) + '9'.center(width) + 
           '10'.center(width) + '11'.center(width) + '12'.center(width) + 
           '13'.center(width) + '14'.center(width) + '15'.center(width)+ '\n')
    f.write(s)

    s = ('{0:^15}'.format('cell size') + '{0:^15}'.format('x') + '{0:^15}'.format('y') + 
        '{0:^15}'.format('z') + '{0:^15}'.format('vx') + '{0:^15}'.format('vy') + 
        '{0:^15}'.format('vz') + '{0:^15}'.format('nH') + '{0:^15}'.format('temperature') + 
        '{0:^15}'.format('SNII mass frac') + '{0:^15}'.format('SNIa mass frac') + 
        '{0:^15}'.format('natom') + '{0:^15}'.format('fion') + '{0:^15}'.format('nion') + 
        '{0:^15}\n'.format('cell id'))
    f.write(s)
   
    s = ('(kpc)'.center(width) + '(kpc)'.center(width) + '(kpc)'.center(width) + 
        '(kpc)'.center(width) + '(km/s)'.center(width) + '(km/s)'.center(width) + 
        '(km/s)'.center(width) + '(cm^-3)'.center(width) + '(K)'.center(width) + 
        ' '.center(width) + ' '.center(width) + '(cm^-3)'.center(width) + 
        ' '.center(width) + '(cm^-3)'.center(width))
    f.write(s)
    f.write('\n')

    return f


def read_los_props(losnum):

    # Open filename
    f = open('lines.props')
    
    f.readline()

    los = losProps() 
    found = 0
    for line in f:
        l = line.split()
        if int(l[0]) == losnum:
           found = 1
           los.xen = float(l[3])
           los.yen = float(l[4])
           los.zen = float(l[5])
           los.losx = float(l[6])
           los.losy = float(l[7])
           los.losz = float(l[8])
           los.a11 = float(l[9])
           los.a12 = float(l[10])
           los.a13 = float(l[11])
           los.a21 = float(l[12])
           los.a22 = float(l[13])
           los.a23 = float(l[14])
           los.a31 = float(l[15])
           los.a32 = float(l[16])
           los.a33 = float(l[17])
           los.Xcom = float(l[18])
           los.Ycom = float(l[19])
           los.Zcom = float(l[20])
           los.VXcom = 0.0
           los.VYcom = 0.0
           los.VZcom = 0.0
           los.x0 = 0.0
           los.y0 = 0.0
           los.z0 = 0.0
           los.vx_obs = 0.0
           los.vy_obs = 0.0
           los.vz_obs = 0.0

    if found==0:
        print('Cannot find LOS {0:d} in lines.props'.format(losnum),flush=True)
        print('Function read_los_props in idcells.py'.format(losnum))
        print('Exitting'.format(losnum))
        sys.exit()
    return los
           
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################


def idcells(run, ions, codeLoc):

    """
    Function to select the physical properties of gas cells
    based on the cellID
    """

    ion_num = len(ions)

    width = 15

    # Read in LOS properties
    los_num, los_b, los_phi = read_lines('lines.info')

    # Loop over ions
    for ion in ions:

        #Read in ion box
        cwd = os.getcwd()
        ionboxfile = '{0:s}/{1:s}_GZa{2:s}.{3:s}.h5'.format(cwd,run.galID,
                                                            run.expn,ion.name)
        ionbox = pd.read_hdf(ionboxfile, 'data')

        # Loop over lines of sight
        for i in range(0,len(los_num)):
            
            # Read in LOS properties (entry points)
            los = read_los_props(i+1)

            # Construct filename that contains list of cells
            cell_file = 'los{0:04d}.cellID.dat'.format(i+1)
            cf = open(cell_file)
            cf.readline()     # Read past header

            # Construct output filename
            outfile = '{0:s}.{1:s}.los{2:04d}.dat'.format(run.galID,
                                                    ion.name,i+1)
            #outfile = galID+'.'+ion+'.los{0:04d}.dat'.format(i+1)

            # Write output header
            l = 0
            b = 90
            of = write_OutfileHdr(outfile,run.expn,los_b[i],los_phi[i],l,b,los)

            # Loop through cells in cellID file
            for line in cf:

                cellnum = int(line)
                ind = cellnum-1
                
                # Cell number corresponds to line of gas file
                cellsize = ionbox['cell_size'][ind]
                x = ionbox['x'][ind]
                y = ionbox['y'][ind]
                z = ionbox['z'][ind]
                vx = ionbox['vx'][ind]
                vy = ionbox['vy'][ind]
                vz = ionbox['vz'][ind]
                nH = ionbox['nH'][ind]
                t = ionbox['temperature'][ind]
                SNII_frac = ionbox['SNII'][ind]
                SNIa_frac = ionbox['SNIa'][ind]
                natom = ionbox['nAtom'][ind]
                fion = ionbox['fIon'][ind]
                nion= ionbox['nIon'][ind]
                cell_id = int(ionbox['ID'][ind])

                if cell_id != cellnum:
                    print('Error')
                    sys.exit()

                # Write to file
                s = (('{0:1.4e}'.format(cellsize)).center(width) + 
                    '{0:.4e}'.format(x).center(width) + '{0:.4e}'.format(y).center(width) + 
                    '{0:.4e}'.format(z).center(width) + '{0:.4e}'.format(vx).center(width) + 
                    '{0:.4e}'.format(vy).center(width ) + '{0:.4e}'.format(vz).center(width) + 
                    '{0:.4e}'.format(nH).center(width) + '{0:.4e}'.format(t).center(width) + 
                    '{0:.4e}'.format(SNII_frac).center(width) + 
                    '{0:.4e}'.format(SNIa_frac).center(width) + 
                    '{0:.4e}'.format(natom).center(width) + '{0:.4e}'.format(fion).center(width) + 
                    '{0:.4e}'.format(nion).center(width) +
                    '{0:d}\n'.format(cell_id).rjust(width-3))
                of.write(s)


            cf.close()
            of.close()

        del ionbox
    # Clean up the directory
    clean_up()



def clean_up():

    # Make a directory for the idcells files
    dirName = 'cellIDs'
    sp.call('mkdir '+dirName, shell=True)

    # Move all cellid files into the directory
    command = 'mv los*cellID.dat ./cellIDs/'
    try:
        sp.check_call(command, shell=True)
    except sp.CalledProcessError:
        pass

