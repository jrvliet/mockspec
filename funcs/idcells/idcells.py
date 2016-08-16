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

    
def write_OutfileHdr(outfile, aexpn, R0, phi, l, b, xen, yen, zen, losx, losy, losz, a11, a12, a13, a21, a22, a23, a31, a32, a33, Xcom, Ycom, Zcom, VXcom, VYcom, VZcom, x0, y0, z0, vx_obs, vy_obs, vz_obs):

    # Write header information to output file

    from math import sqrt

    f = open(outfile,'w')
    
    str = 'aexpn   {0:.3f} Galaxy center in kpc:  {1:.6f}  {2:.6f}  {3:.6f} Peculiar velocity center in km/s:    {4:.3f}   {5:.3f}    {6:.3f}\n'.format(float(aexpn), Xcom, Ycom, Zcom, VXcom, VYcom, VZcom)
    f.write(str)
    
    str = 'Galactic Coordinates (l:b) in degrees:     {0:.1f}    {1:.1f}  Observer Position (R(kpc):phi(deg)):     {2:.2f}  {3:.2f}\n'.format(l,b,R0,phi)
    f.write(str)
    
    str = 'Observer Velocity (vx vy vz) in km/s in Box Frame:        {0:.3E}      {1:.3E}      {2:.3E}'.format(vx_obs, vy_obs, vz_obs)
    f.write(str)
    f.write('\n')

    str = 'Box Entry Point (x0:y0:z0) in kpc:      {0:.6f}     {1:.6f}       {2:.6f}  LOS Direction Cosines (l:m:n):    {3:.6E}   {4:.6E}  {5:.6E}'.format(xen,yen,zen,losx,losy,losz)
    f.write(str)
    f.write('\n')

    str = 'Z prime axis and Lvec/L (a31 a32 a33): {0: .7E} {1: .7E} {2: .7E}\n'.format(a31,a32,a33)
    f.write(str)

    str = 'X prime axis            (a11 a12 a13): {0: .7E} {1: .7E} {2: .7E}\n'.format(a11,a12,a13)
    f.write(str)

    str = 'Y prime axis            (a21 a22 a23): {0: .7E} {1: .7E} {2: .7E}\n'.format(a21,a22,a23)
    f.write(str)
    
    width = 15

    str = '1'.center(width) + '2'.center(width) + '3'.center(width) + '4'.center(width) + '5'.center(width) + '6'.center(width) + '7'.center(width) + '8'.center(width) + '9'.center(width) + '10'.center(width) + '11'.center(width) + '12'.center(width) + '13'.center(width) + '14'.center(width) + '15'.center(width)+ '\n'
    f.write(str)

    str = '{0:^15}'.format('cell size') + '{0:^15}'.format('x') + '{0:^15}'.format('y') + '{0:^15}'.format('z') + '{0:^15}'.format('vx') + '{0:^15}'.format('vy') + '{0:^15}'.format('vz') + '{0:^15}'.format('nH') + '{0:^15}'.format('temperature') + '{0:^15}'.format('SNII mass frac') + '{0:^15}'.format('SNIa mass frac') + '{0:^15}'.format('natom') + '{0:^15}'.format('fion') + '{0:^15}'.format('nion') + '{0:^15}\n'.format('cell id')

    f.write(str)
   
    str = '(kpc)'.center(width) + '(kpc)'.center(width) + '(kpc)'.center(width) + '(kpc)'.center(width) + '(km/s)'.center(width) + '(km/s)'.center(width) + '(km/s)'.center(width) + '(cm^-3)'.center(width) + '(K)'.center(width) + ' '.center(width) + ' '.center(width) + '(cm^-3)'.center(width) + ' '.center(width) + '(cm^-3)'.center(width)
    f.write(str)
    f.write('\n')

    return f


def read_los_props(losnum):

    # Open filename
    f = open('lines.props')
    
    f.readline()

    for line in f:
        l = line.split()
        if int(l[0]) == losnum:
           xen = float(l[3])
           yen = float(l[4])
           zen = float(l[5])
           losx = float(l[6])
           losy = float(l[7])
           losz = float(l[8])
           a11 = float(l[9])
           a12 = float(l[10])
           a13 = float(l[11])
           a21 = float(l[12])
           a22 = float(l[13])
           a23 = float(l[14])
           a31 = float(l[15])
           a32 = float(l[16])
           a33 = float(l[17])
           Xcom = float(l[18])
           Ycom = float(l[19])
           Zcom = float(l[20])
           VXcom = 0.0
           VYcom = 0.0
           VZcom = 0.0
           x0 = 0.0
           y0 = 0.0
           z0 = 0.0
           vx_obs = 0.0
           vy_obs = 0.0
           vz_obs = 0.0

    return xen, yen, zen, losx, losy, losz, a11, a12, a13, a21, a22, a23, a31, a32, a33, Xcom, Ycom, Zcom, VXcom, VYcom, VZcom, x0, y0, z0, vx_obs, vy_obs, vz_obs
           


###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################


def idcells(galID, aexpn, ion_list, codeLoc):

    """
    Function to select the physical properties of gas cells
    based on the cellID
    """

    ion_num = len(ion_list)

    width = 15

    # Read in LOS properties
    los_num, los_b, los_phi = read_lines('lines.info')

    # Loop over ions
    for ion in ion_list:

        #Read in ion box
        cwd = os.getcwd()
        ionboxfile = '{0:s}/{1:s}_GZa{2:s}.{3:s}.h5'.format(cwd,galID,aexpn,ion)
        ionbox = pd.read_hdf(ionboxfile, 'data')
        print ionboxfile
        print ionbox.shape

        # Loop over lines of sight
        for i in range(0,len(los_num)):
            
            # Read in LOS properties (entry points)
            xen, yen, zen, losx, losy, losz, a11, a12, a13, a21, a22, a23, a31, a32, a33, Xcom, Ycom, Zcom, VXcom, VYcom, VZcom, x0, y0, z0, vx_obs, vy_obs, vz_obs = read_los_props(i+1)

            # Construct filename that contains list of cells
            cell_file = 'los{0:04d}.cellID.dat'.format(i+1)
            cf = open(cell_file)
            cf.readline()     # Read past header

            # Construct output filename
            outfile = galID+'.'+ion+'.los{0:04d}.dat'.format(i+1)

            # Write output header
            l = 0
            b = 90
            of = write_OutfileHdr(outfile, aexpn, los_b[i], los_phi[i], l, b, xen, yen, zen, losx, losy, losz, a11, a12, a13, a21, a22, a23, a31, a32, a33, Xcom, Ycom, Zcom, VXcom, VYcom, VZcom, x0, y0, z0, vx_obs, vy_obs, vz_obs)

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
                    print 'Error'
                    sys.exit()

                # Write to file
                str = ('{0:1.4e}'.format(cellsize)).center(width) + '{0:.4e}'.format(x).center(width) + '{0:.4e}'.format(y).center(width) + '{0:.4e}'.format(z).center(width) + '{0:.4e}'.format(vx).center(width) + '{0:.4e}'.format(vy).center(width ) + '{0:.4e}'.format(vz).center(width) + '{0:.4e}'.format(nH).center(width) + '{0:.4e}'.format(t).center(width) + '{0:.4e}'.format(SNII_frac).center(width) + '{0:.4e}'.format(SNIa_frac).center(width) + '{0:.4e}'.format(natom).center(width) + '{0:.4e}'.format(fion).center(width) + '{0:.4e}'.format(nion).center(width) + '{0:d}\n'.format(cell_id).rjust(width-3)
                of.write(str)


            cf.close()
            of.close()

    # Clean up the directory
    clean_up()



def clean_up():

    # Make a directory for the idcells files
    dirName = 'cellIDs'
    sp.call('mkdir '+dirName, shell=True)

    # Move all cellid files into the directory
    command = 'mv los*cellID.dat ./cellIDs/'
    sp.call(command, shell=True)


