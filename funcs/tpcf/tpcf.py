
'''
Module to calcuate the TPCF between pixel velocities within absorbers 
defined by sysanal

Steps:
1) Identify each absorbing region through regabs/sysabs outputs
2) Select the velocity of each pixel in the regions
3) Calculate the velocity seperation between each possible pixel pairs
4) Bin the seperations into bins 10-20 km/s wide
5) Normalize the histograms
6) Bootstrap to get errors
'''

from __future__ import print_function,division
import pandas as pd
import numpy as np
import subprocess as sp
import itertools as it 

import sys

def absorber_tpcf(run,ion):

    '''
    Calculates TPCF within each absorber
    '''


    # Set up a dataframe with regabs information
    absorbers = regabs(run,ion)
    print(len(absorbers))
    absorbers.to_csv('absorbers.csv',index=False)

    # Set up a dataframe with pixel velocity information
    # Each column is a seperate absorber
    pixVel = velocities(run,ion,absorbers)
    
    # Get the seperation between each possible pair of 
    # pixel velocties
    velDiff = seperations(run,ion,pixVel)

def regabs(run,ion):

    '''
    Gets all information about the absorbing regions from regabs 
    and sysabs files
    Returns a dataframe
    '''

    loc = '{0:s}/i{1:d}/{2:s}/'.format(run.rootLoc,int(run.incline),ion.name)
    sysName = '{0:s}.{1:s}.los{2:04d}.sysabs'
    regName = '{0:s}.{1:s}.los{2:04d}.regabs'
    
    los,regnum,zabs,vneg,vpos,ew = [],[],[],[],[],[]

    for i in range(1,run.nlos):

        sName = loc+sysName.format(run.galID,ion.name,i)
        rName = loc+regName.format(run.galID,ion.name,i)

        try:
            with open(rName) as f:
                f.readline()
                for j,line in enumerate(f):
                    l = line.split()
                    los.append(i)
                    regnum.append(int(l[0]))
                    zabs.append(float(l[1]))
                    vneg.append(float(l[2]))
                    vpos.append(float(l[3]))
                    ew.append(float(l[4]))
        except IOError:
            try:
                 with open(sName) as f:
                    f.readline()
                    for j,line in enumerate(f):
                        l = line.split()
                        los.append(i)
                        regnum.append(0)
                        zabs.append(float(l[0]))
                        vneg.append(float(l[1]))
                        vpos.append(float(l[2]))
                        ew.append(float(l[3]))           
            except IOError:
                pass
                    

    header = 'los regnum zabs vneg vpos ew'.split()
    df = pd.DataFrame(index=range(len(los)),columns=header)
    fields = [los,regnum,zabs,vneg,vpos,ew]
    for col,field in zip(header,fields):
        df[col] = field
                
    return df
    

def velocities(run,ion,absorbers):

    '''
    Gets all the pixel velocities in each spectrum that lie within the 
    velocity limits of the region
    Returns a dataframe with each column corresponding to an absorbing
    region
    '''

    
    loc = '{0:s}/i{1:d}/{2:s}/'.format(run.rootLoc,int(run.incline),ion.name)

    specFile = '{0:s}.{1:s}.los{2:04d}.{3:s}.spec'
    specHeader = 'lambda velocity flux dum1 dum2 dum3'.split()

    # Get the transition name
    # Use the transition for the ion that is turned on with the strongest 
    # oscillator strength
    transitions = pd.read_csv(loc+'Mockspec.transitions',sep='\s+')
    transition = transitions.loc[transitions[transitions['ion']==ion.name]['fosc'].idxmax()]['trani']
    
    pixVel = pd.DataFrame()

    catch = 0
    system = 0
    for i in range(1,run.nlos+1):

        sfile = loc+specFile.format(run.galID,ion.name,i,transition)
        try:
            spec = pd.read_csv(sfile,sep='\s+',names=specHeader)
        except IOError:
            continue
        
        regions = absorbers[absorbers['los']==i]
        for j in range(len(regions)):
            vlo = regions['vneg'].iloc[j]
            vhi = regions['vpos'].iloc[j]
            selection = (spec['velocity']>=vlo) & (spec['velocity']<=vhi)
            vels = spec[selection]['velocity']
            vels.reset_index(inplace=True,drop=True)
            label = '{0:d}.{1:d}'.format(i,j)
            vels.name = label
            vels = pd.DataFrame(vels)
            
            pixVel = pd.concat([pixVel,vels],axis=1)
    
    pixVel.to_csv('pixVel.csv',index=False)
            
    return pixVel
    

def seperations(run,rion,pixVel):
    '''
    Calculates the velocity seperation for all possible combinations
    of pixel velocities in each column of pixVel
    Returns a dataframe were each column is an absorber containing all the
    seperations within the absorber
    '''

    velDiff = pd.DataFrame(columns=pixVel.columns)


    for col in pixVel.columns:
        
        velDiff[col] = pixVel[col].apply(pixSep)
        comb = it.combinations(pixVel[col],2)

def pix_sep(data):
    comb = np.array(it.combinations(data.dropna(),2))
    return np.diff(comb)
    

if __name__ == '__main__':

    from testclass import *

    run = runProps()
    run.galID = 'D9o2'
    run.expn = '1.002'
    run.rootLoc = '/mnt/cluster/abs/cgm/dwarfs/D9o2/a1.002/'
    run.rootLoc = '/Users/jacob/research/dwarfs/D9o2/a1.002/'
    run.incline = 90
    
    ion = ionProps()
    ion.name = 'MgII'

    print(ion.name)
    absorber_tpcf(run,ion)




