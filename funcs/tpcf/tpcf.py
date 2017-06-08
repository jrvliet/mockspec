
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
    print('Absorber Size = {0:d}'.format(len(absorbers)))

    # Set up a dataframe with pixel velocity information
    # Each column is a seperate absorber
    pixVel = velocities(run,ion,absorbers)
    pixVel.to_csv('pixVel.csv',index=False)
    print('PixVel Size = {0:d}'.format(len(pixVel)))
    
    # Get the seperation between each possible pair of 
    # pixel velocties
    velDiff = seperations(run,ion,pixVel)
    velDiff.to_csv('velDiff.csv',index=False)
    print('velDiff Size = {0:d}'.format(len(velDiff)))


    # Bin the data to create TPCF
    # Set up bins
    maxDiff = velDiff.fillna(0).values.max()
    print('Maxdiff = {0:.4f}'.format(maxDiff))
    binSize = 10.
    nbins = int(np.ceil(maxDiff/binSize))
    print('nbins = {0:d}'.format(nbins))
    endPoint = binSize*(nbins+1)
    print('endPoint = {0:.4f}'.format(endPoint))
    bins = np.arange(0,endPoint,binSize)
    labels = [(bins[i]+bins[i+1])/2. for i in range(nbins)]

    # Bin the velDiff dataframe
    velBins = velDiff.apply(lambda x: 
                pd.cut(x,bins,labels=labels,include_lowest=True))
    c = pd.Series(velBins.values.flatten()).value_counts().sort_index()
    
    fig,ax = plt.subplots(1,1,figsize=(7,7))
    #ax.scatter(c.index,c,marker='s',s=5)
    ax.plot(c.index,c,marker='s',ms=3)
    fig.savefig('tpcf_bin.png',bbox_inches='tight',dpi=300)
    plt.close(fig)

    # Plot
    fig,ax = plt.subplots(1,1,figsize=(7,7))
    ax.hist(pd.Series(velDiff.values.flatten()).dropna(),bins=50,histtype='step')
    fig.savefig('tpcf_hist.png',bbox_inches='tight',dpi=300)
    plt.close(fig)


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
    
            
    return pixVel
    

def seperations(run,rion,pixVel):
    '''
    Calculates the velocity seperation for all possible combinations
    of pixel velocities in each column of pixVel
    Returns a dataframe were each column is an absorber containing all the
    seperations within the absorber
    '''

    velDiff = pd.DataFrame(columns=pixVel.columns)

    print(len(pixVel.columns))
    for col in pixVel.columns:
        
        # Create pairs of all possible combinations
        comb = np.array(list(it.combinations(pixVel[col].dropna(),2)))

        # Get the difference between the pairs
        diff = np.abs(np.diff(comb))
        diff = pd.Series(diff.flatten(),name=col)
        velDiff = pd.concat([velDiff,diff],axis=1,ignore_index=True)

    return velDiff






if __name__ == '__main__':

    from testclass import *
    import matplotlib.pyplot as plt

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




