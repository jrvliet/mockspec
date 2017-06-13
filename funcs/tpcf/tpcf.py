
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

class tpcfProps(object):

    def __init__ (self):
       
        self.ewLo = 0.
        self.ewHi = 10.
        self.dLo = 0.
        self.dHi = 200.
        self.fraction = 0.1
        self.binSize = 10.
        self.bootNum = 10.

 
#def absorber_tpcf(run,ion):
def absorber_tpcf(run,ion,tpcfProp):

    '''
    Calculates TPCF within each absorber
    '''

    # Set up a dataframe with regabs information
    print('Reading Absorbers')
    absorbers = regabs(run,ion,tpcfProp)
    absorbersName = '{0:s}_{1:s}_{2:s}_absorbers.csv'.format(
                    run.galID,run.expn,ion.name)
    absorbers.to_csv(absorbersName,index=False)

    # Set up a dataframe with pixel velocity information
    # Each column is a seperate absorber
    print('Creating pixel velocity')
    pixVel = velocities(run,ion,absorbers)
    pixVelName = '{0:s}_{1:s}_{2:s}_pixVel.csv'.format(
                    run.galID,run.expn,ion.name)
    pixVel.to_csv(pixVelName,index=False)
    
    # Get the seperation between each possible pair of 
    # pixel velocties
    print('Calculating velocity seperations')
    velDiff = seperations(run,ion,pixVel)
    velDiffName = '{0:s}_{1:s}_{2:s}_velDiff.csv'.format(
                    run.galID,run.expn,ion.name)
    velDiff.to_csv(velDiffName,index=False)

    # Bin the data to create TPCF
    # Set up bins
    maxDiff = velDiff.fillna(0).values.max()
    binSize = tpcfProp.binSize
    nbins = int(np.ceil(maxDiff/binSize))
    endPoint = binSize*(nbins+1)
    bins = np.arange(0,endPoint,binSize)
    labels = [(bins[i]+bins[i+1])/2. for i in range(nbins)]

    # Bin the velDiff dataframe
    print('Binning')
    velBins = velDiff.apply(lambda x: 
                pd.cut(x,bins,labels=labels,include_lowest=True))
    c = pd.Series(velBins.values.flatten()).value_counts().sort_index()
    c = c/c.sum()
    
    # Bootstrap for errors
    print('Bootstrapping')
    cMean,cStd = bootstrap(bins,labels,tpcfProp,velDiff)

    return c,cMean,cStd

def bootstrap(bins,labels,tpcfProp,velDiff):
    
    '''
    Calculates errors using bootstrap, randomly selecting 
    fraction of total lines from velDiff a number of times
    '''

    # Create dataframe to hold all resulting TPCFs
    # Each row will be a TPCF from a sample of absorbers
    df = pd.DataFrame()
    
    # Loop over all bootstrap instances
    for i in range(tpcfProp.bootNum):

        # Pick out a sample of absorbers from velDiff
        sample = velDiff.sample(frac=tpcfProp.fraction,axis=1)

        # Bin the samples
        velBins = sample.apply(lambda x: 
                    pd.cut(x,bins,labels=labels,include_lowest=True))
        c = pd.Series(velBins.values.flatten()).value_counts().sort_index()
        c = c/c.sum()
            
        df = pd.concat([df,c],ignore_index=True,axis=1)
    
    df.to_csv('bootstrap.csv',index=False)
    return df.mean(axis=1),df.std(axis=1)
        

    

def cutBins(x,bins,labels):

    # Add noise to bin edges
    noiseMu = 0.
    noiseSigma = 1.2
    noise = np.random.normal(noiseMu,noiseSigma,len(bins))
    bins = bins+noise
    
    return pd.cut(x,bins,labels=labels,include_lowest=True)
    
    

def regabs(run,ion,tpcfProp):
    '''
    Gets all information about the absorbing regions from regabs 
    and sysabs files
    Returns a dataframe
    '''
    
    absheader = 'los D phi region zabs v- v+ EW_r'.split()

    fname = '{0:s}/i{1:d}/{2:s}/{3:s}.{2:s}.a{4:s}.i{1:d}.ALL.regabs.h5'.format(
            run.rootLoc,int(run.incline),ion.name,run.galID,run.expn)
    try:
        regdf = pd.read_hdf(fname,'data')
        reglos = regdf['los'].unique()

        regdf = regdf[absheader]

    except IOError:
        regdf = pd.DataFrame(columns=absheader)
        reglos = []
    
    fname = '{0:s}/i{1:d}/{2:s}/{3:s}.{2:s}.a{4:s}.i{1:d}.ALL.sysabs.h5'.format(
            run.rootLoc,int(run.incline),ion.name,run.galID,run.expn)
    try:
        alldf = pd.read_hdf(fname,'data')
    except IOError:
        print('Cannot open {0:s} in regabs in tpcf.py'.format(fname))
        sys.exit()
    alldf['region'] = 1.

    
    
    selection = ((alldf['EW_r']>=tpcfProp.ewLo) & 
                 (alldf['EW_r']<=tpcfProp.ewHi) &
                 (alldf['EW_r']>0) &
                 (alldf['D']>=tpcfProp.dLo) & 
                 (alldf['D']<=tpcfProp.dHi) &
                 (~alldf['los'].isin(reglos)))
    df = pd.concat([regdf,alldf[absheader][selection]],ignore_index=True)

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
    noiseMu = 0
    noiseSigma = 1.2
    for i in range(1,run.nlos+1):

        sfile = loc+specFile.format(run.galID,ion.name,i,transition)
        try:
            spec = pd.read_csv(sfile,sep='\s+',names=specHeader)
        except IOError:
            continue
        
        regions = absorbers[absorbers['los']==i]
        for j in range(len(regions)):
            vlo = regions['v-'].iloc[j]
            vhi = regions['v+'].iloc[j]
            selection = (spec['velocity']>=vlo) & (spec['velocity']<=vhi)
            vels = spec[selection]['velocity']
            vels.reset_index(inplace=True,drop=True)
            label = '{0:d}.{1:d}'.format(i,j)
            vels.name = label

            # Add random noise to the velocities
            vels = vels+np.random.normal(noiseMu,noiseSigma,len(vels))

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
    run.rootLoc = '/Users/jacob/research/dwarfs/D9o2/a1.002/'
    run.rootLoc = '/mnt/cluster/abs/cgm/dwarfs/D9o2/a1.002/'
    run.incline = 90
    

    tpcfProp = tpcfProps()
    tpcfProp.ewLo = 0.
    tpcfProp.ewHi = 5.
    tpcfProp.dLo = 5.0
    tpcfProp.dHi = 200.
    tpcfProp.fraction = 0.15
    tpcfProp.binSize = 10.
    tpcfProp.bootNum = 1000

    ionP = ionProps()
    ions = 'HI MgII CIV OVI'.split()
    fig,axes = plt.subplots(2,2,figsize=(10,10))
    
    for ion,ax in zip(ions,axes.flatten()):
        print('\nIon = {0:s}'.format(ion))
        ionP.name = ion

        c,cMean,cStd = absorber_tpcf(run,ionP,tpcfProp)

        ax.plot(c.index,c,marker='s',color='b',label='Full')
        ax.plot(c.index,cMean,marker='o',color='g',label='Mean')
        ax.plot(c.index,cMean-cStd,color='r',label='Err Down')
        ax.plot(c.index,cMean+cStd,color='r',label='Err Up')
        ax.legend(loc='upper right',frameon=True)
        
        ax.set_title(ion)
    
    fig.tight_layout()
    fig.savefig('tpcf_bin.png',bbox_inches='tight',dpi=300)
    plt.close(fig)

       



