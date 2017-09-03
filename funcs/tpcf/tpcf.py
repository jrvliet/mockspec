
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
import joblib as jl
import os
import numba as nb
import tempfile

import sys
import decimal


class tpcfProps(object):
    
    '''
    Class to decribe settings for TPCF
    '''

    def __init__ (self):

        self.ewLo = 0.
        self.ewHi = 10.
        self.dLo = 0.
        self.dHi = 200.
        self.azLo = 0.
        self.azHi = 90.
        self.binSize = 10.
        self.bootNum = 1000


def dfSize(df):

    size = df.values.nbytes + df.index.nbytes + df.columns.nbytes
    x = decimal.Decimal(size)
    return x.normalize().to_eng_string()

def absorber_tpcf(run,ions,tpcfProp):
    
    print('\n\n',flush=True)
    # Check if TPCF directory exits
    tpcfDir = '{0:s}/i{1:d}/tpcf/'.format(run.runLoc,int(run.incline))
    if not os.path.isdir(tpcfDir):
        command = 'mkdir {0:s}'.format(tpcfDir)
        sp.call(command,shell=True)
    
    for ion in ions:
        print('\nIon = {0:s}'.format(ion.name),flush=True)
        tpcf_ion_loop(run,ion,tpcfProp,tpcfDir)
    print('\n\nDone with absorber_tpcf function',flush=True)
    #jl.Parallel(n_jobs=run.ncores,verbose=5)(
    #    jl.delayed(tpcf_ion_loop)(run,ion,tpcfProp,tpcfDir) for ion in ions)

def tpcf_ion_loop(run,ion,tpcfProp,tpcfDir):

    '''
    Calculates TPCF within each absorber
    '''

    # Set up a dataframe with regabs information
    print('Reading Absorbers',flush=True)
    absorbers = regabs(run,ion,tpcfProp)
    absorbersName = '{0:s}/{1:s}_{2:s}_{3:s}_absorbers.csv'.format(
                    tpcfDir,run.galID,run.expn,ion.name)
    absorbers.to_csv(absorbersName,index=False)
    print('Absorbers Size = {0:s}'.format(dfSize(absorbers)),flush=True)

    # Set up a dataframe with pixel velocity information
    # Each column is a seperate absorber
    print('Creating pixel velocity',flush=True)
    pixVel = velocities(run,ion,absorbers)
    pixVelName = '{0:s}/{1:s}_{2:s}_{3:s}_pixVel.csv'.format(
                    tpcfDir,run.galID,run.expn,ion.name)
    pixVel.to_csv(pixVelName,index=False)
    print('Pixvel Size = {0:s}'.format(dfSize(pixVel)),flush=True)
    
    # Get the seperation between each possible pair of 
    # pixel velocties
    print('Calculating velocity seperations',flush=True)
    velDiffMem,velDiffShape = seperations(run,ion,pixVel,tpcfDir)

    # Bin the data to create TPCF
    # Set up bins
    print(velDiffMem,velDiffShape,flush=True)
    bins,labels = velocity_bins(velDiffMem,velDiffShape,tpcfProp)
       # Bin the velDiff dataframe
    print('Binning',flush=True)
    c = cut_bins(velDiffMem,velDiffShape,bins,labels)
    
    # Bootstrap for errors
    #boot = pd.DataFrame(index=range(tpcfProp.bootNum))
    #boot = np.empty((tpcfProp.bootNum,len(c)))
    #boot[:] = np.nan

    print('Bootstrapping',flush=True)
    path = tempfile.mkdtemp()
    mempath = os.path.join(path,'bootstrap.mmap')
    boot = np.memmap(mempath,dtype='float',
                    shape=(tpcfProp.bootNum,len(c)),mode='w+')
    jl.Parallel(n_jobs=run.ncores,verbose=5)(
        jl.delayed(bstrap)(velDiffMem,velDiffShape,bins,labels,boot,i)
        for i in range(tpcfProp.bootNum))
    
    # Generate the dataframe containing the final results
    #tpcfFull = pd.DataFrame(index=c.index)
    tpcfFull = pd.DataFrame(index=labels)
    tpcfFull.index.name = 'Velocity'
    print('\n\nShape of tpcfFull = ',tpcfFull.shape,flush=True)
    print('Length of c = ',len(c),flush=True)
    tpcfFull['Full'] = c
    m = np.nanmean(boot,axis=0)
    s = np.nanstd(boot,axis=0)
    print('Mean = ',m,flush=True)
    print('Std = ',s,flush=True)
    tpcfFull['Mean'] = np.pad(m,(0,len(c)-len(m)),'constant')
    tpcfFull['Std'] = np.pad(s,(0,len(c)-len(s)),'constant')

    tpcfFull[tpcfFull==0] = np.nan
    outName = '{0:s}/{1:s}_{2:s}_{3:s}_tpcf.csv'.format(
                tpcfDir,run.galID,run.expn,ion.name)
    tpcfFull.to_csv(outName)
    print('Finished ion loop for {0:s}'.format(ion.name),flush=True)

def velocity_bins(velDiffName,velDiffShape,tpcfProp):
    print('velDiffName = ',velDiffName,flush=True)
    print('VelDiffShape = ',velDiffShape,flush=True)
    velDiff = np.memmap(velDiffName,dtype='float',
                        mode='r',shape=velDiffShape)
    #maxDiff = velDiff.fillna(0).values.max()
    #maxDiff = velDiff.max()
    maxDiff = np.nanmax(velDiff)
    print('maxDiff = ',maxDiff)
    binSize = tpcfProp.binSize
    nbins = int(np.ceil(maxDiff/binSize))
    endPoint = binSize*(nbins+1)
    bins = np.arange(0,endPoint,binSize)
    labels = [(bins[i]+bins[i+1])/2. for i in range(nbins)]
    lastLabel = labels[-1]+(labels[1]-labels[0])
    labels.append(lastLabel)
    print('\n\nBins : \n\tNumber = {0:d}\n\tLen(bins) = {1:d}\n\tLen(labels) = {2:d}'.format(
            nbins,len(bins),len(labels)),flush=True)
    
    return bins,labels



@nb.jit
def bstrap(velDiffMem,velDiffShape,bins,labels,boot,i):
    
    #sample = velDiff.sample(frac=1.,replace=True)
    resample = 1
    b = cut_bins(velDiffMem,velDiffShape,bins,labels,resample)
    boot[i,:len(b)] = b



@nb.jit
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
        

def cut_bins(velDiffMem,velDiffShape,bins,labels,resample=0):
    sample = np.memmap(velDiffMem,dtype='float',
                        mode='r',shape=velDiffShape)
    if resample!=0:
        sample = sample[:,np.random.choice(velDiffShape[1],
                                           velDiffShape[1],
                                           replace=True)]
        
    flat = sample.flatten()
    flat = flat[~np.isnan(flat)]
    c = np.sort(np.bincount(np.digitize(flat,bins)))[::-1]
    #velBins = sample.apply(lambda x: 
    #            pd.cut(x,bins,labels=labels,include_lowest=True))
    #c = pd.Series(velBins.values.flatten()).value_counts().sort_index()
    c = c/c.sum()
    return c

    
def az(phi):
    if phi<90:
        return phi
    elif phi<180:
        return 180.-phi
    elif phi<270:
        return phi-180.
    else:
        return 360.-phi 

def regabs(run,ion,tpcfProp):
    '''
    Gets all information about the absorbing regions from regabs 
    and sysabs files
    Returns a dataframe
    '''
    
    absheader = 'los D phi region zabs v- v+ EW_r'.split()

    # Try to read in regabs ALL file
    fname = '{0:s}/i{1:d}/{2:s}/{3:s}.{2:s}.a{4:s}.i{1:d}.ALL.regabs.h5'.format(
            run.runLoc,int(run.incline),ion.name,run.galID,run.expn)
    try:
        regabs = pd.read_hdf(fname,'data')
        reglos = regabs['los'].unique()
        regabs = regdf[absheader]
        if 'azimuthal' not in regabs.columns:
            regabs['azimuthal'] = regabs['phi'].apply(az)
        noregabs = 0
    except IOError:
        print('No ALL.regabs file for {0:s}'.format(ion.name),flush=True)
        regabs = pd.DataFrame(columns=absheader)
        regabs = []
        noregabs = 1
    
    # Try to read in sysabs ALL file
    fname = '{0:s}/i{1:d}/{2:s}/{3:s}.{2:s}.a{4:s}.i{1:d}.ALL.sysabs.h5'.format(
            run.runLoc,int(run.incline),ion.name,run.galID,run.expn)
    try:
        sysabs = pd.read_hdf(fname,'data')
        sysabs['region'] = 1.
        if 'azimuthal' not in sysabs.columns:
            sysabs['azimuthal'] = sysabs['phi'].apply(az)
    except IOError:
        print('Cannot open {0:s} in regabs in tpcf.py'.format(fname),flush=True)
        sys.exit()

    
    if noregabs==0:
        nonRegSysabs = sysabs[~sysabs['los'].isin(regabs['los'])]
        full = pd.concat([nonRegSysabs,regabs],ignore_index=True)
        full = full.sort_values('los')
        full.reset_index(drop=True,inplace=True)
    else:
        full = sysabs

    selection = ((full['EW_r']>=tpcfProp.ewLo) & 
                 (full['EW_r']<=tpcfProp.ewHi) &
                 (full['EW_r']>0) &
                 (full['D']>=tpcfProp.dLo) & 
                 (full['D']<=tpcfProp.dHi) &
                 (full['azimuthal']>=tpcfProp.azLo) & 
                 (full['azimuthal']<=tpcfProp.azHi) &
                 (~full['los'].isin(reglos)))
    print('For ion {0:s}, len(alldf) = {1:d}, len(selection) = {2:d}'.format(
            ion.name,len(alldf),selection.sum()),flush=True)
    #df = pd.concat([regdf,alldf[absheader][selection]],ignore_index=True)

    return full[absheader][selection]

def velocities(run,ion,absorbers):

    '''
    Gets all the pixel velocities in each spectrum that lie within the 
    velocity limits of the region
    Returns a dataframe with each column corresponding to an absorbing
    region
    '''

    
    loc = '{0:s}/i{1:d}/{2:s}/'.format(run.runLoc,int(run.incline),ion.name)
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
    

def seperations(run,ion,pixVel,tpcfDir):
    '''
    Calculates the velocity seperation for all possible combinations
    of pixel velocities in each column of pixVel
    Returns a dataframe were each column is an absorber containing all the
    seperations within the absorber
    '''

    velDiff = pd.DataFrame(columns=pixVel.columns)

    for col in pixVel.columns:
        
        # Create pairs of all possible combinations
        comb = np.array(list(it.combinations(pixVel[col].dropna(),2)))

        # Get the difference between the pairs
        diff = np.abs(np.diff(comb))
        diff = pd.Series(diff.flatten(),name=col)
        #velDiff = pd.concat([velDiff,diff],axis=1,ignore_index=True)
        velDiff[col] = diff

    # Convert the velDiff array to memmap object
    velDiffPath = tempfile.mkdtemp()
    velMemPath = os.path.join(velDiffPath,'velDiff.mmap')
    velDiffMem = np.memmap(velMemPath,dtype='float',
                shape=velDiff.shape,mode='w+')
    print('velDiffMem Location: ',velMemPath,flush=True)
    velDiffMem[:] = velDiff.values[:]

    velDiffName = '{0:s}/{1:s}_{2:s}_{3:s}_velDiff.csv'.format(
                    tpcfDir,run.galID,run.expn,ion.name)
    velDiff.to_csv(velDiffName,index=False)
    print('Veldiff Size = {0:s}'.format(dfSize(velDiff)),flush=True)
    return velMemPath,velDiff.shape




if __name__ == '__main__':

    from testclass import *
    from files import *
    import matplotlib.pyplot as plt

    run = runProps()
    run.galID = 'vela2b-25'
    run.expn = '0.490'
    run.runLoc = '/Users/jacob/research/dwarfs/D9o2/a1.002/'
    run.runLoc = '/mnt/cluster/abs/cgm/vela2b/vela25/a0.490/'
    #run.runLoc = '/home/sims/vela2b/vela25/a0.490/'
    run.incline = 90

    tpcfProp = tpcfProps()
    tpcfProp.ewLo = 0.
    tpcfProp.ewHi = 10.
    tpcfProp.dLo = 5.0
    tpcfProp.dHi = 200.
    tpcfProp.fraction = 0.15
    tpcfProp.binSize = 20.
    tpcfProp.bootNum = 1000

    ions = 'HI MgII CIV OVI'.split()
    ionPs  = []
    for ion in ions:
        ionP = ionProps()
        ionP.name = ion
        ionPs.append(ionP)
    
    absorber_tpcf(run,ionPs,tpcfProp)
    #c,cMean,cStd = absorber_tpcf(run,ionPs,tpcfProp)

#    fig,axes = plt.subplots(2,2,figsize=(10,10))
#    for ion,ax in zip(ions,axes.flatten()):
#        print('\nIon = {0:s}'.format(ion),flush=True)
#        ionP.name = ion
#
#
#        ax.step(c.index,c,marker='s',color='b',label='Full',where='mid')
#        ax.step(c.index,cMean,marker='o',color='g',label='Mean',where='mid')
#        ax.step(c.index,cMean-cStd,color='r',label='Err Down',where='mid')
#        ax.step(c.index,cMean+cStd,color='r',label='Err Up',where='mid')
#        ax.legend(loc='upper right',frameon=True)
#        
#        ax.set_title(ion)
    
#    fig.tight_layout()
#    fig.savefig('tpcf_bin.png',bbox_inches='tight',dpi=300)
#    plt.close(fig)

       



