#
#  mockspec.py
# 
#  Full pipeline for performing synthetic QSO observations
#  of ART simulations of galaxy CGM

# General libraries
from __future__ import print_function
#import numpy as np
#import tables as tb
import os
#import pandas as pd
import subprocess as sp
import sys
import joblib as jl
import time 
import datetime as dt

# Code specific libraries
import files as fi
import ratesControl as rc
import funcs.locatecells.locatecells as lc
import funcs.idcells.idcells as ic
import runFunctions as rf 
#import funcs.sigcells.sigcells as sc
import funcs.plotting.analysis_control as ac
import funcs.compileFiles.mkHDF5 as hdf
import funcs.tpcf.tpcf as cf



def ionLoop(run,ion):

    if sum([run.runLos7,run.runSpecsynth,run.runSysanal,
            run.runCullabs,run.runLocateCells])==0:
        print('\t\tNothing to do',flush=True)
        return 1

    ionloc = '{0:s}/i{1:d}/{2:s}'.format(run.runLoc,
                int(run.incline),ion.name)
    print('\n\n\tIon: ', ion.name,flush=True)
    # Move into the ion directory
    os.chdir(ionloc)

    ##### 
    #  
    #  los7
    #  Determine path length through each cell and create the
    #  .losdata and .lines files
    #
    #####
    if run.runLos7==1:
        print('\n\tDetermining cell path length and applying rough cut',
                flush=True)
        rf.los7(codeLoc)
    else:
        print('\tSkipping los7...',flush=True)


    ##### 
    #  
    #  specsynth
    #  Generate the synthetic spectra
    #
    #####
    if run.runSpecsynth==1:
        print('\n\tGenerating spectra',flush=True)
        rf.specsynth(codeLoc)
    else:
        print('\tSkipping specsynth...',flush=True)

    ##### 
    #  
    #  sysanal
    #  Analyze the spectra
    #
    #####
    if run.runSysanal==1:
        print('\n\tAnalyzing spectra...',flush=True)
        rf.sysanal(codeLoc)
    else:
        print('\tSkipping sysanal...',flush=True)


    ##### 
    #  
    #  cullabs
    #  Create sysabs file
    #
    #####
    if run.runCullabs==1:
        print('\n\tCreating sysabs',flush=True)
        rf.cullabs(codeLoc)
        hdf.sysabs_to_hdf5(run,ion,codeLoc)
        hdf.regabs_to_hdf5(run,ion,codeLoc)
    else:
        print('\tSkipping sysabs...',flush=True)

    ##### 
    #  
    #  locatecells
    #  Identify the cells that are significant contributers to the absorption
    #
    #####
    if run.runLocateCells==1:
        print('\n\tIdentifying significant cells',flush=True)
        lc.locateSigCells(run,ion,codeLoc)
        #lc.locateSigCells(galID, expn, ion, sigcellsCut, codeLoc, incline)
#        lc.sigCells(galID, expn, ion, sigcellsCut, codeLoc)
        
        # Conversion to HDF5 no longer needed as locatecells 
        # outputs an HDF5 file
        # hdf.abscells_to_hdf5(codeLoc,run,ion)
    else:
        print('\tSkipping locatecells...',flush=True)

    

    # Move back up to the parent directory
    os.chdir('..')




startTime = time.time()

# Get the location where the code lives
pathname = os.path.dirname(sys.argv[0])
codeLoc = os.path.abspath(pathname)
summaryLoc = codeLoc+'/summaries/'
requiredLoc = codeLoc+'/controls/'

#  Read in the control file
print('\n\nReading in control file...',flush=True)
run,ions,tpcfProp = fi.read_control_file()

# Genearte the name of the gasfile
gasfile = '{0:s}_GZa{1:s}.txt'.format(run.galID,run.expn)
print('\nRun Parameters:',flush=True)
print('\tGasfile:      ', gasfile,flush=True)
print('\tIons:         ', [i.name for i in ions],flush=True)
print('\tNLOS:         ', run.nlos,flush=True)
print('\tMax impact:   ', run.maximpact,flush=True)
print('\tIncline:      ', run.incline,flush=True)
print('\tEWCut:        ', run.ewcut,flush=True)
print('\tSNR:          ', run.snr,flush=True)
print('\tNCores:       ', run.ncores,flush=True)
print('\tRun Loc:      ', run.runLoc,flush=True)
print('\tSigcells Cut: {0:.0%}'.format( 
        run.sigcellsCut/100.),flush=True)

print('\nRun Flags:',flush=True)
print('\tRates:       {0:d}'.format(run.runRates,flush=True),flush=True)
print('\tGenLOS:      {0:d}'.format(run.runGenLOS,flush=True),flush=True)
print('\tCellfinder:  {0:d}'.format(run.runCellfinder,flush=True),flush=True)
print('\tIDCells:     {0:d}'.format(run.runIdcells,flush=True),flush=True)
print('\tLos7:        {0:d}'.format(run.runLos7,flush=True),flush=True)
print('\tSpecsynth:   {0:d}'.format(run.runSpecsynth,flush=True),flush=True)
print('\tSysanal:     {0:d}'.format(run.runSysanal,flush=True),flush=True)
print('\tCullabs:     {0:d}'.format(run.runCullabs,flush=True),flush=True)
print('\tLocateCells: {0:d}'.format(run.runLocateCells,flush=True),flush=True)
print('\tPlotting:    {0:d}'.format(run.runPlotting,flush=True),flush=True)
sys.stdout.flush()

# Setup timing file
timef,timeStr = fi.timing_setup(startTime,run,ions,tpcfProp)

# Test the summary location
sumFile = '{0:s}/rotmat_a{1:s}.txt'.format(os.getcwd(), run.expn)

# Generate gal_props.dat file
print('\n\nGenerating gal_props.dat...',flush=True)
#fi.setup_galprops( galID, expn, requiredLoc, summaryLoc )
fi.setup_galprops( run, requiredLoc, summaryLoc )

# Generate galaxy.props file, needed for analysis codes
print('\n\nGenerating galaxy.props...',flush=True)
#fi.setup_galaxy_props(sumFile, galID, expn, incline)
fi.setup_galaxy_props(run, sumFile)

##### 
#  
#  Run rates
#
#####
if run.runRates==1:
    print('\nRates:',flush=True)
    print('\t Setting up rates.inp...',flush=True)
    #rc.setup_rates_control( gasfile, run.expn, ions, requiredLoc)
    print('\t Setting up rates.outfiles...',flush=True)
    #rc.setup_rates_outputs(run, ions, codeLoc, requiredLoc) 
    print('\t Setting up rates data location...',flush=True)
    #rc.setup_rates_data(codeLoc)
    print('\t Running rates...',flush=True)
    #rc.run_rates(codeLoc)
    print('\t Converting to hdf5...',flush=True)
    hdf.gasbox_to_hdf5(codeLoc, ions, run)
else:
    print('Skipping rates',flush=True)

fi.record_time(timef,timeStr,'Rates',time.time(),startTime)
print('\n\nTime elapse = {0:s}\n'.format(str(dt.timedelta(
        seconds=time.time()-startTime))),flush=True)


##### 
#  
#  Set up inclination directory
#
#####
mainLoc = os.getcwd()
#incLoc = fi.setup_inclination_dir(incline, ions, runRates, galID, expn)
incLoc = fi.setup_inclination_dir(run, ions)
os.chdir(incLoc)

fi.record_time(timef,timeStr,'Inclination Dir',time.time(),startTime)


##### 
#  
#  Generate the lines of sight
#
#####
if run.runGenLOS==1:
    print('\nGenerating lines of sight...',flush=True)
    rf.genLOS( codeLoc, mainLoc, run )
else:
    print('Skipping genLOS...',flush=True)
fi.record_time(timef,timeStr,'genLOS',time.time(),startTime)
print('\n\nTime elapse = {0:s}\n'.format(str(dt.timedelta(
        seconds=time.time()-startTime))),flush=True)

##### 
#  
#  Run lines of sight through box with cellfinder
#
#####
if run.runCellfinder==1:
    print('\nRunning LOS through box...',flush=True)
    rf.runCellfinder(codeLoc, run.ncores)
else: 
    print('Skipping cellfinder...',flush=True)
fi.record_time(timef,timeStr,'cellfinder',time.time(),startTime)
print('\n\nTime elapse = {0:s}\n'.format(str(dt.timedelta(
        seconds=time.time()-startTime))),flush=True)

##### 
#  
#  Identify cells in the ion boxes
#
#####
if run.runIdcells==1:
    print('\nIdentifying probed cell properites...',flush=True)
    ic.idcells(run, ions, codeLoc)
else:
    print('Skipping IDcells...',flush=True)
fi.record_time(timef,timeStr,'idCells',time.time(),startTime)
print('\n\nTime elapse = {0:s}\n'.format(str(dt.timedelta(
        seconds=time.time()-startTime))),flush=True)

# Setup Mockspec 
fi.setup_mockspec(ions, run, requiredLoc)
fi.record_time(timef,timeStr,'setup mockspec',time.time(),startTime)
print('\n\nTime elapse = {0:s}\n'.format(str(dt.timedelta(
        seconds=time.time()-startTime))),flush=True)


# Setup the ion directory
#ionloc = fi.setup_ion_dir(ion, galID, expn, codeLoc) 
for ion in ions:
    ionloc = fi.setup_ion_dir(ion, run, codeLoc) 
fi.record_time(timef,timeStr,'ion dir',time.time(),startTime)
print('\n\nTime elapse = {0:s}\n'.format(str(dt.timedelta(
        seconds=time.time()-startTime))),flush=True)


# Start looping over ions
jl.Parallel(n_jobs=run.ncores,verbose=5)(
        jl.delayed(ionLoop)(run,ion) for ion in ions)
fi.record_time(timef,timeStr,'ion loop',time.time(),startTime)
print('\n\nTime elapse = {0:s}\n'.format(str(dt.timedelta(
        seconds=time.time()-startTime))),flush=True)
    
##### 
#  
#  Rename files
#  Add inclination angle to filenames
#
#####
# This is no longer needed as the cullabs functions output
# files with the correct names
#if run.runCullabs==1 or run.runLocateCells==1:
#    for ion in ions:
#        print('\n\t Renaming files...',flush=True)
#        fi.rename(run,ion)
#else:
#    print('What was the point...',flush=True)

# Generate summary files

# Summary is currently broken, fix later
run.runSummaries = 0
if run.runSummaries==1:
    print('\n\nGenerating summary files',flush=True)
    hdf.genSummaries(run, ions)    
fi.record_time(timef,timeStr,'summaries',time.time(),startTime)
print('\n\nTime elapse = {0:s}\n'.format(str(dt.timedelta(
        seconds=time.time()-startTime))),flush=True)

# Generate TPCFs
if run.runTPCF==1:
    print('\n\nGenerating TPCFs',flush=True)
    cf.absorber_tpcf(run,ions,tpcfProp)
fi.record_time(timef,timeStr,'tpcf',time.time(),startTime)
print('\n\nTime elapse = {0:s}\n'.format(str(dt.timedelta(
        seconds=time.time()-startTime))),flush=True)


# Create basic plots
if run.runPlotting==1:
    print('\n\nGenerating plots',flush=True)
    ac.make_plots(ions)
fi.record_time(timef,timeStr,'Done',time.time(),startTime)
timef.close()

