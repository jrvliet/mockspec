#
#  mockspec.py
# 
#  Full pipeline for performing synthetic QSO observations
#  of ART simulations of galaxy CGM

# General libraries
#import numpy as np
#import tables as tb
import os
#import pandas as pd
import subprocess as sp
import sys
import joblib as jl

# Code specific libraries
import files as fi
import ratesControl as rc
import funcs.locatecells.locatecells as lc
import funcs.idcells.idcells as ic
import runFunctions as rf 
import funcs.sigcells.sigcells as sc
import funcs.plotting.analysis_control as ac
import funcs.compileFiles.mkHDF5 as hdf



def ionLoop(run,ion):

    if sum([run.runLos7,run.runSpecsynth,run.runSysanal,run.runCullabs,run.runLocateCells])==0:
        print '\t\tNothing to do'
        return 1

    ionloc = '{0:s}/i{1:d}/{2:s}'.format(run.runLoc,int(run.incline),ion.name)
    print '\n\n\tIon: ', ion.name
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
        print '\n\tDetermining cell path length and applying rough cut'
        rf.los7(codeLoc)
    else:
        print '\tSkipping los7...'


    ##### 
    #  
    #  specsynth
    #  Generate the synthetic spectra
    #
    #####
    if run.runSpecsynth==1:
        print '\n\tGenerating spectra'
        rf.specsynth(codeLoc)
    else:
        print '\tSkipping specsynth...'

    ##### 
    #  
    #  sysanal
    #  Analyze the spectra
    #
    #####
    if run.runSysanal==1:
        print '\n\tAnalyzing spectra...'
        rf.sysanal(codeLoc)
    else:
        print '\tSkipping sysanal...'


    ##### 
    #  
    #  cullabs
    #  Create sysabs file
    #
    #####
    if run.runCullabs==1:
        print '\n\tCreating sysabs'
        rf.cullabs(codeLoc)
        hdf.sysabs_to_hdf5(run,ion,codeLoc)
        hdf.regabs_to_hdf5(run,ion,codeLoc)
    else:
        print '\tSkipping sysabs...'

    ##### 
    #  
    #  locatecells
    #  Identify the cells that are significant contributers to the absorption
    #
    #####
    if run.runLocateCells==1:
        print '\n\tIdentifying significant cells'
        lc.locateSigCells(run,ion,codeLoc)
        #lc.locateSigCells(galID, expn, ion, sigcellsCut, codeLoc, incline)
#        lc.sigCells(galID, expn, ion, sigcellsCut, codeLoc)
#        hdf.abscells_to_hdf5(codeLoc)
    else:
        print '\tSkipping locatecells...'


    # Move back up to the parent directory
    os.chdir('..')





# Get the location where the code lives
pathname = os.path.dirname(sys.argv[0])
codeLoc = os.path.abspath(pathname)
summaryLoc = codeLoc+'/summaries/'
requiredLoc = codeLoc+'/controls/'

#  Read in the control file
print '\n\nReading in control file...'
run, ions = fi.read_control_file()
#props, flags, ions, xh, instruments = fi.read_control_file()
#galID, expn, nlos, maximpact, incline, ewcut, snr, ncores, rootLoc, sigcellsCut = props
#runRates, runGenLOS, runCellfinder, runIdcells, runLos7, runSpecsynth, runSysanal, runCullabs, runLocateCells, runSummaries, runPlotting = flags

# Genearte the name of the gasfile
#gasfile = galID+'_GZa'+expn+'.txt'
gasfile = '{0:s}_GZa{1:s}.txt'.format(run.galID,run.expn)
print '\nRun Parameters:'
print '\tGasfile:      ', gasfile
print '\tIons:         ', [i.name for i in ions]
print '\tNLOS:         ', run.nlos
print '\tMax impact:   ', run.maximpact
print '\tIncline:      ', run.incline
print '\tEWCut:        ', run.ewcut
print '\tSNR:          ', run.snr
print '\tNCores:       ', run.ncores
print '\tRun Loc:      ', run.runLoc
print '\tSigcells Cut: {0:.0%}'.format( run.sigcellsCut/100. )

print '\nRun Flags:'
print '\tRates:       {0:d}'.format(run.runRates)
print '\tGenLOS:      {0:d}'.format(run.runGenLOS)
print '\tCellfinder:  {0:d}'.format(run.runCellfinder)
print '\tIDCells:     {0:d}'.format(run.runIdcells)
print '\tLos7:        {0:d}'.format(run.runLos7)
print '\tSpecsynth:   {0:d}'.format(run.runSpecsynth)
print '\tSysanal:     {0:d}'.format(run.runSysanal)
print '\tCullabs:     {0:d}'.format(run.runCullabs)
print '\tLocateCells: {0:d}'.format(run.runLocateCells)
print '\tPlotting:    {0:d}'.format(run.runPlotting)

# Test the summary location
sumFile = '{0:s}/rotmat_a{1:s}.txt'.format(os.getcwd(), run.expn)

# Generate gal_props.dat file
print '\n\nGenerating gal_props.dat...'
#fi.setup_galprops( galID, expn, requiredLoc, summaryLoc )
fi.setup_galprops( run, requiredLoc, summaryLoc )

# Generate galaxy.props file, needed for analysis codes
print '\n\nGenerating galaxy.props...'
#fi.setup_galaxy_props(sumFile, galID, expn, incline)
fi.setup_galaxy_props(run, sumFile)

##### 
#  
#  Run rates
#
#####
if run.runRates==1:
    print '\nRates:'
    print '\t Setting up rates.inp...'
    rc.setup_rates_control( gasfile, run.expn, ions, requiredLoc)
    print '\t Setting up rates.outfiles...'
    rc.setup_rates_outputs(run, ions, codeLoc, requiredLoc) 
    print '\t Setting up rates data location...'
    rc.setup_rates_data(codeLoc)
    print '\t Running rates...'
    rc.run_rates(codeLoc)
    print '\t Converting to hdf5...'
    hdf.gasbox_to_hdf5(codeLoc, ions, run)
else:
    print 'Skipping rates'



##### 
#  
#  Set up inclination directory
#
#####
mainLoc = os.getcwd()
#incLoc = fi.setup_inclination_dir(incline, ions, runRates, galID, expn)
incLoc = fi.setup_inclination_dir(run, ions)
os.chdir(incLoc)



##### 
#  
#  Generate the lines of sight
#
#####
if run.runGenLOS==1:
    print '\nGenerating lines of sight...'
    rf.genLOS( codeLoc, mainLoc, run )
else:
    print 'Skipping genLOS...'

##### 
#  
#  Run lines of sight through box with cellfinder
#
#####
if run.runCellfinder==1:
    print '\nRunning LOS through box...'
    rf.runCellfinder(codeLoc, run.ncores)
else: 
    print 'Skipping cellfinder...'


##### 
#  
#  Identify cells in the ion boxes
#
#####
if run.runIdcells==1:
    print '\nIdentifying probed cell properites...'
    ic.idcells(run, ions, codeLoc)
else:
    print 'Skipping IDcells...'

# Setup Mockspec 
fi.setup_mockspec(ions, run, requiredLoc)



# Setup the ion directory
#ionloc = fi.setup_ion_dir(ion, galID, expn, codeLoc) 
for ion in ions:
    ionloc = fi.setup_ion_dir(ion, run, codeLoc) 

# Start looping over ions
jl.Parallel(n_jobs=run.ncores,verbose=5)(
        jl.delayed(ionLoop)(run,ion) for ion in ions)
    
##### 
#  
#  Rename files
#  Add inclination angle to filenames
#
#####
if run.runCullabs==1 or run.runLocateCells==1:
    for ion in ions:
        print '\n\t Renaming files...'
        fi.rename(run,ion)
else:
    print 'What was the point...'



# Generate summary files
if run.runSummaries==1:
    print '\n\nGenearting summary files'
    hdf.genSummaries(galID, expn, incline, ions, nlos)    


# Create basic plots
if run.runPlotting==1:
    print '\n\nGenerating plots'
    ac.make_plots(ions)


